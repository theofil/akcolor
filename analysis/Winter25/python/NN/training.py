import sys
import json
import numpy as np
import uproot
import awkward as ak
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from model import MLP

# ── paths ──────────────────────────────────────────────────────────────────────
NN_DIR  = Path(__file__).resolve().parent
FIGS    = NN_DIR / 'paper' / 'figs'
SCRIPT  = Path(__file__).stem
MODELS  = NN_DIR / 'models'
DATA_DIR = NN_DIR.parents[1]   # Winter25/

BKG_FILE = DATA_DIR / 'QCDHtoInv.root'
SIG_FILE = DATA_DIR / 'VBFHtoInv.root'

arch = json.load(open(NN_DIR / 'architecture.json'))
FEATURE_SETS = arch['feature_sets']
MJJ_CUT      = arch['mjj_cut']

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Device: {device}')

# ── data loading ───────────────────────────────────────────────────────────────
# column layout: [mjj, |dYjj|, |dPhijj|, ptjj, |SPVA[0]|, |SPVA[1]|, PVM[0], PVM[1]]

def load_file(path, label):
    print(f'Loading {path.name}...')
    with uproot.open(str(path)) as f:
        flat   = f['events'].arrays(['mjj', 'dYjj', 'dPhijj', 'ptjj'], library='np')
        jagged = f['events'].arrays(['jetSPVA', 'jetPVM'], library='ak')

    spva = ak.pad_none(jagged['jetSPVA'], 2, axis=1)
    pvm  = ak.pad_none(jagged['jetPVM'],  2, axis=1)
    spva0 = np.abs(ak.to_numpy(ak.fill_none(spva[:, 0], np.nan)))
    spva1 = np.abs(ak.to_numpy(ak.fill_none(spva[:, 1], np.nan)))
    pvm0  = ak.to_numpy(ak.fill_none(pvm[:, 0], np.nan))
    pvm1  = ak.to_numpy(ak.fill_none(pvm[:, 1], np.nan))

    X = np.column_stack([
        flat['mjj'],
        np.abs(flat['dYjj']),
        np.abs(flat['dPhijj']),
        flat['ptjj'],
        spva0, spva1, pvm0, pvm1,
    ])

    # mjj cut + drop any NaN row
    mask = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X).any(axis=1)
    X = X[mask].astype(np.float32)
    print(f'  {X.shape[0]} events after mjj > {MJJ_CUT} and NaN removal')
    return X

X_bkg = load_file(BKG_FILE, 'bkg')
X_sig = load_file(SIG_FILE, 'sig')

X = np.concatenate([X_bkg, X_sig], axis=0)
y = np.concatenate([np.zeros(len(X_bkg), dtype=np.float32),
                    np.ones (len(X_sig), dtype=np.float32)])

# ── train / val / inference split (60 / 20 / 20) ──────────────────────────────
rng = np.random.default_rng(42)
idx = rng.permutation(len(X))
n_train = int(0.6 * len(X))
n_val   = int(0.2 * len(X))

idx_train = idx[:n_train]
idx_val   = idx[n_train:n_train + n_val]
idx_inf   = idx[n_train + n_val:]

np.savez(NN_DIR / 'split_indices.npz',
         idx_train=idx_train, idx_val=idx_val, idx_inf=idx_inf)
print(f'\nSplit — train: {len(idx_train)}  val: {len(idx_val)}  inference: {len(idx_inf)}')


# ── helpers ────────────────────────────────────────────────────────────────────
def normalize(X_tr, X_va, X_te):
    mu  = X_tr.mean(0)
    std = X_tr.std(0) + 1e-8
    return ((X_tr-mu)/std, (X_va-mu)/std, (X_te-mu)/std, mu, std)


def save_fig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{SCRIPT}_{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {name}')


# ── training loop ──────────────────────────────────────────────────────────────
def train_model(key):
    cfg  = FEATURE_SETS[key]
    cols = cfg['cols']
    name = cfg['name']
    print(f'\n{"─"*60}')
    print(f'Model {key}  ({name})  features: {cfg["labels"]}')
    print(f'{"─"*60}')

    Xk = X[:, cols]
    Xk_tr, Xk_va, Xk_te, mu, std = normalize(
        Xk[idx_train], Xk[idx_val], Xk[idx_inf])
    y_tr = y[idx_train]
    y_va = y[idx_val]

    np.savez(NN_DIR / f'norm_{key}.npz', mu=mu, std=std)

    loader_tr = DataLoader(
        TensorDataset(torch.from_numpy(Xk_tr), torch.from_numpy(y_tr)),
        batch_size=arch['batch_size'], shuffle=True, drop_last=True)
    loader_va = DataLoader(
        TensorDataset(torch.from_numpy(Xk_va), torch.from_numpy(y_va)),
        batch_size=arch['batch_size'] * 4)

    model = MLP(len(cols), arch['hidden_layers'],
                arch['dropout'], arch['batch_norm']).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=arch['learning_rate'])
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=5, factor=0.5, verbose=False)
    criterion = nn.BCELoss()

    hist = {'train_loss': [], 'val_loss': [], 'val_auc': []}
    best_val  = np.inf
    no_improve = 0

    for epoch in range(arch['epochs']):
        # ── train ──
        model.train()
        t_loss = 0.0
        for xb, yb in loader_tr:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(xb), yb)
            loss.backward()
            optimizer.step()
            t_loss += loss.item() * len(xb)
        t_loss /= len(loader_tr.dataset)

        # ── validate ──
        model.eval()
        v_loss, preds = 0.0, []
        with torch.no_grad():
            for xb, yb in loader_va:
                xb, yb = xb.to(device), yb.to(device)
                p = model(xb)
                v_loss += criterion(p, yb).item() * len(xb)
                preds.extend(p.cpu().numpy())
        v_loss /= len(loader_va.dataset)
        v_auc = roc_auc_score(y_va, preds)

        scheduler.step(v_loss)
        hist['train_loss'].append(t_loss)
        hist['val_loss'].append(v_loss)
        hist['val_auc'].append(v_auc)

        if epoch % 10 == 0 or epoch < 3:
            lr_now = optimizer.param_groups[0]['lr']
            print(f'  ep {epoch:3d}  train={t_loss:.4f}  val={v_loss:.4f}'
                  f'  AUC={v_auc:.4f}  lr={lr_now:.2e}')

        if v_loss < best_val:
            best_val   = v_loss
            no_improve = 0
            torch.save(model.state_dict(), MODELS / f'model_{key}.pt')
        else:
            no_improve += 1
            if no_improve >= arch['patience']:
                print(f'  early stop at epoch {epoch}  (best val={best_val:.4f})')
                break

    # ── training history plots ─────────────────────────────────────────────────
    epochs = range(len(hist['train_loss']))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].plot(epochs, hist['train_loss'], label='train', color='steelblue', lw=2)
    axes[0].plot(epochs, hist['val_loss'],   label='val',   color='crimson',   lw=2)
    axes[0].set_xlabel('Epoch', fontsize=20)
    axes[0].set_ylabel('BCE Loss', fontsize=20)
    axes[0].set_title(f'Model {key} ({name}) — Loss', fontsize=22)
    axes[0].legend(fontsize=14)
    axes[0].tick_params(axis='both', labelsize=16)
    axes[0].grid(True, linestyle=':', alpha=0.4)

    axes[1].plot(epochs, hist['val_auc'], color='green', lw=2)
    axes[1].set_xlabel('Epoch', fontsize=20)
    axes[1].set_ylabel('AUC', fontsize=20)
    axes[1].set_title(f'Model {key} — Validation AUC', fontsize=22)
    axes[1].set_ylim(0.5, 1.0)
    axes[1].tick_params(axis='both', labelsize=16)
    axes[1].grid(True, linestyle=':', alpha=0.4)

    # loss ratio (val/train) — should stay near 1; drifting up = overtraining
    ratio = np.array(hist['val_loss']) / np.array(hist['train_loss'])
    axes[2].plot(epochs, ratio, color='darkorange', lw=2)
    axes[2].axhline(1.0, color='black', linestyle='--', lw=1, alpha=0.5)
    axes[2].set_xlabel('Epoch', fontsize=20)
    axes[2].set_ylabel('val loss / train loss', fontsize=20)
    axes[2].set_title(f'Model {key} — Overtraining monitor', fontsize=22)
    axes[2].tick_params(axis='both', labelsize=16)
    axes[2].grid(True, linestyle=':', alpha=0.4)

    plt.tight_layout()
    save_fig(fig, f'history_{key}')

    print(f'  best val loss={best_val:.4f}  final AUC={hist["val_auc"][-1]:.4f}')
    return hist


for key in ('a', 'b', 'c'):
    train_model(key)

print('\nAll models trained. Run inference_ROC.py to compare.')
