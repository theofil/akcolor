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

DISCO_DIR    = Path(__file__).resolve().parent
NN_DIR       = DISCO_DIR.parent
DATA_DIR     = NN_DIR.parents[1]
MODELS_DISCO = DISCO_DIR / 'models'
MODELS_DISCO.mkdir(exist_ok=True)

sys.path.insert(0, str(NN_DIR))
from model import MLP

arch         = json.load(open(NN_DIR / 'architecture.json'))
MJJ_CUT      = arch['mjj_cut']
COLS_D       = [6, 7, 8, 9, 10, 11, 12, 13]
LAMBDA       = 10.0
N_DISCO_MIN  = 2

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Device: {device}')


def load_file(path):
    print(f'Loading {path.name}...')
    with uproot.open(str(path)) as f:
        flat   = f['events'].arrays(['mjj', 'dYjj', 'dPhijj', 'ptjj', 'c21', 'c12'], library='np')
        jagged = f['events'].arrays(['jetPVM', 'jetPt', 'jetEta', 'jetPhi', 'jetM', 'jetPVA'], library='ak')

    def pad2(arr): return ak.pad_none(arr, 2, axis=1)
    def col(arr, i): return ak.to_numpy(ak.fill_none(pad2(arr)[:, i], np.nan))

    X = np.column_stack([
        flat['mjj'], flat['dYjj'], flat['dPhijj'], flat['ptjj'],
        col(jagged['jetPVM'], 0), col(jagged['jetPVM'], 1),
        col(jagged['jetPt'],  0), col(jagged['jetPt'],  1),
        col(jagged['jetEta'], 0), col(jagged['jetEta'], 1),
        col(jagged['jetPhi'], 0), col(jagged['jetPhi'], 1),
        col(jagged['jetM'],   0), col(jagged['jetM'],   1),
        flat['c21'], flat['c12'],
        col(jagged['jetPVA'], 0), col(jagged['jetPVA'], 1),
    ])
    mask = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X).any(axis=1)
    X = X[mask].astype(np.float32)
    print(f'  {X.shape[0]} events after mjj > {MJJ_CUT} and NaN removal')
    return X


print('Loading data...')
X_bkg = load_file(DATA_DIR / 'QCDHtoInv.root')
X_sig = load_file(DATA_DIR / 'VBFHtoInv.root')
X = np.concatenate([X_bkg, X_sig], axis=0)
y = np.concatenate([np.zeros(len(X_bkg), dtype=np.float32),
                    np.ones (len(X_sig), dtype=np.float32)])

split     = np.load(NN_DIR / 'split_indices.npz')
idx_train = split['idx_train']
idx_val   = split['idx_val']

norm     = np.load(NN_DIR / 'norm_d.npz')
mu, std  = norm['mu'], norm['std']

Xk_tr  = (X[idx_train][:, COLS_D] - mu) / std
Xk_va  = (X[idx_val  ][:, COLS_D] - mu) / std
y_tr   = y[idx_train]
y_va   = y[idx_val]
mjj_tr = X[idx_train][:, 0]

print(f'Train: {len(idx_train)}  Val: {len(idx_val)}')

loader_tr = DataLoader(
    TensorDataset(
        torch.from_numpy(Xk_tr),
        torch.from_numpy(y_tr),
        torch.from_numpy(mjj_tr),
    ),
    batch_size=arch['batch_size'], shuffle=True, drop_last=True)

loader_va = DataLoader(
    TensorDataset(torch.from_numpy(Xk_va), torch.from_numpy(y_va)),
    batch_size=arch['batch_size'] * 4)


def dcor2(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    A = (x.unsqueeze(0) - x.unsqueeze(1)).abs()
    B = (y.unsqueeze(0) - y.unsqueeze(1)).abs()
    A = A - A.mean(1, keepdim=True) - A.mean(0, keepdim=True) + A.mean()
    B = B - B.mean(1, keepdim=True) - B.mean(0, keepdim=True) + B.mean()
    dCov  = (A * B).mean()
    dVarX = (A * A).mean()
    dVarY = (B * B).mean()
    denom = (dVarX * dVarY).clamp(min=1e-20).sqrt()
    dc    = dCov / denom
    return torch.where((dVarX < 1e-10) | (dVarY < 1e-10),
                       torch.zeros_like(dc), dc)


model     = MLP(8, arch['hidden_layers'], arch['dropout'], arch['batch_norm']).to(device)
model.load_state_dict(torch.load(NN_DIR / 'models' / 'model_d.pt',
                                 map_location=device, weights_only=True))
print('Warm-started from model_d.pt')

optimizer = torch.optim.Adam(model.parameters(), lr=arch['learning_rate'])
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, patience=5, factor=0.5, verbose=False)
criterion = nn.BCELoss()

hist = {'train_bce': [], 'train_disco': [], 'train_combined': [],
        'val_bce': [], 'val_auc': []}
best_val_bce = np.inf
no_improve   = 0

for epoch in range(arch['epochs']):
    model.train()
    t_bce = t_disco = 0.0

    for xb, yb, mjj_b in loader_tr:
        xb, yb, mjj_b = xb.to(device), yb.to(device), mjj_b.to(device)
        optimizer.zero_grad()
        scores   = model(xb)
        bce_loss = criterion(scores, yb)

        qcd_mask = (yb == 0)
        n_qcd    = int(qcd_mask.sum())
        if n_qcd >= N_DISCO_MIN:
            disco_val  = LAMBDA * dcor2(scores[qcd_mask], mjj_b[qcd_mask].detach())
        else:
            disco_val  = torch.tensor(0.0, device=device)

        (bce_loss + disco_val).backward()
        optimizer.step()

        bs       = len(xb)
        t_bce   += bce_loss.item()  * bs
        t_disco += disco_val.item() * bs

    t_bce   /= len(loader_tr.dataset)
    t_disco /= len(loader_tr.dataset)

    model.eval()
    v_bce = 0.0
    preds = []
    with torch.no_grad():
        for xb, yb in loader_va:
            xb, yb = xb.to(device), yb.to(device)
            p = model(xb)
            v_bce += criterion(p, yb).item() * len(xb)
            preds.extend(p.cpu().numpy())
    v_bce /= len(loader_va.dataset)
    v_auc  = roc_auc_score(y_va, preds)

    scheduler.step(v_bce)

    hist['train_bce'].append(t_bce)
    hist['train_disco'].append(t_disco)
    hist['train_combined'].append(t_bce + t_disco)
    hist['val_bce'].append(v_bce)
    hist['val_auc'].append(v_auc)

    if epoch % 10 == 0 or epoch < 3:
        lr_now = optimizer.param_groups[0]['lr']
        print(f'  ep {epoch:3d}  bce={t_bce:.4f}  disco={t_disco:.4f}'
              f'  val_bce={v_bce:.4f}  AUC={v_auc:.4f}  lr={lr_now:.2e}')

    if v_bce < best_val_bce:
        best_val_bce = v_bce
        no_improve   = 0
        torch.save(model.state_dict(), MODELS_DISCO / 'model_d_disco.pt')
    else:
        no_improve += 1
        if no_improve >= arch['patience']:
            print(f'  early stop at epoch {epoch}  (best val_bce={best_val_bce:.4f})')
            break

print(f'\nDone. Best val_bce={best_val_bce:.4f}')

# ── training history plot ──────────────────────────────────────────────────────
epochs = range(len(hist['train_bce']))

fig, axes = plt.subplots(1, 4, figsize=(20, 4))

axes[0].plot(epochs, hist['train_bce'], label='train', color='steelblue', lw=2)
axes[0].plot(epochs, hist['val_bce'],   label='val',   color='crimson',   lw=2)
axes[0].set_xlabel('Epoch', fontsize=16)
axes[0].set_ylabel('BCE Loss', fontsize=16)
axes[0].set_title('BCE Loss', fontsize=18)
axes[0].legend(fontsize=12)
axes[0].grid(True, linestyle=':', alpha=0.4)
axes[0].tick_params(labelsize=13)

axes[1].plot(epochs, hist['train_disco'], color='darkorange', lw=2)
axes[1].set_xlabel('Epoch', fontsize=16)
axes[1].set_ylabel(fr'$\lambda \cdot dCor^2$', fontsize=16)
axes[1].set_title('DisCo penalty', fontsize=18)
axes[1].grid(True, linestyle=':', alpha=0.4)
axes[1].tick_params(labelsize=13)

axes[2].plot(epochs, hist['train_combined'], color='purple', lw=2, label='train combined')
axes[2].plot(epochs, hist['val_bce'],        color='crimson', lw=2, linestyle='--', label='val BCE')
axes[2].set_xlabel('Epoch', fontsize=16)
axes[2].set_ylabel('Loss', fontsize=16)
axes[2].set_title('Combined loss', fontsize=18)
axes[2].legend(fontsize=12)
axes[2].grid(True, linestyle=':', alpha=0.4)
axes[2].tick_params(labelsize=13)

axes[3].plot(epochs, hist['val_auc'], color='green', lw=2)
axes[3].set_xlabel('Epoch', fontsize=16)
axes[3].set_ylabel('AUC', fontsize=16)
axes[3].set_title('Validation AUC', fontsize=18)
axes[3].set_ylim(0.5, 1.0)
axes[3].grid(True, linestyle=':', alpha=0.4)
axes[3].tick_params(labelsize=13)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(DISCO_DIR / f'training_history_disco.{ext}', bbox_inches='tight')
plt.close(fig)
print(f'Saved training_history_disco.pdf/png to {DISCO_DIR}')
