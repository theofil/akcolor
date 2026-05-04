import sys
import os
import json
import numpy as np
import uproot
import awkward as ak
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from model import MLP

# ── paths ──────────────────────────────────────────────────────────────────────
NN_DIR   = Path(__file__).resolve().parent
FIGS     = NN_DIR / 'paper' / 'figs'
SCRIPT   = Path(__file__).stem
MODELS   = NN_DIR / 'models'
DATA_DIR = NN_DIR.parents[1]   # Winter25/

BKG_FILE = DATA_DIR / 'QCDHtoInv.root'
SIG_FILE = DATA_DIR / 'VBFHtoInv.root'

arch = json.load(open(NN_DIR / 'architecture.json'))
FEATURE_SETS = arch['feature_sets']
MJJ_CUT      = arch['mjj_cut']

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# ── data loading (mirrors training.py exactly) ─────────────────────────────────
def load_file(path):
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

    mask = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X).any(axis=1)
    return X[mask].astype(np.float32)


print('Loading data...')
X_bkg = load_file(BKG_FILE)
X_sig = load_file(SIG_FILE)

X = np.concatenate([X_bkg, X_sig], axis=0)
y = np.concatenate([np.zeros(len(X_bkg), dtype=np.float32),
                    np.ones (len(X_sig), dtype=np.float32)])

# ── restore inference split ────────────────────────────────────────────────────
split = np.load(NN_DIR / 'split_indices.npz')
idx_inf = split['idx_inf']

X_inf = X[idx_inf]
y_inf = y[idx_inf]
print(f'Inference set: {len(y_inf)} events  '
      f'({(y_inf==1).sum()} sig, {(y_inf==0).sum()} bkg)')


# ── per-model inference ────────────────────────────────────────────────────────
def run_inference(key):
    cfg  = FEATURE_SETS[key]
    cols = cfg['cols']

    norm = np.load(NN_DIR / f'norm_{key}.npz')
    mu, std = norm['mu'], norm['std']

    Xk = (X_inf[:, cols] - mu) / std
    Xt = torch.from_numpy(Xk.astype(np.float32)).to(device)

    model = MLP(len(cols), arch['hidden_layers'],
                arch['dropout'], arch['batch_norm']).to(device)
    model.load_state_dict(torch.load(MODELS / f'model_{key}.pt',
                                     map_location=device,
                                     weights_only=True))
    model.eval()
    with torch.no_grad():
        scores = model(Xt).cpu().numpy()

    return scores


scores = {}
aucs   = {}
fprs   = {}
tprs   = {}

for key in ('a', 'b', 'c'):
    cfg = FEATURE_SETS[key]
    print(f"Running inference for model {key} ({cfg['name']})...")
    s = run_inference(key)
    scores[key] = s
    aucs[key]   = roc_auc_score(y_inf, s)
    fpr, tpr, _ = roc_curve(y_inf, s)
    fprs[key]   = fpr
    tprs[key]   = tpr
    print(f"  AUC = {aucs[key]:.4f}")

# standalone discriminants
for key, s in [('mjj',   X_inf[:, 0]),
               ('spva0', -X_inf[:, 4])]:    # -|SPVA[0]|: lower pull = more VBF-like
    scores[key] = s
    aucs[key]   = roc_auc_score(y_inf, s)
    fpr, tpr, _ = roc_curve(y_inf, s)
    fprs[key]   = fpr
    tprs[key]   = tpr
    print(f"  standalone {key}  AUC = {aucs[key]:.4f}")


# ── helpers ────────────────────────────────────────────────────────────────────
def save_fig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{SCRIPT}_{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {name}')


# ── ROC comparison plot ────────────────────────────────────────────────────────
colors = {'a': 'steelblue', 'b': 'crimson', 'c': 'darkgreen',
          'mjj': 'darkorange', 'spva0': 'purple'}
abbrev = {'kinematic': 'kin', 'pull': 'pull', 'combined': 'comb'}
labels = {key: f"NN {abbrev[FEATURE_SETS[key]['name']]}  AUC = {aucs[key]:.2f}"
          for key in ('a', 'b', 'c')}
labels['b']     = f"NN pull  AUC = {aucs['b']:.2f}"
labels['mjj']   = f"$m_{{jj}}$  AUC = {aucs['mjj']:.2f}"
labels['spva0'] = f"$|\\theta_s[j_1]|$  AUC = {aucs['spva0']:.2f}"

fig, ax = plt.subplots(figsize=(6, 5))
for key in ('a', 'b', 'c'):
    ax.plot(fprs[key], tprs[key], lw=2, color=colors[key], label=labels[key])
for key in ('mjj', 'spva0'):
    ax.plot(fprs[key], tprs[key], lw=2, color=colors[key],
            linestyle='--', label=labels[key])
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4)
ax.set_xlabel('False Positive Rate', fontsize=20)
ax.set_ylabel('True Positive Rate', fontsize=20)
ax.set_title('VBF h vs QCD h', fontsize=22)
ax.legend(fontsize=14, loc='lower right')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.tick_params(axis='both', labelsize=16)
ax.grid(True, linestyle=':', alpha=0.4)
plt.tight_layout()
save_fig(fig, 'comparison')


# ── score distributions ────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
bins = np.linspace(0, 1, 51)

for ax, key in zip(axes, ('a', 'b', 'c')):
    cfg  = FEATURE_SETS[key]
    s    = scores[key]
    mask_bkg = y_inf == 0
    mask_sig = y_inf == 1

    ax.hist(s[mask_bkg], bins=bins, density=True,
            histtype='step', color='black',   lw=2, label='QCD (bkg)')
    ax.hist(s[mask_sig], bins=bins, density=True,
            histtype='step', color='crimson', lw=2, label='VBF (sig)', linestyle='--')
    ax.set_xlabel('NN score', fontsize=20)
    ax.set_ylabel('Normalised events', fontsize=20)
    ax.set_title(f"Model ({key}) {cfg['name']}  AUC={aucs[key]:.4f}", fontsize=22)
    ax.legend(fontsize=14)
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(True, linestyle=':', alpha=0.4)

plt.tight_layout()
save_fig(fig, 'score_distributions')


# ── launch viewer ──────────────────────────────────────────────────────────────
viewer = Path.home() / 'qplots' / 'viewer.py'
print(f'\nDone. Launching viewer at {FIGS}')
os.execv(sys.executable, [sys.executable, str(viewer), str(FIGS)])
