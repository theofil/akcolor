import sys
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
DATA_DIR = NN_DIR.parents[1]

BKG_FILE = DATA_DIR / 'QCDWtoInv.root'
SIG_FILE = DATA_DIR / 'VBFWtoInv.root'

arch = json.load(open(NN_DIR / 'architecture.json'))
FEATURE_SETS = arch['feature_sets']
MJJ_CUT      = arch['mjj_cut']

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# ── data loading (mirrors inference_ROC_Z.py) ──────────────────────────────────
# column layout:
#   0-3:   mjj, dYjj, dPhijj, ptjj
#   4-5:   PVM[0], PVM[1]
#   6-13:  jetPt[0..1], jetEta[0..1], jetPhi[0..1], jetM[0..1]
#   14-15: c21, c12
#   16-17: PVA[0], PVA[1]
def load_file(path):
    with uproot.open(str(path)) as f:
        flat   = f['events'].arrays(['mjj', 'dYjj', 'dPhijj', 'ptjj', 'c21', 'c12'], library='np')
        jagged = f['events'].arrays(['jetSPVA', 'jetPVM', 'jetPt', 'jetEta', 'jetPhi', 'jetM', 'jetPVA'], library='ak')

    def pad2(arr): return ak.pad_none(arr, 2, axis=1)
    def col(arr, i): return ak.to_numpy(ak.fill_none(pad2(arr)[:, i], np.nan))

    spva0_abs = np.abs(col(jagged['jetSPVA'], 0))

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
    return X[mask].astype(np.float32), spva0_abs[mask].astype(np.float32)


print('Loading W data...')
X_bkg, spva0_bkg = load_file(BKG_FILE)
X_sig, spva0_sig = load_file(SIG_FILE)

X     = np.concatenate([X_bkg,     X_sig],     axis=0)
spva0 = np.concatenate([spva0_bkg, spva0_sig], axis=0)
y = np.concatenate([np.zeros(len(X_bkg), dtype=np.float32),
                    np.ones (len(X_sig), dtype=np.float32)])

# no W training split — run on full dataset
X_inf     = X
spva0_inf = spva0
y_inf     = y
print(f'Full W dataset: {len(y_inf)} events  '
      f'({(y_inf==1).sum()} sig, {(y_inf==0).sum()} bkg)')


# ── per-model inference (Z-trained models applied to W) ────────────────────────
def run_inference(key):
    cfg  = FEATURE_SETS[key]
    cols = cfg['cols']

    norm = np.load(NN_DIR / f'norm_Z_{key}.npz')
    mu, std = norm['mu'], norm['std']

    Xk = (X_inf[:, cols] - mu) / std
    Xt = torch.from_numpy(Xk.astype(np.float32)).to(device)

    model = MLP(len(cols), arch['hidden_layers'],
                arch['dropout'], arch['batch_norm']).to(device)
    model.load_state_dict(torch.load(MODELS / f'model_Z_{key}.pt',
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

for key in ('a', 'c', 'd', 'e', 'f'):
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
               ('spva0', -spva0_inf)]:
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
colors = {'a': 'steelblue', 'c': 'darkgreen', 'd': 'darkorchid',
          'e': 'sienna', 'f': 'teal', 'mjj': 'darkorange', 'spva0': 'purple'}
abbrev = {'kinematic': 'kin', 'pull': 'pull', 'combined': 'comb', 'raw': 'raw', 'wild': 'wild',
          'pull_new': 'pull'}
labels = {key: f"NN {abbrev[FEATURE_SETS[key]['name']]}  AUC = {aucs[key]:.2f}"
          for key in ('a', 'c', 'd', 'e', 'f')}
labels['mjj']   = f"$m_{{jj}}$  AUC = {aucs['mjj']:.2f}"
labels['spva0'] = f"$|\\theta_s[j_1]|$  AUC = {aucs['spva0']:.2f}"

fig, ax = plt.subplots(figsize=(6, 5))
for key in ('a', 'c', 'd', 'e', 'f'):
    ax.plot(fprs[key], tprs[key], lw=2, color=colors[key], label=labels[key])
for key in ('mjj', 'spva0'):
    ax.plot(fprs[key], tprs[key], lw=2, color=colors[key],
            linestyle='--', label=labels[key])
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4)
ax.set_xlabel('False Positive Rate', fontsize=20)
ax.set_ylabel('True Positive Rate', fontsize=20)
ax.set_title('VBF W vs QCD W', fontsize=22)
ax.legend(fontsize=14, loc='lower right')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.tick_params(axis='both', labelsize=16)
ax.grid(True, linestyle=':', alpha=0.4)
plt.tight_layout()
save_fig(fig, 'comparison')

print(f'\nDone. Figures saved to {FIGS}')
