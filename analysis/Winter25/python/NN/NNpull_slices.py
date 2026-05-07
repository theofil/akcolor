import sys
import json
import numpy as np
import uproot
import awkward as ak
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from model import MLP
from style import apply_style, FIG_SIZE

NN_DIR   = Path(__file__).resolve().parent
FIGS     = NN_DIR / 'paper' / 'figs'
SCRIPT   = Path(__file__).stem
MODELS   = NN_DIR / 'models'
DATA_DIR = NN_DIR.parents[1]

FIGS.mkdir(exist_ok=True)

arch         = json.load(open(NN_DIR / 'architecture.json'))
FEATURE_SETS = arch['feature_sets']
MJJ_CUT      = arch['mjj_cut']

KEY  = 'f'
cfg  = FEATURE_SETS[KEY]
COLS = cfg['cols']

norm     = np.load(NN_DIR / f'norm_{KEY}.npz')
mu, std  = norm['mu'], norm['std']

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

model = MLP(len(COLS), arch['hidden_layers'], arch['dropout'], arch['batch_norm']).to(device)
model.load_state_dict(torch.load(MODELS / f'model_{KEY}.pt',
                                 map_location=device, weights_only=True))
model.eval()

MJJ_SLICES   = [(400, 550), (550, 700), (700, 900), (900, None)]
SLICE_COLORS = ['black', 'steelblue', 'darkgreen', 'crimson']

samples = [
    {'label': 'QCD h', 'file': DATA_DIR / 'QCDHtoInv.root', 'suffix': 'QCD'},
    {'label': 'VBF h', 'file': DATA_DIR / 'VBFHtoInv.root', 'suffix': 'VBF'},
]

bins        = np.linspace(0, 1, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width   = bins[1] - bins[0]


def load_and_score(path, mjj_min, mjj_max):
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

    base_mask = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X).any(axis=1)
    mjj_mask  = flat['mjj'] > mjj_min
    if mjj_max is not None:
        mjj_mask = mjj_mask & (flat['mjj'] < mjj_max)
    mask = base_mask & mjj_mask

    Xk = (X[mask][:, COLS] - mu) / std
    Xt = torch.from_numpy(Xk.astype(np.float32)).to(device)
    with torch.no_grad():
        scores = model(Xt).cpu().numpy().ravel()
    return scores


def norm_hist(values):
    counts, _ = np.histogram(values, bins=bins)
    norm = counts.sum() * bin_width
    if norm > 0:
        return counts / norm, np.sqrt(counts) / norm
    return np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)


def slice_title(lo, hi):
    if hi is None:
        return fr'$m_{{jj}} > {lo}$ GeV'
    return fr'${lo} < m_{{jj}} < {hi}$ GeV'


for s in samples:
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for color, (lo, hi) in zip(SLICE_COLORS, MJJ_SLICES):
        tag = slice_title(lo, hi)
        scores = load_and_score(s['file'], lo, hi)
        h, herr = norm_hist(scores)
        ax.hist(bin_centers, bins=bins, weights=h, histtype='step',
                linewidth=2.5, color=color, label=tag)
        ax.errorbar(bin_centers, h, yerr=herr, fmt='none',
                    color=color, capsize=3, capthick=1, elinewidth=1)

    apply_style(ax,
                xlabel='NN pull score',
                ylabel='Normalized density',
                title=s['label'],
                xlim=(0, 1),
                legend_loc='upper left')

    ax.legend(loc='upper left', frameon=False, fontsize=11, ncol=1)
    plt.tight_layout()
    out = f'{SCRIPT}_mjj_{s["suffix"]}'
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{out}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}.pdf/png to {FIGS}')
