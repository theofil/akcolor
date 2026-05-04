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

BKG_FILE = DATA_DIR / 'QCDHtoInv.root'
SIG_FILE = DATA_DIR / 'VBFHtoInv.root'

arch         = json.load(open(NN_DIR / 'architecture.json'))
FEATURE_SETS = arch['feature_sets']
MJJ_CUT      = arch['mjj_cut']

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# region boundaries (same as ABCD_spva.py / ABCD_yields.py)
NN_LOW_MAX    = 0.65
NN_HIGH_MIN   = 0.65
SPVA_LOW_MAX  = 0.9
SPVA_HIGH_MIN = 2.0


def load_file(path):
    with uproot.open(str(path)) as f:
        flat   = f['events'].arrays(['mjj', 'dYjj', 'dPhijj', 'ptjj', 'kWeight'], library='np')
        jagged = f['events'].arrays(['jetSPVA', 'jetPVM'], library='ak')
    spva  = ak.pad_none(jagged['jetSPVA'], 2, axis=1)
    pvm   = ak.pad_none(jagged['jetPVM'],  2, axis=1)
    spva0 = np.abs(ak.to_numpy(ak.fill_none(spva[:, 0], np.nan)))
    spva1 = np.abs(ak.to_numpy(ak.fill_none(spva[:, 1], np.nan)))
    pvm0  = ak.to_numpy(ak.fill_none(pvm[:, 0], np.nan))
    pvm1  = ak.to_numpy(ak.fill_none(pvm[:, 1], np.nan))
    X    = np.column_stack([
        flat['mjj'], np.abs(flat['dYjj']), np.abs(flat['dPhijj']), flat['ptjj'],
        spva0, spva1, pvm0, pvm1,
    ])
    w    = flat['kWeight']
    mask = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X).any(axis=1)
    return X[mask].astype(np.float32), w[mask]


def get_scores(key, X):
    cfg  = FEATURE_SETS[key]
    cols = cfg['cols']
    norm = np.load(NN_DIR / f'norm_{key}.npz')
    mu, std = norm['mu'], norm['std']
    Xk  = (X[:, cols] - mu) / std
    Xt  = torch.from_numpy(Xk.astype(np.float32)).to(device)
    model = MLP(len(cols), arch['hidden_layers'], arch['dropout'], arch['batch_norm']).to(device)
    model.load_state_dict(torch.load(MODELS / f'model_{key}.pt',
                                     map_location=device, weights_only=True))
    model.eval()
    with torch.no_grad():
        return model(Xt).cpu().numpy()


def region_masks(sc, spva0):
    return {
        'A': (sc < NN_LOW_MAX)  & (spva0 > SPVA_HIGH_MIN),
        'B': (sc > NN_HIGH_MIN) & (spva0 > SPVA_HIGH_MIN),
        'C': (sc < NN_LOW_MAX)  & (spva0 < SPVA_LOW_MAX),
        'D': (sc > NN_HIGH_MIN) & (spva0 < SPVA_LOW_MAX),
    }


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{SCRIPT}_{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {name}')


print('Loading data...')
X_bkg, w_bkg = load_file(BKG_FILE)
X_sig, w_sig = load_file(SIG_FILE)

print('Running NN inference...')
sc_bkg = get_scores('a', X_bkg)
sc_sig = get_scores('a', X_sig)

masks_bkg = region_masks(sc_bkg, X_bkg[:, 4])
masks_sig = region_masks(sc_sig, X_sig[:, 4])

bins       = np.linspace(MJJ_CUT, 1000, 31)
bc         = (bins[:-1] + bins[1:]) / 2

region_titles = {
    'A': r'A: low NN, high $|\theta_s|$',
    'B': r'B: high NN, high $|\theta_s|$',
    'C': r'C: low NN, low $|\theta_s|$',
    'D': r'D: high NN, low $|\theta_s|$',
}

samples = [
    {'label': 'QCD $h$', 'mjj': X_bkg[:, 0], 'w': w_bkg, 'masks': masks_bkg,
     'color': 'black', 'linestyle': '-'},
    {'label': 'VBF $h$', 'mjj': X_sig[:, 0], 'w': w_sig, 'masks': masks_sig,
     'color': 'black', 'linestyle': ':'},
]

print('Plotting mjj per ABCD region...')
fig, axes = plt.subplots(2, 2, figsize=(2 * FIG_SIZE[0], 2 * FIG_SIZE[1]))

for ax, r in zip(axes.flat, ('A', 'B', 'C', 'D')):
    for s in samples:
        mask     = s['masks'][r]
        w        = s['w'][mask]
        mjj      = s['mjj'][mask]
        hist, _  = np.histogram(mjj, bins=bins, weights=w)
        hist_err = np.sqrt(np.histogram(mjj, bins=bins, weights=w**2)[0])
        ax.hist(bc, bins=bins, weights=hist, histtype='step',
                linewidth=3, color=s['color'], linestyle=s['linestyle'],
                label=s['label'])
        ax.errorbar(bc, hist, yerr=hist_err, fmt='none',
                    color=s['color'], capsize=3, capthick=1, elinewidth=1)

    apply_style(ax,
                xlabel=r'$m_{jj}$ [GeV]',
                ylabel='Events / pb',
                title=region_titles[r],
                xlim=(MJJ_CUT, 1000),
                log_y=True,
                legend_loc='upper right')

plt.tight_layout()
save(fig, 'regions')

viewer = Path.home() / 'qplots' / 'viewer.py'
print(f'\nFigures saved to {FIGS.resolve()}')
os.execv(sys.executable, [sys.executable, str(viewer), str(FIGS.resolve())])
