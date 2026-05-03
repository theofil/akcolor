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

NN_DIR   = Path(__file__).resolve().parent
FIGS     = NN_DIR / 'figs'
MODELS   = NN_DIR / 'models'
DATA_DIR = NN_DIR.parents[1]

BKG_FILE = DATA_DIR / 'QCDHtoInv.root'
SIG_FILE = DATA_DIR / 'VBFHtoInv.root'

arch         = json.load(open(NN_DIR / 'architecture.json'))
FEATURE_SETS = arch['feature_sets']
MJJ_CUT      = arch['mjj_cut']

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# region boundaries
NN_LOW_MAX    = 0.65
NN_HIGH_MIN   = 0.65
SPVA_LOW_MAX  = 0.9
SPVA_HIGH_MIN = 2.0


def load_file(path):
    with uproot.open(str(path)) as f:
        flat   = f['events'].arrays(['mjj', 'dYjj', 'dPhijj', 'ptjj'], library='np')
        jagged = f['events'].arrays(['jetSPVA', 'jetPVM'], library='ak')
    spva  = ak.pad_none(jagged['jetSPVA'], 2, axis=1)
    pvm   = ak.pad_none(jagged['jetPVM'],  2, axis=1)
    spva0 = np.abs(ak.to_numpy(ak.fill_none(spva[:, 0], np.nan)))
    spva1 = np.abs(ak.to_numpy(ak.fill_none(spva[:, 1], np.nan)))
    pvm0  = ak.to_numpy(ak.fill_none(pvm[:, 0], np.nan))
    pvm1  = ak.to_numpy(ak.fill_none(pvm[:, 1], np.nan))
    X = np.column_stack([
        flat['mjj'], np.abs(flat['dYjj']), np.abs(flat['dPhijj']), flat['ptjj'],
        spva0, spva1, pvm0, pvm1,
    ])
    mask = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X).any(axis=1)
    return X[mask].astype(np.float32)


def get_scores(key, X):
    cfg  = FEATURE_SETS[key]
    cols = cfg['cols']
    norm = np.load(NN_DIR / f'norm_{key}.npz')
    mu, std = norm['mu'], norm['std']
    Xk = (X[:, cols] - mu) / std
    Xt = torch.from_numpy(Xk.astype(np.float32)).to(device)
    model = MLP(len(cols), arch['hidden_layers'], arch['dropout'], arch['batch_norm']).to(device)
    model.load_state_dict(torch.load(MODELS / f'model_{key}.pt', map_location=device))
    model.eval()
    with torch.no_grad():
        return model(Xt).cpu().numpy()


def norm_hist(values, bins):
    bw = bins[1] - bins[0]
    counts, _ = np.histogram(values, bins=bins)
    norm = counts.sum() * bw
    if norm > 0:
        return counts / norm, np.sqrt(counts) / norm
    return np.zeros_like(counts, float), np.zeros_like(counts, float)


def save_fig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {name}')


print('Loading data...')
X_bkg = load_file(BKG_FILE)
X_sig = load_file(SIG_FILE)

print('Running kinematic NN inference...')
sc_bkg = get_scores('a', X_bkg)
sc_sig = get_scores('a', X_sig)

mjj_bkg   = X_bkg[:, 0]
spva0_bkg = X_bkg[:, 4]
mjj_sig   = X_sig[:, 0]
spva0_sig = X_sig[:, 4]


def region_masks(sc, spva0):
    low_nn   = sc     < NN_LOW_MAX
    high_nn  = sc     > NN_HIGH_MIN
    low_spva = spva0  < SPVA_LOW_MAX
    high_spva= spva0  > SPVA_HIGH_MIN
    return {
        'A': low_nn  & high_spva,
        'B': high_nn & high_spva,
        'C': low_nn  & low_spva,
        'D': high_nn & low_spva,
    }


masks_bkg = region_masks(sc_bkg, spva0_bkg)
masks_sig = region_masks(sc_sig, spva0_sig)

region_style = {
    'A': {'color': 'steelblue',  'label': r'A: low NN, high $|\theta_s|$'},
    'B': {'color': 'crimson',    'label': r'B: high NN, high $|\theta_s|$'},
    'C': {'color': 'darkgreen',  'label': r'C: low NN, low $|\theta_s|$'},
    'D': {'color': 'darkorange', 'label': r'D: high NN, low $|\theta_s|$'},
}

print('\nEvent counts per region:')
print(f'  {"Region":<6}  {"QCD H":>8}  {"VBF H":>8}')
for r in ('A', 'B', 'C', 'D'):
    print(f'  {r:<6}  {masks_bkg[r].sum():>8}  {masks_sig[r].sum():>8}')

bins        = np.linspace(400, 1000, 36)
bin_centers = (bins[:-1] + bins[1:]) / 2

for sample, mjj, masks, tag in [
    ('QCD H', mjj_bkg, masks_bkg, 'QCD'),
    ('VBF H', mjj_sig, masks_sig, 'VBF'),
]:
    fig, ax = plt.subplots(figsize=(6, 6))

    for r, style in region_style.items():
        h, herr = norm_hist(mjj[masks[r]], bins)
        ax.hist(bin_centers, bins=bins, weights=h, histtype='step',
                linewidth=3, color=style['color'], label=style['label'])
        ax.errorbar(bin_centers, h, yerr=herr, fmt='none',
                    color=style['color'], capsize=3, capthick=1, elinewidth=1)

    ax.set_xlabel(r'$m_{jj}$ [GeV]', fontsize=20)
    ax.set_ylabel('Normalized density [1/GeV]', fontsize=20)
    ax.set_title(sample, fontsize=22)
    ax.set_xlim(400, 1000)
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.legend(loc='upper right', frameon=False, fontsize=14)

    plt.tight_layout()
    save_fig(fig, f'ABCD_spva_{tag}')

viewer = Path.home() / 'qplots' / 'viewer.py'
print(f'\nDone. Figures saved to {FIGS.resolve()}')
os.execv(sys.executable, [sys.executable, str(viewer), str(FIGS.resolve())])
