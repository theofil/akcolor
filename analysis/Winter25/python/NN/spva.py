import sys
import uproot
import awkward as ak
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from style import apply_style, FIG_SIZE

NN_DIR   = Path(__file__).resolve().parent
FIGS     = NN_DIR / 'paper' / 'figs'
SCRIPT   = Path(__file__).stem
DATA_DIR = NN_DIR.parents[1]

FIGS.mkdir(exist_ok=True)

samples = [
    {'name': 'QCDHtoInv', 'label': 'QCD h', 'file': str(DATA_DIR / 'QCDHtoInv.root'), 'group': 'QCD'},
    {'name': 'QCDZtoInv', 'label': 'QCD Z', 'file': str(DATA_DIR / 'QCDZtoInv.root'), 'group': 'QCD'},
    {'name': 'QCDWtoInv', 'label': 'QCD W', 'file': str(DATA_DIR / 'QCDWtoInv.root'), 'group': 'QCD'},
    {'name': 'VBFHtoInv', 'label': 'VBF h', 'file': str(DATA_DIR / 'VBFHtoInv.root'), 'group': 'VBF'},
    {'name': 'VBFZtoInv', 'label': 'VBF Z', 'file': str(DATA_DIR / 'VBFZtoInv.root'), 'group': 'VBF'},
    {'name': 'VBFWtoInv', 'label': 'VBF W', 'file': str(DATA_DIR / 'VBFWtoInv.root'), 'group': 'VBF'},
]

mjj_cut     = 400
bins        = np.linspace(0, np.pi, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width   = bins[1] - bins[0]

colors      = {'H': 'black', 'Z': '#1f77b4', 'W': '#d62728'}
line_styles = {'QCD': '-', 'VBF': ':'}


def load_spva(file_path, mjj_min, mjj_max=None):
    with uproot.open(file_path) as f:
        tree    = f['events']
        mjj     = tree['mjj'].array(library='ak')
        jetSPVA = tree['jetSPVA'].array(library='ak')
    mask = mjj > mjj_min
    if mjj_max is not None:
        mask = mask & (mjj < mjj_max)
    spva0 = jetSPVA[mask, 0]
    spva0 = ak.to_numpy(ak.fill_none(spva0, np.nan))
    return np.abs(spva0[~np.isnan(spva0)])


def norm_hist(values, bins):
    counts, _ = np.histogram(values, bins=bins)
    norm = counts.sum() * bin_width
    if norm > 0:
        return counts / norm, np.sqrt(counts) / norm
    return np.zeros_like(counts, dtype=float), np.zeros_like(counts, dtype=float)


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{SCRIPT}_{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved figs/{name}.pdf + .png')


print(f'Plotting per-sample SPVA (mjj > {mjj_cut} GeV)...')
fig, ax = plt.subplots(figsize=FIG_SIZE)

for s in samples:
    label = s['label']
    key   = 'H' if 'h' in label.lower() else 'Z' if 'Z' in label else 'W'
    spva  = load_spva(s['file'], mjj_cut)
    hist, hist_err = norm_hist(spva, bins)
    ax.hist(bin_centers, bins=bins, weights=hist, histtype='step',
            linewidth=3, color=colors[key], linestyle=line_styles[s['group']], label=label)
    ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                color=colors[key], capsize=3, capthick=1, elinewidth=1)

apply_style(ax,
            xlabel=r'$|\theta_s|$ [rad]',
            ylabel='Normalized density [1/rad]',
            title='',
            xlim=(0, np.pi), ylim=(0.25, 0.47),
            legend_loc='upper right')
plt.tight_layout()
save(fig, f'leading_jet_mjj_gt_{mjj_cut}')

print(f'\nAll figures saved to {FIGS.resolve()}')
