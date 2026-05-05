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
    {'label': 'QCD h', 'file': str(DATA_DIR / 'QCDHtoInv.root'), 'group': 'QCD'},
    {'label': 'QCD Z', 'file': str(DATA_DIR / 'QCDZtoInv.root'), 'group': 'QCD'},
    {'label': 'QCD W', 'file': str(DATA_DIR / 'QCDWtoInv.root'), 'group': 'QCD'},
    {'label': 'VBF h', 'file': str(DATA_DIR / 'VBFHtoInv.root'), 'group': 'VBF'},
    {'label': 'VBF Z', 'file': str(DATA_DIR / 'VBFZtoInv.root'), 'group': 'VBF'},
    {'label': 'VBF W', 'file': str(DATA_DIR / 'VBFWtoInv.root'), 'group': 'VBF'},
]

mjj_cut     = 400
bins        = np.linspace(0, 0.06, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width   = bins[1] - bins[0]

colors      = {'H': 'black', 'Z': '#1f77b4', 'W': '#d62728'}
line_styles = {'QCD': '-', 'VBF': ':'}


def load_pvm(file_path, mjj_min, mjj_max=None):
    with uproot.open(file_path) as f:
        tree   = f['events']
        mjj    = tree['mjj'].array(library='ak')
        jetPVM = tree['jetPVM'].array(library='ak')
    mask = mjj > mjj_min
    if mjj_max is not None:
        mask = mask & (mjj < mjj_max)
    pvm0 = ak.to_numpy(ak.fill_none(jetPVM[mask, 0], np.nan))
    return pvm0[~np.isnan(pvm0)]


def norm_hist(values):
    counts, _ = np.histogram(values, bins=bins)
    norm = counts.sum() * bin_width
    if norm > 0:
        return counts / norm, np.sqrt(counts) / norm
    return np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{SCRIPT}_{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved figs/{SCRIPT}_{name}.pdf + .png')


# ── Figure 1: per-sample PVM, mjj > 400 GeV ────────────────────────────────
print(f'Plotting per-sample PVM (mjj > {mjj_cut} GeV)...')
fig, ax = plt.subplots(figsize=FIG_SIZE)

for s in samples:
    label = s['label']
    key   = 'H' if 'h' in label.lower() else 'Z' if 'Z' in label else 'W'
    pvm   = load_pvm(s['file'], mjj_cut)
    hist, hist_err = norm_hist(pvm)
    ax.hist(bin_centers, bins=bins, weights=hist, histtype='step',
            linewidth=3, color=colors[key], linestyle=line_styles[s['group']], label=label)
    ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                color=colors[key], capsize=3, capthick=1, elinewidth=1)

apply_style(ax,
            xlabel=r'$|\vec{t}\,|$',
            ylabel='Normalized density',
            title='',
            xlim=(0, 0.02),
            log_y=True,
            legend_loc='upper right')
plt.tight_layout()
save(fig, f'leading_jet_mjj_gt_{mjj_cut}')


# ── Figure 2: PVM in mjj slices for QCD h vs VBF h ─────────────────────────
MJJ_SLICES   = [(400, 550), (550, 700), (700, 900), (900, None)]
SLICE_COLORS = ['black', 'steelblue', 'darkgreen', 'crimson']

samples_slices = [
    {'label': 'QCD h', 'file': DATA_DIR / 'QCDHtoInv.root', 'ls': '-'},
    {'label': 'VBF h', 'file': DATA_DIR / 'VBFHtoInv.root', 'ls': '--'},
]


def slice_title(lo, hi):
    if hi is None:
        return fr'$m_{{jj}} > {lo}$ GeV'
    return fr'${lo} < m_{{jj}} < {hi}$ GeV'


print('Plotting PVM in mjj slices...')
fig, ax = plt.subplots(figsize=FIG_SIZE)

for color, (lo, hi) in zip(SLICE_COLORS, MJJ_SLICES):
    tag = slice_title(lo, hi)
    for s in samples_slices:
        pvm = load_pvm(str(s['file']), lo, hi)
        h, herr = norm_hist(pvm)
        ax.hist(bin_centers, bins=bins, weights=h, histtype='step',
                linewidth=2.5, color=color, linestyle=s['ls'],
                label=f"{s['label']}, {tag}")
        ax.errorbar(bin_centers, h, yerr=herr, fmt='none',
                    color=color, capsize=3, capthick=1, elinewidth=1)

apply_style(ax,
            xlabel=r'$|\vec{t}\,|$',
            ylabel='Normalized density',
            title='',
            xlim=(0, 0.02),
            log_y=True,
            legend_loc='upper right')
ax.legend(loc='upper right', frameon=False, fontsize=11, ncol=1)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(FIGS / f'{SCRIPT}_slices_mjj.{ext}', bbox_inches='tight')
plt.close(fig)
print(f'  saved figs/{SCRIPT}_slices_mjj.pdf + .png')

print(f'\nAll figures saved to {FIGS.resolve()}')
