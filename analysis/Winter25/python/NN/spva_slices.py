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

MJJ_SLICES = [
    (400,  550),
    (550,  700),
    (700,  900),
    (900,  None),
]

SLICE_COLORS = ['black', 'steelblue', 'darkgreen', 'crimson']

samples = [
    {'label': 'QCD h', 'file': DATA_DIR / 'QCDHtoInv.root', 'ls': '-'},
    {'label': 'VBF h', 'file': DATA_DIR / 'VBFHtoInv.root', 'ls': '--'},
]

bins        = np.linspace(0, np.pi, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width   = bins[1] - bins[0]


def load(path, mjj_min, mjj_max):
    with uproot.open(str(path)) as f:
        mjj     = f['events']['mjj'].array(library='ak')
        jetSPVA = f['events']['jetSPVA'].array(library='ak')
    mask = mjj > mjj_min
    if mjj_max is not None:
        mask = mask & (mjj < mjj_max)
    spva0 = ak.to_numpy(ak.fill_none(jetSPVA[mask, 0], np.nan))
    return np.abs(spva0[~np.isnan(spva0)])


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


fig, ax = plt.subplots(figsize=FIG_SIZE)

for color, (lo, hi) in zip(SLICE_COLORS, MJJ_SLICES):
    tag = slice_title(lo, hi)
    for s in samples:
        spva = load(s['file'], lo, hi)
        h, herr = norm_hist(spva)
        ax.hist(bin_centers, bins=bins, weights=h, histtype='step',
                linewidth=2.5, color=color, linestyle=s['ls'],
                label=f"{s['label']}, {tag}")
        ax.errorbar(bin_centers, h, yerr=herr, fmt='none',
                    color=color, capsize=3, capthick=1, elinewidth=1)

apply_style(ax,
            xlabel=r'$|\theta_s|$ [rad]',
            ylabel='Normalized density [1/rad]',
            title='',
            xlim=(0, np.pi), ylim=(0.25, 0.47),
            legend_loc='upper right')

ax.legend(loc='upper right', frameon=False, fontsize=11, ncol=1)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(FIGS / f'{SCRIPT}_mjj.{ext}', bbox_inches='tight')
plt.close(fig)
print(f'Saved {SCRIPT}_mjj.pdf/png to {FIGS}')
