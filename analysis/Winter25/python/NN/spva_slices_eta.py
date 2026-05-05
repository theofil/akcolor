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

MJJ_CUT = 400

ETA_SLICES = [
    (0.0, 1.5),
    (1.5, 2.5),
    (2.5, None),   # up to acceptance edge ~3.0
]

SLICE_COLORS = ['black', 'steelblue', 'crimson']

samples = [
    {'label': 'QCD h', 'file': DATA_DIR / 'QCDHtoInv.root', 'ls': '-'},
    {'label': 'VBF h', 'file': DATA_DIR / 'VBFHtoInv.root', 'ls': '--'},
]

bins        = np.linspace(0, np.pi, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width   = bins[1] - bins[0]


def load(path, eta_lo, eta_hi):
    with uproot.open(str(path)) as f:
        mjj     = f['events']['mjj'].array(library='ak')
        jetSPVA = f['events']['jetSPVA'].array(library='ak')
        jetEta  = f['events']['jetEta'].array(library='ak')

    def pad_col(arr, i):
        return ak.to_numpy(ak.fill_none(ak.pad_none(arr, i + 1, axis=1)[:, i], np.nan))

    spva0  = pad_col(jetSPVA, 0)
    eta0   = np.abs(pad_col(jetEta, 0))
    mjj_np = ak.to_numpy(mjj)
    mask = (mjj_np > MJJ_CUT) & (eta0 >= eta_lo) & ~np.isnan(spva0) & ~np.isnan(eta0)
    if eta_hi is not None:
        mask &= eta0 < eta_hi
    return np.abs(spva0[mask])


def norm_hist(values):
    counts, _ = np.histogram(values, bins=bins)
    norm = counts.sum() * bin_width
    if norm > 0:
        return counts / norm, np.sqrt(counts) / norm
    return np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)


def slice_label(lo, hi):
    if hi is None:
        return fr'$|\eta(j_1)| > {lo}$'
    return fr'${lo} \leq |\eta(j_1)| < {hi}$'


fig, ax = plt.subplots(figsize=FIG_SIZE)

for color, (lo, hi) in zip(SLICE_COLORS, ETA_SLICES):
    tag = slice_label(lo, hi)
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
            xlim=(0, np.pi),
            ylim=(0.22, 0.45),
            legend_loc='upper right')

ax.legend(loc='upper right', frameon=False, fontsize=11, ncol=1)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(FIGS / f'{SCRIPT}_eta.{ext}', bbox_inches='tight')
plt.close(fig)
print(f'Saved {SCRIPT}_eta.pdf/png to {FIGS}')
