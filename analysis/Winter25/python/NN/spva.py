import os
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
DATA_DIR = NN_DIR.parents[1]   # Winter25/

FIGS.mkdir(exist_ok=True)

samples = [
    {'name': 'QCDHtoInv', 'label': 'QCD h', 'file': str(DATA_DIR / 'QCDHtoInv.root'), 'group': 'QCD'},
    {'name': 'QCDZtoInv', 'label': 'QCD Z', 'file': str(DATA_DIR / 'QCDZtoInv.root'), 'group': 'QCD'},
    {'name': 'QCDWtoInv', 'label': 'QCD W', 'file': str(DATA_DIR / 'QCDWtoInv.root'), 'group': 'QCD'},
    {'name': 'VBFHtoInv', 'label': 'VBF h', 'file': str(DATA_DIR / 'VBFHtoInv.root'), 'group': 'VBF'},
    {'name': 'VBFZtoInv', 'label': 'VBF Z', 'file': str(DATA_DIR / 'VBFZtoInv.root'), 'group': 'VBF'},
    {'name': 'VBFWtoInv', 'label': 'VBF W', 'file': str(DATA_DIR / 'VBFWtoInv.root'), 'group': 'VBF'},
]

mjj_cut   = 400
bins      = np.linspace(0, np.pi, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width = bins[1] - bins[0]

colors     = {'H': 'black', 'Z': '#1f77b4', 'W': '#d62728'}
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


# ── plot 1: per-sample leading-jet |θ_s| for mjj > mjj_cut ───────────────────
print(f'Plotting per-sample SPVA (mjj > {mjj_cut} GeV)...')
fig, ax = plt.subplots(figsize=FIG_SIZE)

for s in samples:
    label = s['label']
    key   = 'H' if 'H' in label else 'Z' if 'Z' in label else 'W'
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
            xlim=(0, np.pi), ylim=(0.25, 0.45),
            legend_loc='upper right')
plt.tight_layout()
save(fig, f'leading_jet_mjj_gt_{mjj_cut}')

# ── plot 2: QCD vs VBF merged, two mjj windows ────────────────────────────────
mjj_low_min, mjj_low_max = 400, 600
mjj_high_min             = 800

mjj_selections = [
    {'label': f'{mjj_low_min}<mjj<{mjj_low_max}', 'min': mjj_low_min, 'max': mjj_low_max,
     'color': 'black', 'alpha': 0.6},
    {'label': f'mjj>{mjj_high_min}',               'min': mjj_high_min, 'max': None,
     'color': 'red',   'alpha': 0.6},
]

print('Plotting merged QCD vs VBF for two mjj windows...')
fig, ax = plt.subplots(figsize=FIG_SIZE)

for sel in mjj_selections:
    for group, linestyle in [('QCD', '-'), ('VBF', ':')]:
        hist_total = np.zeros(len(bins) - 1)
        n_entries  = 0
        for s in samples:
            if s['group'] != group:
                continue
            spva = load_spva(s['file'], sel['min'], sel['max'])
            counts, _ = np.histogram(spva, bins=bins)
            hist_total += counts
            n_entries  += len(spva)
        norm = hist_total.sum() * bin_width
        if norm > 0:
            hist     = hist_total / norm
            hist_err = np.sqrt(hist_total) / norm
        else:
            hist     = np.zeros_like(hist_total)
            hist_err = np.zeros_like(hist_total)
        label = f'{group} ({sel["label"]})'
        ax.hist(bin_centers, bins=bins, weights=hist, histtype='step', linewidth=3,
                color=sel['color'], linestyle=linestyle, alpha=sel['alpha'], label=label)
        ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                    color=sel['color'], alpha=sel['alpha'], capsize=3, capthick=1, elinewidth=1)
        print(f'  {label}: {n_entries} leading-jet entries')

apply_style(ax,
            xlabel=r'$|\theta_s|$ [rad]',
            ylabel='Normalized density [1/rad]',
            title='',
            xlim=(0, np.pi), ylim=(0.25, 0.45),
            legend_loc='upper right')
plt.tight_layout()
save(fig, 'QCD_vs_VBF_merged')

viewer = Path.home() / 'qplots' / 'viewer.py'
print(f'\nAll figures saved to {FIGS.resolve()}')
print('Launching viewer...')
os.execv(sys.executable, [sys.executable, str(viewer), str(FIGS.resolve())])
