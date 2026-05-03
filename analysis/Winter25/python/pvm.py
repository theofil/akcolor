import os
import sys
import uproot
import awkward as ak
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

FIGS = Path('figs')
FIGS.mkdir(exist_ok=True)

samples = [
    {'name': 'QCDHtoInv', 'label': 'QCD Hjj', 'file': '../QCDHtoInv.root', 'group': 'QCD'},
    {'name': 'QCDZtoInv', 'label': 'QCD Zjj', 'file': '../QCDZtoInv.root', 'group': 'QCD'},
    {'name': 'QCDWtoInv', 'label': 'QCD Wjj', 'file': '../QCDWtoInv.root', 'group': 'QCD'},
    {'name': 'VBFHtoInv', 'label': 'VBF Hqq', 'file': '../VBFHtoInv.root', 'group': 'VBF'},
    {'name': 'VBFZtoInv', 'label': 'VBF Zqq', 'file': '../VBFZtoInv.root', 'group': 'VBF'},
    {'name': 'VBFWtoInv', 'label': 'VBF Wqq', 'file': '../VBFWtoInv.root', 'group': 'VBF'},
]

mjj_cut = 400
bins = np.linspace(0, 0.02, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width = bins[1] - bins[0]

colors = {'H': 'black', 'Z': '#1f77b4', 'W': '#d62728'}
line_styles = {'QCD': '-', 'VBF': ':'}


def load_pvm(file_path, mjj_min, mjj_max=None):
    with uproot.open(file_path) as f:
        tree = f['events']
        mjj   = tree['mjj'].array(library='ak')
        jetPVM = tree['jetPVM'].array(library='ak')
    mask = mjj > mjj_min
    if mjj_max is not None:
        mask = mask & (mjj < mjj_max)
    pvm0 = jetPVM[mask, 0]
    pvm0 = ak.to_numpy(ak.fill_none(pvm0, np.nan))
    return pvm0[~np.isnan(pvm0)]


def norm_hist(values, bins):
    counts, _ = np.histogram(values, bins=bins)
    norm = counts.sum() * bin_width
    if norm > 0:
        return counts / norm, np.sqrt(counts) / norm
    return np.zeros_like(counts, dtype=float), np.zeros_like(counts, dtype=float)


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved figs/{name}.pdf + .png')


# ── plot 1: per-sample leading-jet PVM for mjj > mjj_cut ─────────────────────
print(f'Plotting per-sample PVM (mjj > {mjj_cut} GeV)...')
fig, ax = plt.subplots(figsize=(6, 6))

for s in samples:
    label = s['label']
    key = 'H' if 'H' in label else 'Z' if 'Z' in label else 'W'
    color = colors[key]
    linestyle = line_styles[s['group']]

    pvm = load_pvm(s['file'], mjj_cut)
    hist, hist_err = norm_hist(pvm, bins)

    ax.hist(bin_centers, bins=bins, weights=hist, histtype='step',
            linewidth=3, color=color, linestyle=linestyle, label=label)
    ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                color=color, capsize=3, capthick=1, elinewidth=1)

ax.set_xlabel(r'PVM$[0]$', fontsize=16)
ax.set_ylabel('Normalized density', fontsize=16)
ax.set_title(f'Leading jet PVM for $m_{{jj}} > {mjj_cut}$ GeV', fontsize=18)
ax.set_xlim(bins[0], bins[-1])
ax.set_yscale('log')
ax.grid(True, linestyle=':', alpha=0.4)
ax.legend(loc='upper right', frameon=False, fontsize=12)
plt.tight_layout()
save(fig, f'PVM_leading_jet_mjj_gt_{mjj_cut}')


viewer = Path.home() / 'qplots' / 'viewer.py'
print(f'\nAll figures saved to {FIGS.resolve()}')
print('Launching viewer...')
os.execv(sys.executable, [sys.executable, str(viewer), str(FIGS.resolve())])
