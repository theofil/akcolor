import os
import sys
import uproot
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

FIGS = Path('figs')
FIGS.mkdir(exist_ok=True)

samples = [
    {'name': 'QCDHtoInv', 'label': 'QCD Hjj', 'file': '../QCDHtoInv.root'},
    {'name': 'QCDZtoInv', 'label': 'QCD Zjj', 'file': '../QCDZtoInv.root'},
    {'name': 'QCDWtoInv', 'label': 'QCD Wjj', 'file': '../QCDWtoInv.root'},
    {'name': 'VBFHtoInv', 'label': 'VBF Hqq', 'file': '../VBFHtoInv.root'},
    {'name': 'VBFZtoInv', 'label': 'VBF Zqq', 'file': '../VBFZtoInv.root'},
    {'name': 'VBFWtoInv', 'label': 'VBF Wqq', 'file': '../VBFWtoInv.root'},
]

branches = ['mjj', 'kWeight', 'ptjj', 'dYjj', 'dPhijj']

data = []
for sample in samples:
    with uproot.open(sample['file']) as f:
        arrays = f['events'].arrays(branches, library='np')
    entry = {'label': sample['label']}
    entry.update(arrays)
    data.append(entry)
    print(f"Loaded {len(arrays['mjj'])} events for {sample['label']}")

print('All samples loaded.\n')

# ── style ──────────────────────────────────────────────────────────────────────
sample_style = {
    'QCD Hjj': {'color': 'black',   'linestyle': '-'},
    'QCD Zjj': {'color': '#1f77b4', 'linestyle': '-'},
    'QCD Wjj': {'color': '#d62728', 'linestyle': '-'},
    'VBF Hqq': {'color': 'black',   'linestyle': ':'},
    'VBF Wqq': {'color': '#d62728', 'linestyle': ':'},
    'VBF Zqq': {'color': '#1f77b4', 'linestyle': ':'},
}

plt.rcParams.update({
    'font.size': 20, 'axes.labelsize': 22,
    'xtick.labelsize': 20, 'ytick.labelsize': 20, 'legend.fontsize': 18,
    'axes.linewidth': 1.2,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True, 'legend.frameon': False,
})


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'  saved figs/{name}.pdf + .png')


def norm_hist(values, weights, bins):
    hist, _ = np.histogram(values, bins=bins, weights=weights)
    hist_err = np.sqrt(np.histogram(values, bins=bins, weights=weights**2)[0])
    total = hist.sum()
    if total > 0:
        hist /= total
        hist_err /= total
    return hist, hist_err


def plot_shapes(var, bins, xlabel, legend_loc, legend_anchor, name,
                mjj_cut=None, log_scale=True, transform=None):
    bin_centers = (bins[:-1] + bins[1:]) / 2
    fig = plt.figure(figsize=(6.944, 6.625), dpi=72)
    ax = fig.add_axes([0.12, 0.08, 0.83, 0.87])
    all_nonzero = []

    for s in data:
        style = sample_style[s['label']]
        mask = s['mjj'] > mjj_cut if mjj_cut is not None else np.ones(len(s['mjj']), dtype=bool)
        vals = s[var][mask] if transform is None else transform(s[var][mask])
        hist, hist_err = norm_hist(vals, s['kWeight'][mask], bins)
        all_nonzero.extend(hist[hist > 0])
        step = np.concatenate([hist, [hist[-1]]])
        ax.step(bins, step, where='post', linewidth=3.5,
                color=style['color'], linestyle=style['linestyle'], label=s['label'])
        ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                    color=style['color'], capsize=3, capthick=1, elinewidth=1.5)

    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel('Normalized shape', fontsize=24)
    ax.set_xlim(bins[0], bins[-1])
    if log_scale:
        ax.set_yscale('log')
        if all_nonzero:
            arr = np.array(all_nonzero)
            ax.set_ylim(max(arr.min() / 2, 1e-8), arr.max() * 1.1)
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.legend(loc=legend_loc, bbox_to_anchor=legend_anchor,
              fontsize=18, handlelength=2, labelspacing=0.2)
    save(fig, name)


# ── plot 1: mjj 0–400, absolute (log) ─────────────────────────────────────────
print('Plotting mjj 0-400...')
bin_edges = np.linspace(0, 400, 41)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
fig = plt.figure(figsize=(6.944, 6.625), dpi=72)
ax = fig.add_axes([0.12, 0.08, 0.83, 0.87])
for s in data:
    style = sample_style[s['label']]
    hist, _ = np.histogram(s['mjj'], bins=bin_edges, weights=s['kWeight'])
    hist_err = np.sqrt(np.histogram(s['mjj'], bins=bin_edges, weights=s['kWeight']**2)[0])
    step = np.concatenate([hist, [hist[-1]]])
    ax.step(bin_edges, step, where='post', linewidth=3.5,
            color=style['color'], linestyle=style['linestyle'], label=s['label'])
    ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                color=style['color'], capsize=3, capthick=1, elinewidth=1.5)
ax.set_xlabel('m(jj) [GeV]', fontsize=24)
ax.set_ylabel('Events / bin', fontsize=24)
ax.set_xlim(0, 400)
ax.set_xticks(np.arange(0, 401, 100))
ax.set_yscale('log')
ax.set_ylim(1e-3, 1e2)
ax.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2])
ax.grid(True, linestyle=':', alpha=0.4)
ax.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98),
          fontsize=18, handlelength=2, labelspacing=0.2)
save(fig, 'mjj_0_400')

print('\nSummary Statistics:')
for s in data:
    print(f"  {s['label']}: {s['kWeight'].sum():.2f} events ({len(s['mjj'])} entries)")

# ── plot 2: mjj 200–1000 normalised ───────────────────────────────────────────
print('\nPlotting mjj 200-1000 normalised...')
plot_shapes('mjj', np.linspace(200, 1000, 81), 'm(jj) [GeV]',
            'upper right', (0.98, 0.98), 'mjj_200_1000_norm', mjj_cut=200)

# ── plot 3: ptjj ──────────────────────────────────────────────────────────────
print('Plotting ptjj...')
plot_shapes('ptjj', np.linspace(0, 1000, 81), 'p$_T$(jj) [GeV]',
            'upper right', (0.98, 0.98), 'ptjj_norm', mjj_cut=200)

# ── plot 4: dYjj ──────────────────────────────────────────────────────────────
print('Plotting dYjj...')
plot_shapes('dYjj', np.linspace(0, 5, 11), '$|\\Delta y(jj)|$',
            'lower right', (0.98, 0.02), 'dYjj_norm', mjj_cut=200, transform=np.abs)

# ── plot 5: dPhijj ────────────────────────────────────────────────────────────
print('Plotting dPhijj...')
plot_shapes('dPhijj', np.linspace(0, np.pi, 11), '$|\\Delta\\phi(jj)|$ [rad]',
            'upper left', (0.02, 0.98), 'dPhijj_norm', mjj_cut=200, transform=np.abs)

viewer = Path.home() / 'qplots' / 'viewer.py'
print(f'\nAll figures saved to {FIGS.resolve()}')
print(f'Launching viewer...')
os.execv(sys.executable, [sys.executable, str(viewer), str(FIGS.resolve())])
