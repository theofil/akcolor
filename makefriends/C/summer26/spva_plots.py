#!/usr/bin/env python3
"""Herwig-only comparison: default color-reconnection probability (0.95, friends/summer26)
vs the precoprob variant (0.25, friends/summer26/precoprob). No Pythia involved -- CRP is
a Herwig-specific parton-shower tunable, so there is nothing to compare it against."""
import pathlib
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import uproot

sys.path.insert(0, '/afs/cern.ch/user/t/theofil/work/akcolor/analysis/Winter25/python/NN')
from style import apply_style, FIG_SIZE

def separation_power(p, q, bin_width):
    """<S^2> = 0.5 * integral[(p-q)^2/(p+q)] dx for two normalized PDFs."""
    denom = p + q
    mask  = denom > 0
    return 0.5 * float(np.sum(np.where(mask, (p - q)**2 / np.where(mask, denom, 1.0), 0.0)) * bin_width)

def kl_divergence(p, q):
    """D_KL(p||q) for normalized probability mass arrays (sum=1)."""
    mask = (p > 0) & (q > 0)
    safe_ratio = np.where(mask, p / np.where(mask, q, 1.0), 1.0)
    return float(np.sum(np.where(mask, p * np.log(safe_ratio), 0.0)))

def counts_from_tree(fpath, branch, bins):
    with uproot.open(str(fpath)) as rf:
        if 'events' not in rf:
            return None
        vals = rf['events'][branch].array(library='np').astype(float)
    counts, _ = np.histogram(vals, bins=bins)
    return counts.astype(float)

DATA_DIR_DEFAULT   = pathlib.Path(__file__).parent.parent / 'friends' / 'summer26'
DATA_DIR_PRECOPROB = DATA_DIR_DEFAULT / 'precoprob'
FIGS               = pathlib.Path(__file__).parent.parent / 'figs' / 'summer26' / 'precoprob'
FIGS.mkdir(parents=True, exist_ok=True)

# ── Sample catalog: Herwig only, one VBF + one QCD sample per boson channel ──────
SAMPLES_BY_PROC = {
    'H': [
        {'name': 'QCDHjj_herwig', 'label': 'QCD H  Herwig', 'group': 'QCD'},
        {'name': 'VBFH_herwig',   'label': 'VBF H  Herwig', 'group': 'VBF'},
    ],
    'W': [
        {'name': 'QCDWjj_herwig', 'label': 'QCD W  Herwig', 'group': 'QCD'},
        {'name': 'VBFW_herwig',   'label': 'VBF W  Herwig', 'group': 'VBF'},
    ],
    'Z': [
        {'name': 'QCDZjj_herwig', 'label': 'QCD Z  Herwig', 'group': 'QCD'},
        {'name': 'VBFZ_herwig',   'label': 'VBF Z  Herwig', 'group': 'VBF'},
    ],
}
# color by group (VBF/QCD), linestyle by CR-probability variant
group_colors  = {'VBF': 'steelblue', 'QCD': 'tomato'}
variant_styles = {'default': '-', 'precoprob': '--'}
variant_labels = {'default': '(default)', 'precoprob': '(modified)'}
DATA_DIRS = {'default': DATA_DIR_DEFAULT, 'precoprob': DATA_DIR_PRECOPROB}


def _density_from_hist(fpath, hist_name, bin_width):
    if not fpath.exists():
        print(f'  Skipping {fpath} (not found)')
        return None, None
    with uproot.open(str(fpath)) as rf:
        if hist_name not in rf:
            print(f'  {hist_name} not found in {fpath}')
            return None, None
        h = rf[hist_name]
        counts, variances = h.values(), h.variances()
    norm = counts.sum() * bin_width
    if norm <= 0:
        return None, None
    return counts / norm, np.sqrt(variances) / norm


def density_comparison(hist_name, bins, xlabel, ylabel, xlim, ylim, out_stub, log_y=False,
                        legend_loc='upper right'):
    """One panel per boson channel: VBF+QCD overlaid, default (solid) vs precoprob (dashed)."""
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_width   = bins[1] - bins[0]
    for proc, samples in SAMPLES_BY_PROC.items():
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        for s in samples:
            for variant, data_dir in DATA_DIRS.items():
                fpath = data_dir / f'{s["name"]}.friend.root'
                hist, hist_err = _density_from_hist(fpath, hist_name, bin_width)
                if hist is None:
                    continue
                ax.hist(bin_centers, bins=bins, weights=hist, histtype='step',
                        linewidth=3, color=group_colors[s['group']], linestyle=variant_styles[variant],
                        label=f'{s["label"]}  {variant_labels[variant]}')
                ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                            color=group_colors[s['group']], capsize=3, capthick=1, elinewidth=1)
        apply_style(ax, xlabel=xlabel, ylabel=ylabel, title='',
                    xlim=xlim, ylim=ylim, log_y=log_y, legend_loc=legend_loc)
        plt.tight_layout()
        outpath = str(FIGS / f'{out_stub}_{proc}.pdf')
        fig.savefig(outpath, bbox_inches='tight')
        print('Saved', outpath)
        plt.close(fig)


def ratio_comparison(hist_name, bins, xlabel, out_stub, ylim=(0.4, 1.6)):
    """precoprob/default ratio, one curve per group (VBF, QCD), per boson channel."""
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_width   = bins[1] - bins[0]
    group_colors_ratio = {'VBF': 'steelblue', 'QCD': 'tomato'}
    for proc, samples in SAMPLES_BY_PROC.items():
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        for s in samples:
            default_hist, default_err = _density_from_hist(
                DATA_DIR_DEFAULT / f'{s["name"]}.friend.root', hist_name, bin_width)
            preco_hist, preco_err = _density_from_hist(
                DATA_DIR_PRECOPROB / f'{s["name"]}.friend.root', hist_name, bin_width)
            if default_hist is None or preco_hist is None:
                continue
            mask = default_hist > 0
            ratio = np.where(mask, preco_hist / np.where(mask, default_hist, 1), np.nan)
            ratio_err = np.where(mask,
                                 ratio * np.sqrt(
                                     (np.where(mask, preco_err / np.where(preco_hist > 0, preco_hist, 1), 0))**2 +
                                     (np.where(mask, default_err / np.where(default_hist > 0, default_hist, 1), 0))**2),
                                 np.nan)
            rms  = np.sqrt(np.nanmean(ratio**2))
            sep2 = separation_power(preco_hist, default_hist, bin_width)
            kld  = kl_divergence(preco_hist * bin_width, default_hist * bin_width)
            ax.hist(bin_centers, bins=bins, weights=np.where(np.isnan(ratio), 0, ratio),
                    histtype='step', linewidth=3, color=group_colors_ratio[s['group']],
                    label=f'{s["label"]}  RMS={rms:.3f}  $\\langle S^2\\rangle$={sep2:.3f}  KL={kld:.3f}')
            ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt='none',
                        color=group_colors_ratio[s['group']], capsize=3, capthick=1, elinewidth=1)
        ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')
        apply_style(ax, xlabel=xlabel, ylabel='precoprob / default', title='',
                    xlim=(bins[0], bins[-1]), ylim=ylim, legend_loc='upper right')
        plt.tight_layout()
        outpath = str(FIGS / f'{out_stub}_{proc}.pdf')
        fig.savefig(outpath, bbox_inches='tight')
        print('Saved', outpath)
        plt.close(fig)


# ── Leading-jet |SPVA| (pull-vector angle) comparison ────────────────────────────
spva_bins = np.linspace(0, np.pi, 21)
density_comparison('hLeadJetSPVA', spva_bins,
                    xlabel=r'$|\theta_s|$ [rad]', ylabel='Normalized density [1/rad]',
                    xlim=(0, np.pi), ylim=(0.25, 0.47), out_stub='jetSPVA')
# jetSPVA ratio intentionally omitted (per instruction)

# ── Leading-jet PVM |t⃗| (pull-vector magnitude) comparison + ratio ──────────────
spvm_bins = np.linspace(0, 0.06, 21)
density_comparison('hLeadJetSPVM', spvm_bins,
                    xlabel=r'$|\vec{t}\,|$', ylabel='Normalized density',
                    xlim=(0, 0.06), ylim=None, out_stub='jetSPVM', log_y=True)
ratio_comparison('hLeadJetSPVM', spvm_bins, xlabel=r'$|\vec{t}\,|$', out_stub='jetSPVM_ratio')

# ── Leading-jet NC (number of constituents) comparison + ratio ──────────────────
nc_bins = np.linspace(0, 100, 26)


def _nc_density(fpath, bin_width):
    counts = counts_from_tree(fpath, 'jetNC0', nc_bins)
    if counts is None:
        return None, None
    norm = counts.sum() * bin_width
    if norm <= 0:
        return None, None
    return counts / norm, np.sqrt(counts) / norm


nc_bin_centers = (nc_bins[:-1] + nc_bins[1:]) / 2
nc_bin_width   = nc_bins[1] - nc_bins[0]

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        for variant, data_dir in DATA_DIRS.items():
            fpath = data_dir / f'{s["name"]}.friend.root'
            if not fpath.exists():
                print(f'  Skipping {fpath} (not found)')
                continue
            hist, hist_err = _nc_density(fpath, nc_bin_width)
            if hist is None:
                print(f'  jetNC0 not found/empty in {fpath}')
                continue
            ax.hist(nc_bin_centers, bins=nc_bins, weights=hist, histtype='step',
                    linewidth=3, color=group_colors[s['group']], linestyle=variant_styles[variant],
                    label=f'{s["label"]}  {variant_labels[variant]}')
            ax.errorbar(nc_bin_centers, hist, yerr=hist_err, fmt='none',
                        color=group_colors[s['group']], capsize=3, capthick=1, elinewidth=1)
    apply_style(ax, xlabel=r'$N_\mathrm{constituents}$', ylabel='Normalized density', title='',
                xlim=(0, 100), legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetNC_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

group_colors_ratio = {'VBF': 'steelblue', 'QCD': 'tomato'}
for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        default_counts = counts_from_tree(DATA_DIR_DEFAULT / f'{s["name"]}.friend.root', 'jetNC0', nc_bins)
        preco_counts   = counts_from_tree(DATA_DIR_PRECOPROB / f'{s["name"]}.friend.root', 'jetNC0', nc_bins)
        if default_counts is None or preco_counts is None:
            print(f'  jetNC0 missing for {s["name"]}, skipping ratio')
            continue
        default_norm = default_counts.sum() * nc_bin_width
        preco_norm   = preco_counts.sum() * nc_bin_width
        if default_norm <= 0 or preco_norm <= 0:
            continue
        default_hist = default_counts / default_norm
        preco_hist   = preco_counts / preco_norm
        default_err  = np.sqrt(default_counts) / default_norm
        preco_err    = np.sqrt(preco_counts) / preco_norm
        mask = default_hist > 0
        ratio = np.where(mask, preco_hist / np.where(mask, default_hist, 1), np.nan)
        ratio_err = np.where(mask,
                             ratio * np.sqrt(
                                 (np.where(mask, preco_err / np.where(preco_hist > 0, preco_hist, 1), 0))**2 +
                                 (np.where(mask, default_err / np.where(default_hist > 0, default_hist, 1), 0))**2),
                             np.nan)
        rms  = np.sqrt(np.nanmean(ratio**2))
        sep2 = separation_power(preco_hist, default_hist, nc_bin_width)
        kld  = kl_divergence(preco_hist * nc_bin_width, default_hist * nc_bin_width)
        ax.hist(nc_bin_centers, bins=nc_bins, weights=np.where(np.isnan(ratio), 0, ratio),
                histtype='step', linewidth=3, color=group_colors_ratio[s['group']],
                label=f'{s["label"]}  RMS={rms:.3f}  $\\langle S^2\\rangle$={sep2:.3f}  KL={kld:.3f}')
        ax.errorbar(nc_bin_centers, ratio, yerr=ratio_err, fmt='none',
                    color=group_colors_ratio[s['group']], capsize=3, capthick=1, elinewidth=1)
    ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')
    apply_style(ax, xlabel=r'$N_\mathrm{constituents}$', ylabel='precoprob / default', title='',
                xlim=(0, 100), ylim=(0.4, 1.6), legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetNC_ratio_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── index.php (same lightweight listing convention as figs/summer26/index.php) ──
index_path = FIGS / 'index.php'
php = '''<?php
$pdfs = glob("*.pdf");
sort($pdfs);
echo "<html><body>";
echo "<h2>akcolor precoprob (CRP 0.95 vs 0.25, Herwig only) plots</h2><ul>";
foreach ($pdfs as $f) {
    echo "<li><a href=\\"$f\\">$f</a></li>";
}
echo "</ul></body></html>";
?>
'''
index_path.write_text(php)
print('Updated', index_path)
