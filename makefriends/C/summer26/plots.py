#!/usr/bin/env python3
import pathlib
import subprocess
import sys
import ROOT

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import uproot
import h5py

sys.path.insert(0, '/afs/cern.ch/user/t/theofil/work/akcolor/analysis/Winter25/python/NN')
from style import apply_style, FIG_SIZE, LABEL_SIZE, TITLE_SIZE, TICK_SIZE

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

def roc_from_histograms(sig_counts, bkg_counts):
    """ROC curve from histogram counts; returns (fpr, tpr, auc) scanning all thresholds."""
    sig_total = sig_counts.sum()
    bkg_total = bkg_counts.sum()
    if sig_total <= 0 or bkg_total <= 0:
        return None, None, None
    sig_cdf = np.concatenate([np.cumsum(sig_counts[::-1])[::-1], [0]]) / sig_total
    bkg_cdf = np.concatenate([np.cumsum(bkg_counts[::-1])[::-1], [0]]) / bkg_total
    tpr = sig_cdf[::-1]
    fpr = bkg_cdf[::-1]
    auc = float(np.trapz(tpr, fpr))
    if auc < 0.5:
        tpr = 1 - tpr[::-1]
        fpr = 1 - fpr[::-1]
        auc = 1 - auc
    return fpr, tpr, auc

def counts_from_tree(fpath, branch, bins, transform=None):
    """Histogram a scalar jet branch from the events tree."""
    with uproot.open(str(fpath)) as rf:
        if 'events' not in rf:
            return None
        vals = rf['events'][branch].array(library='np').astype(float)
    if transform is not None:
        vals = transform(vals)
    counts, _ = np.histogram(vals, bins=bins)
    return counts.astype(float)

TDR_STYLE = pathlib.Path.home() / 'qplot' / 'setTDRStyle.C'
DATA_DIR   = pathlib.Path(__file__).parent.parent / 'friends' / 'summer26'
FIGS       = pathlib.Path(__file__).parent.parent / 'figs' / 'summer26'

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.LoadMacro(str(TDR_STYLE))
ROOT.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

FIGS.mkdir(parents=True, exist_ok=True)

SKIP_HISTS   = {'hLeadJetSPVA', 'hLeadJetSPVM', 'hJetPt1', 'hNJets', 'jetShapeRaw'}
HIST_CLASSES = {'TH1F', 'TH2F', 'TH2Poly'}

for root_file in sorted(DATA_DIR.glob('*friend.root')):
    stem = root_file.stem.replace('.friend', '')
    print('Processing', root_file)

    f = ROOT.TFile(str(root_file))
    if not f or f.IsZombie():
        print('  Could not open', root_file)
        continue

    for key in f.GetListOfKeys():
        classname = key.GetClassName()
        if classname not in HIST_CLASSES:
            continue
        obj = key.ReadObj()
        hname = obj.GetName()
        if hname in SKIP_HISTS or hname.startswith('jetSPVA_'):
            continue
        obj.SetTitle('')
        is2d = classname in ('TH2F', 'TH2Poly')
        c = ROOT.TCanvas(f'c_{stem}_{hname}', '', 600 if is2d else 800, 600)
        if is2d:
            c.SetRightMargin(0.15)
            obj.Draw('COLZ')
            if classname == 'TH2Poly':
                obj.GetXaxis().SetTitle('#Delta y_{i}')
                obj.GetYaxis().SetTitle('#Delta#varphi_{i}')
                obj.GetXaxis().SetNdivisions(505)
                obj.GetYaxis().SetNdivisions(505)
        else:
            obj.Draw('HIST')
        outpath = str(FIGS / f'{hname}_{stem}.pdf')
        c.SaveAs(outpath)
        c.Close()
        print('  Saved', outpath)

    f.Close()

# ── Erase stale PDFs ──────────────────────────────────────────────────────────
for stale in FIGS.glob('*.friend.pdf'):
    stale.unlink()
    print('Removed', stale)
for pattern in ('hJetPt1_*.pdf', 'hNJets_*.pdf', 'jetShapeRaw_*.pdf'):
    for stale in FIGS.glob(pattern):
        stale.unlink()
        print('Removed', stale)
for stale_name in (
    'hLeadJetSPVA_comparison.pdf',
    'hLeadJetSPVA_comparison.png',
    'mjj.pdf', 'mjj_Z.pdf', 'dYjj.pdf', 'dPhijj.pdf', 'Ptjj.pdf',
    'jetPt.pdf', 'jetEta.pdf', 'jetNC.pdf', 'jetNC_ratio.pdf',
    'jet_roc.pdf', 'jetSPVA.pdf', 'jetSPVA_ratio.pdf',
    'jetSPVM.pdf', 'jetSPVM_ratio.pdf',
):
    p = FIGS / stale_name
    if p.exists():
        p.unlink()
        print('Removed', p)

# ── Sample catalog ─────────────────────────────────────────────────────────────
SAMPLES_BY_PROC = {
    'H': [
        {'name': 'QCDHjj_mg5_pythia', 'label': 'QCD H  Pythia', 'group': 'QCD', 'gen': 'mg5'},
        {'name': 'QCDHjj_herwig',     'label': 'QCD H  Herwig',      'group': 'QCD', 'gen': 'herwig'},
        {'name': 'VBFH_mg5_pythia',   'label': 'VBF H  Pythia',  'group': 'VBF', 'gen': 'mg5'},
        {'name': 'VBFH_herwig',       'label': 'VBF H  Herwig',      'group': 'VBF', 'gen': 'herwig'},
    ],
    'W': [
        {'name': 'QCDWjj_mg5_pythia', 'label': 'QCD W  Pythia', 'group': 'QCD', 'gen': 'mg5'},
        {'name': 'QCDWjj_herwig',     'label': 'QCD W  Herwig',      'group': 'QCD', 'gen': 'herwig'},
        {'name': 'VBFW_mg5_pythia',   'label': 'VBF W  Pythia',  'group': 'VBF', 'gen': 'mg5'},
        {'name': 'VBFW_herwig',       'label': 'VBF W  Herwig',      'group': 'VBF', 'gen': 'herwig'},
    ],
    'Z': [
        {'name': 'QCDZjj_mg5_pythia', 'label': 'QCD Z  Pythia', 'group': 'QCD', 'gen': 'mg5'},
        {'name': 'QCDZjj_herwig',     'label': 'QCD Z  Herwig',      'group': 'QCD', 'gen': 'herwig'},
        {'name': 'VBFZ_mg5_pythia',   'label': 'VBF Z  Pythia',  'group': 'VBF', 'gen': 'mg5'},
        {'name': 'VBFZ_herwig',       'label': 'VBF Z  Herwig',      'group': 'VBF', 'gen': 'herwig'},
    ],
}
# color by group (VBF/QCD), linestyle by generator
group_colors     = {'VBF': 'steelblue', 'QCD': 'tomato'}
gen_styles       = {'mg5': '--', 'herwig': '-'}
# for ratio plots: color by generator
gen_colors_ratio = {'mg5': 'black', 'herwig': 'royalblue'}
gen_labels_ratio = {'mg5': 'Pythia', 'herwig': 'Herwig'}

# ── Leading-jet |SPVA| comparison ────────────────────────────────────────────
bins        = np.linspace(0, np.pi, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width   = bins[1] - bins[0]

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'hLeadJetSPVA' not in rf:
                print(f'  hLeadJetSPVA not found in {fpath}')
                continue
            h         = rf['hLeadJetSPVA']
            counts    = h.values()
            variances = h.variances()
        norm = counts.sum() * bin_width
        if norm > 0:
            hist     = counts / norm
            hist_err = np.sqrt(variances) / norm
        else:
            hist     = np.zeros_like(counts, dtype=float)
            hist_err = np.zeros_like(counts, dtype=float)
        ax.hist(bin_centers, bins=bins, weights=hist, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
        ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                    color=group_colors[s['group']], capsize=3, capthick=1, elinewidth=1)
    apply_style(ax,
                xlabel=r'$|\theta_s|$ [rad]',
                ylabel='Normalized density [1/rad]',
                title='',
                xlim=(0, np.pi), ylim=(0.25, 0.47),
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetSPVA_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Leading-jet |SPVA| VBF/QCD ratio ─────────────────────────────────────────
for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for gen in ('mg5', 'herwig'):
        gen_samps = [s for s in samples if s['gen'] == gen]
        qcd_s = next((s for s in gen_samps if s['group'] == 'QCD'), None)
        vbf_s = next((s for s in gen_samps if s['group'] == 'VBF'), None)
        if qcd_s is None or vbf_s is None:
            continue
        qcd_counts = qcd_vars = vbf_counts = vbf_vars = None
        for s, store in [(qcd_s, 'qcd'), (vbf_s, 'vbf')]:
            fpath = DATA_DIR / f'{s["name"]}.friend.root'
            if not fpath.exists():
                print(f'  Skipping {fpath} (not found)')
                break
            with uproot.open(str(fpath)) as rf:
                if 'hLeadJetSPVA' not in rf:
                    print(f'  hLeadJetSPVA not found in {fpath}')
                    break
                h = rf['hLeadJetSPVA']
                if store == 'qcd':
                    qcd_counts, qcd_vars = h.values(), h.variances()
                else:
                    vbf_counts, vbf_vars = h.values(), h.variances()
        if qcd_counts is None or vbf_counts is None:
            continue
        qcd_norm = qcd_counts.sum() * bin_width
        vbf_norm = vbf_counts.sum() * bin_width
        if qcd_norm <= 0 or vbf_norm <= 0:
            continue
        qcd_hist = qcd_counts / qcd_norm
        vbf_hist = vbf_counts / vbf_norm
        qcd_err  = np.sqrt(qcd_vars) / qcd_norm
        vbf_err  = np.sqrt(vbf_vars) / vbf_norm
        mask = qcd_hist > 0
        ratio     = np.where(mask, vbf_hist / np.where(mask, qcd_hist, 1), np.nan)
        ratio_err = np.where(mask,
                             ratio * np.sqrt(
                                 (np.where(mask, vbf_err / np.where(vbf_hist > 0, vbf_hist, 1), 0))**2 +
                                 (np.where(mask, qcd_err / np.where(qcd_hist > 0, qcd_hist, 1), 0))**2),
                             np.nan)
        rms  = np.sqrt(np.nanmean(ratio**2))
        sep2 = separation_power(vbf_hist, qcd_hist, bin_width)
        kld  = kl_divergence(vbf_hist * bin_width, qcd_hist * bin_width)
        ax.hist(bin_centers, bins=bins, weights=np.where(np.isnan(ratio), 0, ratio),
                histtype='step', linewidth=3,
                color=gen_colors_ratio[gen], linestyle=gen_styles[gen],
                label=f'{gen_labels_ratio[gen]}  RMS={rms:.3f}  $\\langle S^2\\rangle$={sep2:.3f}  KL={kld:.3f}')
        ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt='none',
                    color=gen_colors_ratio[gen], capsize=3, capthick=1, elinewidth=1)
    ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')
    apply_style(ax,
                xlabel=r'$|\theta_s|$ [rad]',
                ylabel='VBF / QCD',
                title='',
                xlim=(0, np.pi), ylim=(0.4, 1.6),
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetSPVA_ratio_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Leading-jet |SPVA| combined W/Z/H (all in one crowded diagram) ──────────
# color = boson (W/Z/H, matching index3.php palette), linestyle = VBF/QCD.
# Filenames carry no trailing _W/_Z/_H boson suffix, so index3.php's
# glob-based W/Z/H triplet grouping does not sweep these in.
boson_color_combined = {'W': '#1a6eb5', 'Z': '#2e8b57', 'H': '#b52020'}
group_ls_combined    = {'VBF': '-', 'QCD': '--'}
combined_ylim        = (0.25, 0.46)

def _read_spva_density(fpath):
    if not fpath.exists():
        print(f'  Skipping {fpath} (not found)')
        return None, None
    with uproot.open(str(fpath)) as rf:
        if 'hLeadJetSPVA' not in rf:
            print(f'  hLeadJetSPVA not found in {fpath}')
            return None, None
        h = rf['hLeadJetSPVA']
        counts, variances = h.values(), h.variances()
    norm = counts.sum() * bin_width
    if norm <= 0:
        return None, None
    return counts / norm, np.sqrt(variances) / norm

# -- Herwig only: 6 curves ----------------------------------------------------
fig, ax = plt.subplots(figsize=FIG_SIZE)
for proc, samples in SAMPLES_BY_PROC.items():
    for group in ('VBF', 'QCD'):
        s = next(s for s in samples if s['group'] == group and s['gen'] == 'herwig')
        hist, hist_err = _read_spva_density(DATA_DIR / f'{s["name"]}.friend.root')
        if hist is None:
            continue
        ax.hist(bin_centers, bins=bins, weights=hist, histtype='step',
                linewidth=3, color=boson_color_combined[proc], linestyle=group_ls_combined[group],
                label=f'{proc} {group} Herwig')
        ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                    color=boson_color_combined[proc], capsize=3, capthick=1, elinewidth=1)
apply_style(ax,
            xlabel=r'$|\theta_s|$ [rad]',
            ylabel='Normalized density [1/rad]',
            title='',
            xlim=(0, np.pi), ylim=combined_ylim,
            legend_loc='upper right')
plt.tight_layout()
outpath = str(FIGS / 'jetSPVA_combined_herwig.pdf')
fig.savefig(outpath, bbox_inches='tight')
print('Saved', outpath)
plt.close(fig)

# -- Both generators: 12 curves, generator encoded via alpha/linewidth -------
gen_alpha_combined = {'mg5': 1.0, 'herwig': 0.5}
gen_lw_combined     = {'mg5': 3.0, 'herwig': 1.8}
gen_label_combined  = {'mg5': 'Pythia', 'herwig': 'Herwig'}

fig, ax = plt.subplots(figsize=FIG_SIZE)
for proc, samples in SAMPLES_BY_PROC.items():
    for group in ('VBF', 'QCD'):
        for gen in ('mg5', 'herwig'):
            s = next(s for s in samples if s['group'] == group and s['gen'] == gen)
            hist, hist_err = _read_spva_density(DATA_DIR / f'{s["name"]}.friend.root')
            if hist is None:
                continue
            ax.hist(bin_centers, bins=bins, weights=hist, histtype='step',
                    linewidth=gen_lw_combined[gen], alpha=gen_alpha_combined[gen],
                    color=boson_color_combined[proc], linestyle=group_ls_combined[group],
                    label=f'{proc} {group} {gen_label_combined[gen]}')
apply_style(ax,
            xlabel=r'$|\theta_s|$ [rad]',
            ylabel='Normalized density [1/rad]',
            title='',
            xlim=(0, np.pi), ylim=combined_ylim,
            legend_loc='upper right')
ax.legend(loc='upper right', frameon=False, fontsize=7, ncol=3)
plt.tight_layout()
outpath = str(FIGS / 'jetSPVA_combined_all.pdf')
fig.savefig(outpath, bbox_inches='tight')
print('Saved', outpath)
plt.close(fig)

# ── Leading-jet PVM |t⃗| comparison ──────────────────────────────────────────
spvm_bins        = np.linspace(0, 0.06, 21)
spvm_bin_centers = (spvm_bins[:-1] + spvm_bins[1:]) / 2
spvm_bin_width   = spvm_bins[1] - spvm_bins[0]

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'hLeadJetSPVM' not in rf:
                print(f'  hLeadJetSPVM not found in {fpath}')
                continue
            h         = rf['hLeadJetSPVM']
            counts    = h.values()
            variances = h.variances()
        norm = counts.sum() * spvm_bin_width
        if norm > 0:
            hist     = counts / norm
            hist_err = np.sqrt(variances) / norm
        else:
            hist     = np.zeros_like(counts, dtype=float)
            hist_err = np.zeros_like(counts, dtype=float)
        ax.hist(spvm_bin_centers, bins=spvm_bins, weights=hist, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
        ax.errorbar(spvm_bin_centers, hist, yerr=hist_err, fmt='none',
                    color=group_colors[s['group']], capsize=3, capthick=1, elinewidth=1)
    apply_style(ax,
                xlabel=r'$|\vec{t}\,|$',
                ylabel='Normalized density',
                title='',
                xlim=(0, 0.06),
                log_y=True,
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetSPVM_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Leading-jet PVM VBF/QCD ratio ────────────────────────────────────────────
for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for gen in ('mg5', 'herwig'):
        gen_samps = [s for s in samples if s['gen'] == gen]
        qcd_s = next((s for s in gen_samps if s['group'] == 'QCD'), None)
        vbf_s = next((s for s in gen_samps if s['group'] == 'VBF'), None)
        if qcd_s is None or vbf_s is None:
            continue
        qcd_counts = qcd_vars = vbf_counts = vbf_vars = None
        for s, store in [(qcd_s, 'qcd'), (vbf_s, 'vbf')]:
            fpath = DATA_DIR / f'{s["name"]}.friend.root'
            if not fpath.exists():
                print(f'  Skipping {fpath} (not found)')
                break
            with uproot.open(str(fpath)) as rf:
                if 'hLeadJetSPVM' not in rf:
                    print(f'  hLeadJetSPVM not found in {fpath}')
                    break
                h = rf['hLeadJetSPVM']
                if store == 'qcd':
                    qcd_counts, qcd_vars = h.values(), h.variances()
                else:
                    vbf_counts, vbf_vars = h.values(), h.variances()
        if qcd_counts is None or vbf_counts is None:
            continue
        qcd_norm = qcd_counts.sum() * spvm_bin_width
        vbf_norm = vbf_counts.sum() * spvm_bin_width
        if qcd_norm <= 0 or vbf_norm <= 0:
            continue
        qcd_hist = qcd_counts / qcd_norm
        vbf_hist = vbf_counts / vbf_norm
        qcd_err  = np.sqrt(qcd_vars) / qcd_norm
        vbf_err  = np.sqrt(vbf_vars) / vbf_norm
        mask = qcd_hist > 0
        ratio     = np.where(mask, vbf_hist / np.where(mask, qcd_hist, 1), np.nan)
        ratio_err = np.where(mask,
                             ratio * np.sqrt(
                                 (np.where(mask, vbf_err / np.where(vbf_hist > 0, vbf_hist, 1), 0))**2 +
                                 (np.where(mask, qcd_err / np.where(qcd_hist > 0, qcd_hist, 1), 0))**2),
                             np.nan)
        rms  = np.sqrt(np.nanmean(ratio**2))
        sep2 = separation_power(vbf_hist, qcd_hist, spvm_bin_width)
        kld  = kl_divergence(vbf_hist * spvm_bin_width, qcd_hist * spvm_bin_width)
        ax.hist(spvm_bin_centers, bins=spvm_bins, weights=np.where(np.isnan(ratio), 0, ratio),
                histtype='step', linewidth=3,
                color=gen_colors_ratio[gen], linestyle=gen_styles[gen],
                label=f'{gen_labels_ratio[gen]}  RMS={rms:.3f}  $\\langle S^2\\rangle$={sep2:.3f}  KL={kld:.3f}')
        ax.errorbar(spvm_bin_centers, ratio, yerr=ratio_err, fmt='none',
                    color=gen_colors_ratio[gen], capsize=3, capthick=1, elinewidth=1)
    ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')
    apply_style(ax,
                xlabel=r'$|\vec{t}\,|$',
                ylabel='VBF / QCD',
                title='',
                xlim=(0, 0.06),
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetSPVM_ratio_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Leading-jet NC (number of constituents) comparison ───────────────────────
nc_bins        = np.linspace(0, 100, 26)
nc_bin_centers = (nc_bins[:-1] + nc_bins[1:]) / 2
nc_bin_width   = nc_bins[1] - nc_bins[0]

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        counts = counts_from_tree(fpath, 'jetNC0', nc_bins)
        if counts is None:
            print(f'  jetNC not found in {fpath}')
            continue
        norm = counts.sum() * nc_bin_width
        if norm > 0:
            hist     = counts / norm
            hist_err = np.sqrt(counts) / norm
        else:
            hist     = np.zeros_like(counts, dtype=float)
            hist_err = np.zeros_like(counts, dtype=float)
        ax.hist(nc_bin_centers, bins=nc_bins, weights=hist, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
        ax.errorbar(nc_bin_centers, hist, yerr=hist_err, fmt='none',
                    color=group_colors[s['group']], capsize=3, capthick=1, elinewidth=1)
    apply_style(ax,
                xlabel=r'$N_\mathrm{constituents}$',
                ylabel='Normalized density',
                title='',
                xlim=(0, 100),
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetNC_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Leading-jet NC VBF/QCD ratio ─────────────────────────────────────────────
for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for gen in ('mg5', 'herwig'):
        gen_samps = [s for s in samples if s['gen'] == gen]
        qcd_s = next((s for s in gen_samps if s['group'] == 'QCD'), None)
        vbf_s = next((s for s in gen_samps if s['group'] == 'VBF'), None)
        if qcd_s is None or vbf_s is None:
            continue
        qcd_counts = vbf_counts = None
        for s, store in [(qcd_s, 'qcd'), (vbf_s, 'vbf')]:
            fpath = DATA_DIR / f'{s["name"]}.friend.root'
            if not fpath.exists():
                print(f'  Skipping {fpath} (not found)')
                break
            c = counts_from_tree(fpath, 'jetNC0', nc_bins)
            if c is None:
                print(f'  jetNC not found in {fpath}')
                break
            if store == 'qcd':
                qcd_counts = c
            else:
                vbf_counts = c
        if qcd_counts is None or vbf_counts is None:
            continue
        qcd_norm = qcd_counts.sum() * nc_bin_width
        vbf_norm = vbf_counts.sum() * nc_bin_width
        if qcd_norm <= 0 or vbf_norm <= 0:
            continue
        qcd_hist = qcd_counts / qcd_norm
        vbf_hist = vbf_counts / vbf_norm
        qcd_err  = np.sqrt(qcd_counts) / qcd_norm
        vbf_err  = np.sqrt(vbf_counts) / vbf_norm
        mask = qcd_hist > 0
        ratio     = np.where(mask, vbf_hist / np.where(mask, qcd_hist, 1), np.nan)
        ratio_err = np.where(mask,
                             ratio * np.sqrt(
                                 (np.where(mask, vbf_err / np.where(vbf_hist > 0, vbf_hist, 1), 0))**2 +
                                 (np.where(mask, qcd_err / np.where(qcd_hist > 0, qcd_hist, 1), 0))**2),
                             np.nan)
        rms  = np.sqrt(np.nanmean(ratio**2))
        sep2 = separation_power(vbf_hist, qcd_hist, nc_bin_width)
        kld  = kl_divergence(vbf_hist * nc_bin_width, qcd_hist * nc_bin_width)
        ax.hist(nc_bin_centers, bins=nc_bins, weights=np.where(np.isnan(ratio), 0, ratio),
                histtype='step', linewidth=3,
                color=gen_colors_ratio[gen], linestyle=gen_styles[gen],
                label=f'{gen_labels_ratio[gen]}  RMS={rms:.3f}  $\\langle S^2\\rangle$={sep2:.3f}  KL={kld:.3f}')
        ax.errorbar(nc_bin_centers, ratio, yerr=ratio_err, fmt='none',
                    color=gen_colors_ratio[gen], capsize=3, capthick=1, elinewidth=1)
    ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')
    apply_style(ax,
                xlabel=r'$N_\mathrm{constituents}$',
                ylabel='VBF / QCD',
                title='',
                xlim=(0, 100),
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetNC_ratio_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── ROC curves: |θs|, |vec{ts}|, mjj, ΔYjj ────────────────────────────────
def _hist_counts(hname):
    def _fn(fpath):
        with uproot.open(str(fpath)) as rf:
            if hname not in rf:
                return None
            return rf[hname].values().astype(float)
    return _fn

def _tree_counts(branch, bins, transform=None):
    def _fn(fpath):
        return counts_from_tree(fpath, branch, bins, transform=transform)
    return _fn

roc_vars = [
    (r'$|\theta_s|$',    '-',   _hist_counts('hLeadJetSPVA')),
    (r'$m_{jj}$',        '-.',  _tree_counts('mjj',  np.arange(0, 2010, 10))),
]
gen_colors_roc = {'mg5': 'black', 'herwig': 'royalblue'}
gen_labels_roc = {'mg5': 'Pythia', 'herwig': 'Herwig'}

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for var_label, ls, counts_fn in roc_vars:
        for gen in ('mg5', 'herwig'):
            gen_samps = [s for s in samples if s['gen'] == gen]
            qcd_s = next((s for s in gen_samps if s['group'] == 'QCD'), None)
            vbf_s = next((s for s in gen_samps if s['group'] == 'VBF'), None)
            if qcd_s is None or vbf_s is None:
                continue
            qcd_fpath = DATA_DIR / f'{qcd_s["name"]}.friend.root'
            vbf_fpath = DATA_DIR / f'{vbf_s["name"]}.friend.root'
            if not qcd_fpath.exists() or not vbf_fpath.exists():
                print(f'  Skipping missing files for {proc} {gen}')
                continue
            qcd_c = counts_fn(qcd_fpath)
            vbf_c = counts_fn(vbf_fpath)
            if qcd_c is None or vbf_c is None:
                continue
            fpr, tpr, auc = roc_from_histograms(vbf_c, qcd_c)
            if fpr is None:
                continue
            ax.plot(fpr, tpr, linewidth=3,
                    color=gen_colors_roc[gen], linestyle=ls,
                    label=f'{var_label} {gen_labels_roc[gen]}  AUC={auc:.3f}')
    ax.plot([0, 1], [0, 1], color='gray', linewidth=1, linestyle='--')
    apply_style(ax,
                xlabel='Background efficiency',
                ylabel='Signal efficiency',
                title='',
                xlim=(0, 1), ylim=(0, 1),
                legend_loc='lower right')
    plt.tight_layout()
    outpath = str(FIGS / f'jet_roc_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Summary ROC: all NNs per channel (no single-observable curves) ───────────
# Reads the roc_data_{proc}.npz files written by each NN dir's inference.py
# (solid = Herwig test split, dashed = MG5+Pythia transfer).
NN_SUMMARY = [
    ('NNkin',  '#7f7f7f'),
    ('NNj',    '#1f77b4'),
    ('NNjB',   '#2ca02c'),
    ('NNjj',   '#ff7f0e'),
    ('NNjjB',  '#9467bd'),
    ('NNjjBj', '#d62728'),
]
NN_BASE = pathlib.Path(__file__).parent

for proc in SAMPLES_BY_PROC:
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    n_curves = 0
    for nn, color in NN_SUMMARY:
        nn_dir = NN_BASE / (nn if proc == 'Z' else f'{nn}_{proc}')
        npz_path = nn_dir / f'roc_data_{proc}.npz'
        if not npz_path.exists():
            print(f'  roc_summary_{proc}: missing {npz_path.name} in {nn_dir.name}, skipping')
            continue
        d = np.load(npz_path)
        ax.plot(d['fpr_hw'],  d['tpr_hw'],  linewidth=2, color=color, linestyle='-',
                label=f'{nn} Herwig  AUC={float(d["auc_hw"]):.3f}')
        ax.plot(d['fpr_mg5'], d['tpr_mg5'], linewidth=2, color=color, linestyle='--',
                label=f'{nn} Pythia  AUC={float(d["auc_mg5"]):.3f}')
        n_curves += 2
    if n_curves == 0:
        plt.close(fig)
        continue
    ax.plot([0, 1], [0, 1], color='gray', linewidth=1, linestyle='--')
    apply_style(ax,
                xlabel='Background efficiency',
                ylabel='Signal efficiency',
                title='',
                xlim=(0, 1), ylim=(0, 1),
                legend_loc='lower right')
    ax.legend(loc='lower right', frameon=False, fontsize=10)  # 12 entries: shrink
    plt.tight_layout()
    outpath = str(FIGS / f'roc_summary_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── mjj distribution (absolute normalization via kWeight) ────────────────────
mjj_bins        = np.arange(0, 2010, 10)
mjj_bin_centers = (mjj_bins[:-1] + mjj_bins[1:]) / 2

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'events' not in rf:
                print(f'  events tree not found in {fpath}')
                continue
            mjj_vals = rf['events']['mjj'].array(library='np').astype(float)
            kw_vals  = rf['events']['kWeight'].array(library='np').astype(float)
        counts, _ = np.histogram(np.clip(mjj_vals, mjj_bins[0], mjj_bins[-1]),
                                 bins=mjj_bins, weights=kw_vals)
        ax.hist(mjj_bin_centers, bins=mjj_bins, weights=counts, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel=r'$m_{jj}$ [GeV]',
                ylabel=r'events / 10 GeV / pb$^{-1}$',
                title='',
                xlim=(0, 2000),
                log_y=True,
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'mjj_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── dYjj distribution (absolute normalization via kWeight) ──────────────────
dYjj_bins        = np.arange(-8, 8.25, 0.25)
dYjj_bin_centers = (dYjj_bins[:-1] + dYjj_bins[1:]) / 2

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'events' not in rf:
                print(f'  events tree not found in {fpath}')
                continue
            dYjj_vals = rf['events']['dYjj'].array(library='np').astype(float)
            kw_vals   = rf['events']['kWeight'].array(library='np').astype(float)
        counts, _ = np.histogram(np.clip(dYjj_vals, dYjj_bins[0], dYjj_bins[-1]),
                                 bins=dYjj_bins, weights=kw_vals)
        ax.hist(dYjj_bin_centers, bins=dYjj_bins, weights=counts, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel=r'$\Delta y_{jj}$',
                ylabel=r'events / 0.25 / pb$^{-1}$',
                title='',
                xlim=(-8, 8),
                log_y=True,
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'dYjj_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── dPhijj distribution (absolute normalization via kWeight) ────────────────
dPhijj_bins        = np.arange(-np.pi, np.pi + 0.05, 0.1)
dPhijj_bin_centers = (dPhijj_bins[:-1] + dPhijj_bins[1:]) / 2

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'events' not in rf:
                print(f'  events tree not found in {fpath}')
                continue
            dPhijj_vals = rf['events']['dPhijj'].array(library='np').astype(float)
            kw_vals     = rf['events']['kWeight'].array(library='np').astype(float)
        counts, _ = np.histogram(np.clip(dPhijj_vals, dPhijj_bins[0], dPhijj_bins[-1]),
                                 bins=dPhijj_bins, weights=kw_vals)
        ax.hist(dPhijj_bin_centers, bins=dPhijj_bins, weights=counts, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel=r'$\Delta\phi_{jj}$ [rad]',
                ylabel=r'events / 0.1 rad / pb$^{-1}$',
                title='',
                xlim=(-np.pi, np.pi),
                log_y=True,
                legend_loc='upper left')
    plt.tight_layout()
    outpath = str(FIGS / f'dPhijj_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Ptjj distribution (absolute normalization via kWeight) ──────────────────
Ptjj_bins        = np.arange(0, 305, 5)
Ptjj_bin_centers = (Ptjj_bins[:-1] + Ptjj_bins[1:]) / 2

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'events' not in rf:
                print(f'  events tree not found in {fpath}')
                continue
            Ptjj_vals = rf['events']['ptjj'].array(library='np').astype(float)
            kw_vals   = rf['events']['kWeight'].array(library='np').astype(float)
        counts, _ = np.histogram(np.clip(Ptjj_vals, Ptjj_bins[0], Ptjj_bins[-1]),
                                 bins=Ptjj_bins, weights=kw_vals)
        ax.hist(Ptjj_bin_centers, bins=Ptjj_bins, weights=counts, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel=r'$p_{\mathrm{T},jj}$ [GeV]',
                ylabel=r'events / 5 GeV / pb$^{-1}$',
                title='',
                xlim=(0, 300),
                log_y=True,
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'Ptjj_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Two-leading-jet pT (both jets superimposed, absolute norm) ───────────────
jetPt_bins        = np.arange(0, 155, 5)
jetPt_bin_centers = (jetPt_bins[:-1] + jetPt_bins[1:]) / 2

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'events' not in rf:
                print(f'  events tree not found in {fpath}')
                continue
            pt0 = rf['events']['jetPt0'].array(library='np').astype(float)
            pt1 = rf['events']['jetPt1'].array(library='np').astype(float)
            kw  = rf['events']['kWeight'].array(library='np').astype(float)
        pt_all = np.concatenate([pt0, pt1])
        kw_all = np.concatenate([kw, kw])
        counts, _ = np.histogram(np.clip(pt_all, jetPt_bins[0], jetPt_bins[-1]),
                                 bins=jetPt_bins, weights=kw_all)
        ax.hist(jetPt_bin_centers, bins=jetPt_bins, weights=counts, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel=r'$p_T^{\,j}$ [GeV]',
                ylabel=r'events / 5 GeV / pb$^{-1}$',
                title='',
                xlim=(0, 150),
                log_y=True,
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetPt_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── Subleading-jet pT (absolute normalization via kWeight) ───────────────────
for group in ['VBF', 'QCD']:
    for proc, samples in SAMPLES_BY_PROC.items():
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        for s in [s for s in samples if s['group'] == group]:
            fpath = DATA_DIR / f'{s["name"]}.friend.root'
            if not fpath.exists():
                print(f'  Skipping {fpath} (not found)')
                continue
            with uproot.open(str(fpath)) as rf:
                if 'events' not in rf:
                    print(f'  events tree not found in {fpath}')
                    continue
                pt1 = rf['events']['jetPt1'].array(library='np').astype(float)
                kw  = rf['events']['kWeight'].array(library='np').astype(float)
            counts, _ = np.histogram(np.clip(pt1, jetPt_bins[0], jetPt_bins[-1]),
                                     bins=jetPt_bins, weights=kw)
            ax.hist(jetPt_bin_centers, bins=jetPt_bins, weights=counts, histtype='step',
                    linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                    label=s['label'])
        apply_style(ax,
                    xlabel=r'$p_T^{\,j_1}$ [GeV]',
                    ylabel=r'events / 5 GeV / pb$^{-1}$',
                    title='',
                    xlim=(0, 150),
                    log_y=False,
                    legend_loc='upper right')
        plt.tight_layout()
        outpath = str(FIGS / f'hJetPt1_{group}_{proc}.pdf')
        fig.savefig(outpath, bbox_inches='tight')
        print('Saved', outpath)
        plt.close(fig)

# ── Two-leading-jet eta (both jets superimposed, absolute norm) ──────────────
jetEta_bins        = np.arange(-5, 5.1, 0.1)
jetEta_bin_centers = (jetEta_bins[:-1] + jetEta_bins[1:]) / 2

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'events' not in rf:
                print(f'  events tree not found in {fpath}')
                continue
            eta0 = rf['events']['jetEta0'].array(library='np').astype(float)
            eta1 = rf['events']['jetEta1'].array(library='np').astype(float)
            kw   = rf['events']['kWeight'].array(library='np').astype(float)
        eta_all = np.concatenate([eta0, eta1])
        kw_all  = np.concatenate([kw, kw])
        counts, _ = np.histogram(np.clip(eta_all, jetEta_bins[0], jetEta_bins[-1]),
                                 bins=jetEta_bins, weights=kw_all)
        ax.hist(jetEta_bin_centers, bins=jetEta_bins, weights=counts, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel=r'$\eta^{\,j}$',
                ylabel=r'events / 0.1 / pb$^{-1}$',
                title='',
                xlim=(-5, 5),
                log_y=True,
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'jetEta_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── hNJets (inclusive jet multiplicity before dijet selection) ────────────────
njets_bins        = np.arange(-0.5, 11.5, 1.0)
njets_bin_centers = np.arange(0, 11, dtype=float)
njets_bin_width   = 1.0

for group in ['VBF', 'QCD']:
    for proc, samples in SAMPLES_BY_PROC.items():
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        for s in [s for s in samples if s['group'] == group]:
            fpath = DATA_DIR / f'{s["name"]}.friend.root'
            if not fpath.exists():
                print(f'  Skipping {fpath} (not found)')
                continue
            with uproot.open(str(fpath)) as rf:
                if 'hNJets' not in rf:
                    print(f'  hNJets not found in {fpath}')
                    continue
                counts = rf['hNJets'].values().astype(float)
                n = len(njets_bin_centers)
                if len(counts) >= n:
                    counts = counts[:n]
                else:
                    counts = np.pad(counts, (0, n - len(counts)))
            ax.hist(njets_bin_centers, bins=njets_bins, weights=counts, histtype='step',
                    linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                    label=s['label'])
        apply_style(ax,
                    xlabel=r'$N_\mathrm{jets}$',
                    ylabel=r'events / 1 / pb$^{-1}$',
                    title='',
                    xlim=(-0.5, 10.5),
                    log_y=True,
                    legend_loc='upper right')
        plt.tight_layout()
        outpath = str(FIGS / f'hNJets_{group}_{proc}.pdf')
        fig.savefig(outpath, bbox_inches='tight')
        print('Saved', outpath)
        plt.close(fig)

# ── event-level NN scores: NNj_jet0 / NNjB / NNjjBj distributions ────────
score_bins        = np.linspace(0, 1, 51)
score_bin_centers = (score_bins[:-1] + score_bins[1:]) / 2
score_bin_width   = score_bins[1] - score_bins[0]

SCORE_COLS = [
    ('NNj_jet0', 'NNj score (jet 0)'),
    ('NNjB',     'NNjB score'),
    ('NNjjBj',   'NNjjBj score'),
]

# inclusive, same selection as mjj
for score_col, score_xlabel in SCORE_COLS:
  for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.h5'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with h5py.File(fpath, 'r') as f:
            if score_col not in f:
                print(f'  {score_col} not found in {fpath}')
                continue
            score = f[score_col][:].astype(float)
            kw    = f['kWeight'][:].astype(float)
        counts, _ = np.histogram(score, bins=score_bins, weights=kw)
        norm = counts.sum() * score_bin_width
        if norm <= 0:
            continue
        hist = counts / norm
        ax.hist(score_bin_centers, bins=score_bins, weights=hist, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel=score_xlabel,
                ylabel='Normalized density',
                title='',
                xlim=(0, 1),
                legend_loc='upper center')
    plt.tight_layout()
    outpath = str(FIGS / f'{score_col}_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# scores inclusive, absolute norm, log-y (same as mjj)
for score_col, score_xlabel in SCORE_COLS:
  for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.h5'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with h5py.File(fpath, 'r') as f:
            if score_col not in f:
                print(f'  {score_col} not found in {fpath}')
                continue
            score = f[score_col][:].astype(float)
            kw    = f['kWeight'][:].astype(float)
        counts, _ = np.histogram(score, bins=score_bins, weights=kw)
        ax.hist(score_bin_centers, bins=score_bins, weights=counts, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel=score_xlabel,
                ylabel=f'events / {score_bin_width:.2f} / pb$^{{-1}}$',
                title='',
                xlim=(0, 1),
                log_y=True,
                legend_loc='upper center')
    plt.tight_layout()
    outpath = str(FIGS / f'{score_col}_{proc}_abs.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── 2D score correlation: NNj_jet0 vs NNj_jet1, VBFZ_herwig vs QCDZjj_herwig ──
for name, label in [('VBFZ_herwig', 'VBF Z  Herwig'), ('QCDZjj_herwig', 'QCD Z  Herwig')]:
    fpath = DATA_DIR / f'{name}.h5'
    if not fpath.exists():
        print(f'  Skipping {fpath} (not found)')
        continue
    with h5py.File(fpath, 'r') as f:
        if 'NNj_jet0' not in f or 'NNj_jet1' not in f:
            print(f'  NNj_jet0/NNj_jet1 not found in {fpath}')
            continue
        jet0 = f['NNj_jet0'][:].astype(float)
        jet1 = f['NNj_jet1'][:].astype(float)
        kw   = f['kWeight'][:].astype(float)
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    _, _, _, im = ax.hist2d(jet0, jet1, bins=score_bins, weights=kw,
                             cmap='viridis', cmin=np.finfo(float).tiny)
    fig.colorbar(im, ax=ax, label='events / pb$^{-1}$')
    ax.set_xlabel('NNj score (jet 0)', fontsize=LABEL_SIZE)
    ax.set_ylabel('NNj score (jet 1)', fontsize=LABEL_SIZE)
    ax.set_title(label, fontsize=TITLE_SIZE)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', labelsize=TICK_SIZE)
    plt.tight_layout()
    outpath = str(FIGS / f'NNj_jet0_vs_jet1_{name}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── 1D NNj_jet0 shape sliced by NNj_jet1 tercile (low/mid/high), Herwig Z only ──
jet1_slices  = [('low', 0.0, 1/3), ('mid', 1/3, 2/3), ('high', 2/3, 1.0)]
slice_colors = {'low': '#440154', 'mid': '#21918c', 'high': '#fde725'}  # viridis triple

for name, label in [('VBFZ_herwig', 'VBF Z  Herwig'), ('QCDZjj_herwig', 'QCD Z  Herwig')]:
    fpath = DATA_DIR / f'{name}.h5'
    if not fpath.exists():
        print(f'  Skipping {fpath} (not found)')
        continue
    with h5py.File(fpath, 'r') as f:
        if 'NNj_jet0' not in f or 'NNj_jet1' not in f:
            print(f'  NNj_jet0/NNj_jet1 not found in {fpath}')
            continue
        jet0 = f['NNj_jet0'][:].astype(float)
        jet1 = f['NNj_jet1'][:].astype(float)
        kw   = f['kWeight'][:].astype(float)
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for slice_name, lo, hi in jet1_slices:
        mask = (jet1 >= lo) & (jet1 <= hi if hi == 1.0 else jet1 < hi)
        if mask.sum() == 0:
            continue
        counts, _ = np.histogram(jet0[mask], bins=score_bins, weights=kw[mask])
        norm = counts.sum() * score_bin_width
        if norm <= 0:
            continue
        hist = counts / norm
        ax.hist(score_bin_centers, bins=score_bins, weights=hist, histtype='step',
                linewidth=3, color=slice_colors[slice_name],
                label=f'NNj jet1 {slice_name}  [{lo:.2f},{hi:.2f}]')
    apply_style(ax,
                xlabel='NNj score (jet 0)',
                ylabel='Normalized density',
                title=label,
                xlim=(0, 1),
                legend_loc='upper center')
    plt.tight_layout()
    outpath = str(FIGS / f'NNj_jet0_slices_{name}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── update index.php ──────────────────────────────────────────────────────────
index_path = FIGS / 'index.php'
php = '''<?php
$pdfs = glob("*.pdf");
$txts = glob("*.txt");
$files = array_merge($pdfs ? $pdfs : array(), $txts ? $txts : array());
sort($files);
echo "<html><body>";
echo "<h2>akcolor plots</h2><ul>";
foreach ($files as $f) {
    echo "<li><a href=\\"$f\\">$f</a></li>";
}
echo "</ul></body></html>";
?>
'''
index_path.write_text(php)

subprocess.run([sys.executable, str(pathlib.Path(__file__).parent / 'checksum.py')], check=True)

# ── generate summer26.html from summer26.md ───────────────────────────────────
import make_html
make_html.main()
print('Updated', index_path)
