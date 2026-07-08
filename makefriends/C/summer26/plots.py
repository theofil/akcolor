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

SKIP_HISTS   = {'hLeadJetSPVA', 'hLeadJetSPVM', 'hJetPt1', 'hNJets'}
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
for pattern in ('hJetPt1_*.pdf', 'hNJets_*.pdf'):
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
    'jetSPVA_SR_Run3H.pdf',
    'jetSPVA_SR_Run3_W.pdf', 'jetSPVA_SR_Run3_Z.pdf', 'jetSPVA_SR_Run3_H.pdf',
):
    p = FIGS / stale_name
    if p.exists():
        p.unlink()
        print('Removed', p)

# ── Sample catalog ─────────────────────────────────────────────────────────────
SAMPLES_BY_PROC = {
    'H': [
        {'name': 'QCDHjj_mg5_pythia', 'label': 'QCD H  MG5+Pythia', 'group': 'QCD', 'gen': 'mg5'},
        {'name': 'QCDHjj_herwig',     'label': 'QCD H  Herwig',      'group': 'QCD', 'gen': 'herwig'},
        {'name': 'VBFH_mg5_pythia',   'label': 'VBF H  MG5+Pythia',  'group': 'VBF', 'gen': 'mg5'},
        {'name': 'VBFH_herwig',       'label': 'VBF H  Herwig',      'group': 'VBF', 'gen': 'herwig'},
    ],
    'W': [
        {'name': 'QCDWjj_mg5_pythia', 'label': 'QCD W  MG5+Pythia', 'group': 'QCD', 'gen': 'mg5'},
        {'name': 'QCDWjj_herwig',     'label': 'QCD W  Herwig',      'group': 'QCD', 'gen': 'herwig'},
        {'name': 'VBFW_mg5_pythia',   'label': 'VBF W  MG5+Pythia',  'group': 'VBF', 'gen': 'mg5'},
        {'name': 'VBFW_herwig',       'label': 'VBF W  Herwig',      'group': 'VBF', 'gen': 'herwig'},
    ],
    'Z': [
        {'name': 'QCDZjj_mg5_pythia', 'label': 'QCD Z  MG5+Pythia', 'group': 'QCD', 'gen': 'mg5'},
        {'name': 'QCDZjj_herwig',     'label': 'QCD Z  Herwig',      'group': 'QCD', 'gen': 'herwig'},
        {'name': 'VBFZ_mg5_pythia',   'label': 'VBF Z  MG5+Pythia',  'group': 'VBF', 'gen': 'mg5'},
        {'name': 'VBFZ_herwig',       'label': 'VBF Z  Herwig',      'group': 'VBF', 'gen': 'herwig'},
    ],
}
# color by group (VBF/QCD), linestyle by generator
group_colors     = {'VBF': 'steelblue', 'QCD': 'tomato'}
gen_styles       = {'mg5': '-', 'herwig': '--'}
# for ratio plots: color by generator
gen_colors_ratio = {'mg5': 'black', 'herwig': 'royalblue'}
gen_labels_ratio = {'mg5': 'MG5+Pythia', 'herwig': 'Herwig'}

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
gen_labels_roc = {'mg5': 'MG5', 'herwig': 'Herwig'}

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
                label=f'{nn} MG5+Py  AUC={float(d["auc_mg5"]):.3f}')
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

# ── mjj distribution Run3 (L = 300 000 pb⁻¹ = 300 fb⁻¹) ────────────────────
L_RUN3 = 300_000.0  # pb^{-1}

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
            mjj_vals  = rf['events']['mjj'].array(library='np').astype(float)
            kw_vals   = rf['events']['kWeight'].array(library='np').astype(float)
            dyjj_vals = rf['events']['dYjj'].array(library='np').astype(float)
        sr = np.abs(dyjj_vals) > 3.0
        mjj_vals, kw_vals = mjj_vals[sr], kw_vals[sr]
        counts_pb, _ = np.histogram(np.clip(mjj_vals, mjj_bins[0], mjj_bins[-1]),
                                    bins=mjj_bins, weights=kw_vals)
        counts_run3 = counts_pb * L_RUN3
        yerr    = np.sqrt(np.maximum(counts_run3, 0))
        yerr_lo = np.minimum(yerr, counts_run3)
        color  = group_colors[s['group']]
        lstyle = gen_styles[s['gen']]
        ax.hist(mjj_bin_centers, bins=mjj_bins, weights=counts_run3,
                histtype='step', linewidth=3,
                color=color, linestyle=lstyle, label=s['label'])
        ax.errorbar(mjj_bin_centers, counts_run3,
                    yerr=[yerr_lo, yerr],
                    fmt='none', ecolor=color, elinewidth=1, capsize=2, alpha=0.6)
    apply_style(ax,
                xlabel=r'$m_{jj}$ [GeV]',
                ylabel=r'events / 10 GeV / 300 fb$^{-1}$',
                title=r'SR: $|\Delta Y_{jj}| > 3$',
                xlim=(0, 2000),
                log_y=True,
                legend_loc='upper right')
    plt.tight_layout()
    outpath = str(FIGS / f'mjj_{proc}_Run3.pdf')
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

# ── NNj score jet0 (inclusive, same selection as mjj) ───────────────────
score_bins        = np.linspace(0, 1, 51)
score_bin_centers = (score_bins[:-1] + score_bins[1:]) / 2
score_bin_width   = score_bins[1] - score_bins[0]

for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.h5'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with h5py.File(fpath, 'r') as f:
            if 'NNj_jet0' not in f:
                print(f'  NNj_jet0 not found in {fpath}')
                continue
            score = f['NNj_jet0'][:].astype(float)
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
                xlabel='NNj score (jet 0)',
                ylabel='Normalized density',
                title='',
                xlim=(0, 1),
                legend_loc='upper center')
    plt.tight_layout()
    outpath = str(FIGS / f'NNj_jet0_{proc}.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── NNj score jet0 in SR (|dYjj| > 3, same selection as mjj_Run3) ───────
for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.h5'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with h5py.File(fpath, 'r') as f:
            if 'NNj_jet0' not in f:
                print(f'  NNj_jet0 not found in {fpath}')
                continue
            score = f['NNj_jet0'][:].astype(float)
            kw    = f['kWeight'][:].astype(float)
            dYjj  = f['dYjj'][:].astype(float)
        sr = np.abs(dYjj) > 3.0
        score, kw = score[sr], kw[sr]
        counts, _ = np.histogram(score, bins=score_bins, weights=kw)
        norm = counts.sum() * score_bin_width
        if norm <= 0:
            continue
        hist = counts / norm
        ax.hist(score_bin_centers, bins=score_bins, weights=hist, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel='NNj score (jet 0)',
                ylabel='Normalized density',
                title=r'SR: $|\Delta Y_{jj}| > 3$',
                xlim=(0, 1),
                legend_loc='upper center')
    plt.tight_layout()
    outpath = str(FIGS / f'NNj_jet0_{proc}_Run3.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── NNj score jet0 inclusive, absolute norm, log-y (same as mjj) ────────
for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.h5'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with h5py.File(fpath, 'r') as f:
            if 'NNj_jet0' not in f:
                print(f'  NNj_jet0 not found in {fpath}')
                continue
            score = f['NNj_jet0'][:].astype(float)
            kw    = f['kWeight'][:].astype(float)
        counts, _ = np.histogram(score, bins=score_bins, weights=kw)
        ax.hist(score_bin_centers, bins=score_bins, weights=counts, histtype='step',
                linewidth=3, color=group_colors[s['group']], linestyle=gen_styles[s['gen']],
                label=s['label'])
    apply_style(ax,
                xlabel='NNj score (jet 0)',
                ylabel=f'events / {score_bin_width:.2f} / pb$^{{-1}}$',
                title='',
                xlim=(0, 1),
                log_y=True,
                legend_loc='upper center')
    plt.tight_layout()
    outpath = str(FIGS / f'NNj_jet0_{proc}_abs.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── NNj score jet0 SR, absolute norm × L_RUN3, log-y (same as mjj_Run3) ─
for proc, samples in SAMPLES_BY_PROC.items():
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.h5'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with h5py.File(fpath, 'r') as f:
            if 'NNj_jet0' not in f:
                print(f'  NNj_jet0 not found in {fpath}')
                continue
            score = f['NNj_jet0'][:].astype(float)
            kw    = f['kWeight'][:].astype(float)
            dYjj  = f['dYjj'][:].astype(float)
        sr = np.abs(dYjj) > 3.0
        score, kw = score[sr], kw[sr]
        counts_pb, _ = np.histogram(score, bins=score_bins, weights=kw)
        counts_run3  = counts_pb * L_RUN3
        yerr    = np.sqrt(np.maximum(counts_run3, 0))
        yerr_lo = np.minimum(yerr, counts_run3)
        color  = group_colors[s['group']]
        lstyle = gen_styles[s['gen']]
        ax.hist(score_bin_centers, bins=score_bins, weights=counts_run3,
                histtype='step', linewidth=3,
                color=color, linestyle=lstyle, label=s['label'])
        ax.errorbar(score_bin_centers, counts_run3,
                    yerr=[yerr_lo, yerr],
                    fmt='none', ecolor=color, elinewidth=1, capsize=2, alpha=0.6)
    apply_style(ax,
                xlabel='NNj score (jet 0)',
                ylabel=f'events / {score_bin_width:.2f} / 300 fb$^{{-1}}$',
                title=r'SR: $|\Delta Y_{jj}| > 3$',
                xlim=(0, 1),
                log_y=True,
                legend_loc='upper center')
    plt.tight_layout()
    outpath = str(FIGS / f'NNj_jet0_{proc}_Run3_abs.pdf')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)

# ── SR_Run3 optimization: per-channel (mjj, |dYjj|, NNj_jet0) cuts ───────
# maximizing S/(S+B). Scan on Herwig samples only, L = 300 fb⁻¹, subject to
# S > 1000 expected signal events and ≥ 10 raw MC background events (guard
# against the B→0 MC-stat floor).
sr_mjj_edges   = np.arange(0, 4010.0, 10.0)
sr_dyjj_edges  = np.arange(0, 6.1, 0.1)
sr_score_edges = np.arange(0, 1.01, 0.01)

def _sr_cumulative_yields(name):
    """Reverse-cumulative (kWeight×L, raw count) grids over
    (mjj, |dYjj|, NNj_jet0) cuts."""
    with h5py.File(DATA_DIR / f'{name}.h5', 'r') as f:
        mjj   = f['mjj'][:].astype(float)
        dY    = np.abs(f['dYjj'][:].astype(float))
        score = f['NNj_jet0'][:].astype(float)
        kw    = f['kWeight'][:].astype(float)
    edges  = [sr_mjj_edges, sr_dyjj_edges, sr_score_edges]
    sample = np.stack([np.clip(v, e[0], e[-1] - 1e-9)
                       for v, e in zip((mjj, dY, score), edges)], axis=1)
    hw, _ = np.histogramdd(sample, bins=edges, weights=kw)
    hn, _ = np.histogramdd(sample, bins=edges)
    def revcum(h):
        for ax in range(h.ndim):
            h = np.flip(np.cumsum(np.flip(h, axis=ax), axis=ax), axis=ax)
        return h
    return revcum(hw) * L_RUN3, revcum(hn)

SR_CUTS = {}
sr_table_lines = [
    'SR_Run3 optimization (Herwig only, L = 300 fb^-1)',
    'maximize S/(S+B) subject to S > 1000 and >= 10 raw MC background events',
    '',
    f'{"chan":>4} {"mjj cut":>8} {"|dYjj| cut":>10} {"NN cut":>7} {"S":>10} '
    f'{"B":>10} {"S/(S+B)":>8} {"raw S":>8} {"raw B":>8}',
]
for proc, samples in SAMPLES_BY_PROC.items():
    hw_samps = [s for s in samples if s['gen'] == 'herwig']
    sig_name = next(s['name'] for s in hw_samps if s['group'] == 'VBF')
    bkg_name = next(s['name'] for s in hw_samps if s['group'] == 'QCD')
    S, rawS = _sr_cumulative_yields(sig_name)
    B, rawB = _sr_cumulative_yields(bkg_name)
    purity = S / np.maximum(S + B, 1e-12)
    allowed = (S > 1000.0) & (rawB >= 10)
    i, j, k = np.unravel_index(np.argmax(np.where(allowed, purity, -1.0)),
                               purity.shape)
    SR_CUTS[proc] = (float(sr_mjj_edges[i]), float(sr_dyjj_edges[j]),
                     float(sr_score_edges[k]))
    sr_table_lines.append(
        f'{proc:>4} {sr_mjj_edges[i]:>8.0f} {sr_dyjj_edges[j]:>10.1f} '
        f'{sr_score_edges[k]:>7.2f} {S[i, j, k]:>10.1f} {B[i, j, k]:>10.1f} '
        f'{purity[i, j, k]:>8.3f} {rawS[i, j, k]:>8.0f} {rawB[i, j, k]:>8.0f}')
sr_table = '\n'.join(sr_table_lines) + '\n'
print(sr_table)
(FIGS / 'SR_optimization.txt').write_text(sr_table)
print('Saved', FIGS / 'SR_optimization.txt')

# ── jetSPVA in SR_Run3 (optimized cuts), per generator, absolute norm, 300 fb⁻¹ ─
spva_bins        = np.linspace(0, np.pi, 21)
spva_bin_centers = (spva_bins[:-1] + spva_bins[1:]) / 2

GEN_LABEL = {'mg5': 'MG5Pythia', 'herwig': 'Herwig'}
GEN_TITLE = {'mg5': 'MG5+Pythia', 'herwig': 'Herwig'}

def _spva_hist(name, jet, mjj_cut, dyjj_cut, score_cut):
    """Return (hist, mc_sq) in SR_Run3 for jet-`jet` |θs|, scaled to L_RUN3."""
    fpath = DATA_DIR / f'{name}.h5'
    if not fpath.exists():
        print(f'  Skipping {fpath} (not found)')
        return None, None
    with h5py.File(fpath, 'r') as f:
        spva  = np.abs(f['jetSPVA'][:, jet].astype(float))
        kw    = f['kWeight'][:].astype(float)
        dYjj  = f['dYjj'][:].astype(float)
        mjj   = f['mjj'][:].astype(float)
        score = f['NNj_jet0'][:].astype(float)
    sr  = (np.abs(dYjj) > dyjj_cut) & (mjj > mjj_cut) & (score > score_cut)
    w   = kw[sr] * L_RUN3
    h,  _ = np.histogram(spva[sr], bins=spva_bins, weights=w)
    hsq,_ = np.histogram(spva[sr], bins=spva_bins, weights=w**2)
    return h, hsq

def _draw_hist_with_errs(ax, centers, bins, hist, mc_sq, color, label):
    ax.hist(centers, bins=bins, weights=hist,
            histtype='step', linewidth=3, color=color, label=label)
    data_err = np.sqrt(np.maximum(hist, 0))
    mc_err   = np.sqrt(mc_sq)
    # MC stat: thick translucent bar, no caps
    ax.errorbar(centers, hist,
                yerr=[np.minimum(mc_err, hist), mc_err],
                fmt='none', ecolor=color, elinewidth=5, capsize=0, alpha=0.25)
    # Data stat: thin bar with caps
    ax.errorbar(centers, hist,
                yerr=[np.minimum(data_err, hist), data_err],
                fmt='none', ecolor=color, elinewidth=1.5, capsize=3, alpha=0.8)

for gen in ('mg5', 'herwig'):
    for proc, samples in SAMPLES_BY_PROC.items():
        gen_samps = [s for s in samples if s['gen'] == gen]
        sig = next(s for s in gen_samps if s['group'] == 'VBF')
        bkg = next(s for s in gen_samps if s['group'] == 'QCD')
        mjj_cut, dyjj_cut, score_cut = SR_CUTS[proc]

        for jet in (0, 1):
            sig_hist, sig_mc_sq = _spva_hist(sig['name'], jet, mjj_cut, dyjj_cut, score_cut)
            bkg_hist, bkg_mc_sq = _spva_hist(bkg['name'], jet, mjj_cut, dyjj_cut, score_cut)
            if sig_hist is None or bkg_hist is None:
                continue

            splusb_hist  = sig_hist  + bkg_hist
            splusb_mc_sq = sig_mc_sq + bkg_mc_sq

            fig, ax = plt.subplots(figsize=FIG_SIZE)
            _draw_hist_with_errs(ax, spva_bin_centers, spva_bins,
                                 sig_hist, sig_mc_sq, 'steelblue',
                                 f'VBF {proc}  ({int(round(sig_hist.sum()))} events)')
            _draw_hist_with_errs(ax, spva_bin_centers, spva_bins,
                                 bkg_hist, bkg_mc_sq, 'tomato',
                                 f'QCD {proc}  ({int(round(bkg_hist.sum()))} events)')
            _draw_hist_with_errs(ax, spva_bin_centers, spva_bins,
                                 splusb_hist, splusb_mc_sq, 'black',
                                 f'S+B  ({int(round(splusb_hist.sum()))} events)')
            apply_style(ax,
                        xlabel=rf'$|\theta_s|$ (jet {jet}) [rad]',
                        ylabel=r'events / bin / 300 fb$^{-1}$',
                        title=(r'SR$_{\mathrm{Run3}}$ (optimized): '
                               rf'$|\Delta Y_{{jj}}|>{dyjj_cut:.1f}$, '
                               rf'$m_{{jj}}>{mjj_cut / 1000:.2f}$ TeV, '
                               rf'NN$>{score_cut:.2f}$'
                               f'  [{GEN_TITLE[gen]}]'),
                        xlim=(0, np.pi),
                        legend_loc='upper right')
            plt.tight_layout()
            jet_tag = '' if jet == 0 else f'_jet{jet}'
            outpath = str(FIGS / f'jetSPVA{jet_tag}_SR_Run3_{proc}_{GEN_LABEL[gen]}.pdf')
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
