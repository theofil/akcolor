#!/usr/bin/env python3
import pathlib
import sys
import ROOT

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

def roc_from_histograms(sig_counts, bkg_counts):
    """ROC curve from histogram counts; returns (fpr, tpr, auc) scanning all thresholds."""
    sig_total = sig_counts.sum()
    bkg_total = bkg_counts.sum()
    if sig_total <= 0 or bkg_total <= 0:
        return None, None, None
    # fraction of events above each bin boundary, from high to low threshold
    sig_cdf = np.concatenate([np.cumsum(sig_counts[::-1])[::-1], [0]]) / sig_total
    bkg_cdf = np.concatenate([np.cumsum(bkg_counts[::-1])[::-1], [0]]) / bkg_total
    tpr = sig_cdf[::-1]   # 0 → 1 as threshold falls
    fpr = bkg_cdf[::-1]
    auc = float(np.trapz(tpr, fpr))
    if auc < 0.5:          # flip cut direction so AUC ≥ 0.5
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
DATA_DIR   = pathlib.Path('friends')
FIGS       = pathlib.Path.home() / 'www/files/akcolor'

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.LoadMacro(str(TDR_STYLE))
ROOT.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

FIGS.mkdir(parents=True, exist_ok=True)

# hLeadJetSPVA and hLeadJetSPVM are plotted separately below; skip individual files
SKIP_HISTS   = {'hLeadJetSPVA', 'hLeadJetSPVM'}
HIST_CLASSES = {'TH1F', 'TH2F', 'TH2Poly'}

for root_file in sorted(DATA_DIR.glob('*friend.root')):
    stem = root_file.stem
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

# ── Erase stale per-sample PDFs ──────────────────────────────────────────────
for stale in ('hLeadJetSPVA_QCDHtoInv.friend.pdf',
              'hLeadJetSPVA_VBFHtoInv.friend.pdf',
              'hLeadJetSPVA_comparison.pdf',
              'hLeadJetSPVA_comparison.png',
              'hLeadJetSPVM_QCDHtoInv.friend.pdf',
              'hLeadJetSPVM_VBFHtoInv.friend.pdf',
              'jetSPVA_a-1p0b2p0_QCDHtoInv.friend.pdf',
              'jetSPVA_a-1p0b2p0_VBFHtoInv.friend.pdf',
              'jetSPVA_a-2p0b2p0_QCDHtoInv.friend.pdf',
              'jetSPVA_a-2p0b2p0_VBFHtoInv.friend.pdf',
              'jetSPVA_a1p0b-1p0_QCDHtoInv.friend.pdf',
              'jetSPVA_a1p0b-1p0_VBFHtoInv.friend.pdf',
              'jetSPVA_a1p0b1p0_QCDHtoInv.friend.pdf',
              'jetSPVA_a1p0b1p0_VBFHtoInv.friend.pdf',
              'jetSPVA_a1p0b2p0_QCDHtoInv.friend.pdf',
              'jetSPVA_a1p0b2p0_VBFHtoInv.friend.pdf',
              'jetSPVA_a2p0b2p0_QCDHtoInv.friend.pdf',
              'jetSPVA_a2p0b2p0_VBFHtoInv.friend.pdf'):
    p = FIGS / stale
    if p.exists():
        p.unlink()
        print('Removed', p)

# ── Leading-jet |SPVA| comparison (matplotlib, matching spva.py style) ────────
spva_samples = [
    {'name': 'QCDHtoInv', 'label': 'QCD h', 'group': 'QCD', 'proc': 'H'},
    {'name': 'QCDZtoInv', 'label': 'QCD Z', 'group': 'QCD', 'proc': 'Z'},
    {'name': 'QCDWtoInv', 'label': 'QCD W', 'group': 'QCD', 'proc': 'W'},
    {'name': 'VBFHtoInv', 'label': 'VBF h', 'group': 'VBF', 'proc': 'H'},
    {'name': 'VBFZtoInv', 'label': 'VBF Z', 'group': 'VBF', 'proc': 'Z'},
    {'name': 'VBFWtoInv', 'label': 'VBF W', 'group': 'VBF', 'proc': 'W'},
]
proc_colors = {'H': 'black', 'Z': '#1f77b4', 'W': '#d62728'}
line_styles = {'QCD': '-', 'VBF': ':'}

bins        = np.linspace(0, np.pi, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width   = bins[1] - bins[0]

fig, ax = plt.subplots(figsize=FIG_SIZE)

for s in spva_samples:
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
            linewidth=3, color=proc_colors[s['proc']], linestyle=line_styles[s['group']],
            label=s['label'])
    ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='none',
                color=proc_colors[s['proc']], capsize=3, capthick=1, elinewidth=1)

apply_style(ax,
            xlabel=r'$|\theta_s|$ [rad]',
            ylabel='Normalized density [1/rad]',
            title='',
            xlim=(0, np.pi), ylim=(0.25, 0.47),
            legend_loc='upper right')
plt.tight_layout()
for ext in ('pdf',):
    outpath = str(FIGS / f'jetSPVA.{ext}')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
plt.close(fig)

# ── Leading-jet |SPVA| VBF/QCD ratio ─────────────────────────────────────────
procs = ['H', 'Z', 'W']
proc_names = {'H': ('QCDHtoInv', 'VBFHtoInv'),
              'Z': ('QCDZtoInv', 'VBFZtoInv'),
              'W': ('QCDWtoInv', 'VBFWtoInv')}
proc_labels = {'H': 'h', 'Z': 'Z', 'W': 'W'}

fig, ax = plt.subplots(figsize=FIG_SIZE)

for proc in procs:
    qcd_name, vbf_name = proc_names[proc]

    qcd_counts, qcd_vars = None, None
    vbf_counts, vbf_vars = None, None

    for sample_name, store in [(qcd_name, 'qcd'), (vbf_name, 'vbf')]:
        fpath = DATA_DIR / f'{sample_name}.friend.root'
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
                         ratio * np.sqrt((np.where(mask, vbf_err / np.where(mask, vbf_hist, 1), 0))**2 +
                                         (np.where(mask, qcd_err / np.where(mask, qcd_hist, 1), 0))**2),
                         np.nan)

    rms  = np.sqrt(np.nanmean(ratio**2))
    sep2 = separation_power(vbf_hist, qcd_hist, bin_width)
    kld  = kl_divergence(vbf_hist * bin_width, qcd_hist * bin_width)
    ax.hist(bin_centers, bins=bins, weights=np.where(np.isnan(ratio), 0, ratio),
            histtype='step', linewidth=3, color=proc_colors[proc],
            label=f'{proc_labels[proc]}  RMS={rms:.3f}  $\\langle S^2\\rangle$={sep2:.3f}  KL={kld:.3f}')
    ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt='none',
                color=proc_colors[proc], capsize=3, capthick=1, elinewidth=1)

ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')

apply_style(ax,
            xlabel=r'$|\theta_s|$ [rad]',
            ylabel='VBF / QCD',
            title='',
            xlim=(0, np.pi), ylim=(0.4, 1.6),
            legend_loc='upper right')
plt.tight_layout()
for ext in ('pdf',):
    outpath = str(FIGS / f'jetSPVA_ratio.{ext}')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
plt.close(fig)

# ── Leading-jet PVM |t⃗| comparison (matplotlib, log-y, matching pvm.py style) ─
spvm_bins        = np.linspace(0, 0.06, 21)
spvm_bin_centers = (spvm_bins[:-1] + spvm_bins[1:]) / 2
spvm_bin_width   = spvm_bins[1] - spvm_bins[0]

fig, ax = plt.subplots(figsize=FIG_SIZE)

for s in spva_samples:
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
            linewidth=3, color=proc_colors[s['proc']], linestyle=line_styles[s['group']],
            label=s['label'])
    ax.errorbar(spvm_bin_centers, hist, yerr=hist_err, fmt='none',
                color=proc_colors[s['proc']], capsize=3, capthick=1, elinewidth=1)

apply_style(ax,
            xlabel=r'$|\vec{t}\,|$',
            ylabel='Normalized density',
            title='',
            xlim=(0, 0.06),
            log_y=True,
            legend_loc='upper right')
plt.tight_layout()
for ext in ('pdf',):
    outpath = str(FIGS / f'jetSPVM.{ext}')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
plt.close(fig)

# ── Leading-jet PVM VBF/QCD ratio ────────────────────────────────────────────
fig, ax = plt.subplots(figsize=FIG_SIZE)

for proc in procs:
    qcd_name, vbf_name = proc_names[proc]

    qcd_counts, qcd_vars = None, None
    vbf_counts, vbf_vars = None, None

    for sample_name, store in [(qcd_name, 'qcd'), (vbf_name, 'vbf')]:
        fpath = DATA_DIR / f'{sample_name}.friend.root'
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
                         ratio * np.sqrt((np.where(mask, vbf_err / np.where(mask, vbf_hist, 1), 0))**2 +
                                         (np.where(mask, qcd_err / np.where(mask, qcd_hist, 1), 0))**2),
                         np.nan)

    rms  = np.sqrt(np.nanmean(ratio**2))
    sep2 = separation_power(vbf_hist, qcd_hist, spvm_bin_width)
    kld  = kl_divergence(vbf_hist * spvm_bin_width, qcd_hist * spvm_bin_width)
    ax.hist(spvm_bin_centers, bins=spvm_bins, weights=np.where(np.isnan(ratio), 0, ratio),
            histtype='step', linewidth=3, color=proc_colors[proc],
            label=f'{proc_labels[proc]}  RMS={rms:.3f}  $\\langle S^2\\rangle$={sep2:.3f}  KL={kld:.3f}')
    ax.errorbar(spvm_bin_centers, ratio, yerr=ratio_err, fmt='none',
                color=proc_colors[proc], capsize=3, capthick=1, elinewidth=1)

ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')

apply_style(ax,
            xlabel=r'$|\vec{t}\,|$',
            ylabel='VBF / QCD',
            title='',
            xlim=(0, 0.06),
            legend_loc='upper right')
plt.tight_layout()
for ext in ('pdf',):
    outpath = str(FIGS / f'jetSPVM_ratio.{ext}')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
plt.close(fig)

# ── Leading-jet NC (number of constituents) comparison ───────────────────────
nc_bins        = np.linspace(0, 100, 26)
nc_bin_centers = (nc_bins[:-1] + nc_bins[1:]) / 2
nc_bin_width   = nc_bins[1] - nc_bins[0]

fig, ax = plt.subplots(figsize=FIG_SIZE)

for s in spva_samples:
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
            linewidth=3, color=proc_colors[s['proc']], linestyle=line_styles[s['group']],
            label=s['label'])
    ax.errorbar(nc_bin_centers, hist, yerr=hist_err, fmt='none',
                color=proc_colors[s['proc']], capsize=3, capthick=1, elinewidth=1)

apply_style(ax,
            xlabel=r'$N_\mathrm{constituents}$',
            ylabel='Normalized density',
            title='',
            xlim=(0, 100),
            legend_loc='upper right')
plt.tight_layout()
for ext in ('pdf',):
    outpath = str(FIGS / f'jetNC.{ext}')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
plt.close(fig)

# ── Leading-jet NC VBF/QCD ratio ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=FIG_SIZE)

for proc in procs:
    qcd_name, vbf_name = proc_names[proc]

    qcd_counts, vbf_counts = None, None

    for sample_name, store in [(qcd_name, 'qcd'), (vbf_name, 'vbf')]:
        fpath = DATA_DIR / f'{sample_name}.friend.root'
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
                         ratio * np.sqrt((np.where(mask, vbf_err / np.where(mask, vbf_hist, 1), 0))**2 +
                                         (np.where(mask, qcd_err / np.where(mask, qcd_hist, 1), 0))**2),
                         np.nan)

    rms  = np.sqrt(np.nanmean(ratio**2))
    sep2 = separation_power(vbf_hist, qcd_hist, nc_bin_width)
    kld  = kl_divergence(vbf_hist * nc_bin_width, qcd_hist * nc_bin_width)
    ax.hist(nc_bin_centers, bins=nc_bins, weights=np.where(np.isnan(ratio), 0, ratio),
            histtype='step', linewidth=3, color=proc_colors[proc],
            label=f'{proc_labels[proc]}  RMS={rms:.3f}  $\\langle S^2\\rangle$={sep2:.3f}  KL={kld:.3f}')
    ax.errorbar(nc_bin_centers, ratio, yerr=ratio_err, fmt='none',
                color=proc_colors[proc], capsize=3, capthick=1, elinewidth=1)

ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')

apply_style(ax,
            xlabel=r'$N_\mathrm{constituents}$',
            ylabel='VBF / QCD',
            title='',
            xlim=(0, 100),
            legend_loc='upper right')
plt.tight_layout()
for ext in ('pdf',):
    outpath = str(FIGS / f'jetNC_ratio.{ext}')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
plt.close(fig)

# ── ROC curves: SPVA, PVM, jetPt, |jetEta|, jetM combined ───────────────────
for stale_roc in ('jetSPVA_roc.pdf', 'jetSPVA_roc.png', 'jetSPVM_roc.pdf', 'jetSPVM_roc.png'):
    p = FIGS / stale_roc
    if p.exists():
        p.unlink()
        print('Removed', p)

# (label, linestyle, counts_fn)  counts_fn(fpath) → np.ndarray or None
pt_bins  = np.linspace(30,  500, 51)
eta_bins = np.linspace(0,     3, 31)
m_bins   = np.linspace(0,   150, 51)

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

nc_roc_bins = np.linspace(0, 100, 26)

roc_vars = [
    (r'$|\theta_s|$',  '-',                      _hist_counts('hLeadJetSPVA')),
    (r'$|\vec{t}\,|$', '--',                      _hist_counts('hLeadJetSPVM')),
    (r'$p_T$',         '-.',                      _tree_counts('jetPt0',  pt_bins)),
    (r'$|\eta|$',      ':',                       _tree_counts('jetEta0', eta_bins, np.abs)),
    (r'$m$',           (0, (3, 1, 1, 1, 1, 1)),  _tree_counts('jetM0',   m_bins)),
    (r'$N_c$',         (0, (5, 1)),               _tree_counts('jetNC0',  nc_roc_bins)),
]

fig, ax = plt.subplots(figsize=FIG_SIZE)

for var_label, ls, counts_fn in roc_vars:
    for proc in procs:
        qcd_name, vbf_name = proc_names[proc]
        qcd_counts = vbf_counts = None
        for sample_name, store in [(qcd_name, 'qcd'), (vbf_name, 'vbf')]:
            fpath = DATA_DIR / f'{sample_name}.friend.root'
            if not fpath.exists():
                print(f'  Skipping {fpath} (not found)')
                break
            c = counts_fn(fpath)
            if c is None:
                print(f'  variable not found in {fpath}')
                break
            if store == 'qcd':
                qcd_counts = c
            else:
                vbf_counts = c
        if qcd_counts is None or vbf_counts is None:
            continue
        fpr, tpr, auc = roc_from_histograms(vbf_counts, qcd_counts)
        if fpr is None:
            continue
        plabel = proc_labels[proc]
        lbl = f'{var_label} {plabel}  AUC={auc:.3f}' if plabel != 'h' else f'{var_label}  AUC={auc:.3f}'
        ax.plot(fpr, tpr, linewidth=3, color=proc_colors[proc], linestyle=ls, label=lbl)

ax.plot([0, 1], [0, 1], color='gray', linewidth=1, linestyle='--')
apply_style(ax,
            xlabel='Background efficiency',
            ylabel='Signal efficiency',
            title='',
            xlim=(0, 1), ylim=(0, 1),
            legend_loc='lower right')
plt.tight_layout()
for ext in ('pdf',):
    outpath = str(FIGS / f'jet_roc.{ext}')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
plt.close(fig)

# ── TH2Poly VBF/QCD ratio2D (ROOT canvas, per process) ───────────────────────
# Discover all TH2Poly histogram names present in any slimfriend file
th2poly_names = set()
for _rf in sorted(DATA_DIR.glob('*friend.root')):
    _f = ROOT.TFile(str(_rf))
    if _f and not _f.IsZombie():
        for _key in _f.GetListOfKeys():
            if _key.GetClassName() == 'TH2Poly':
                th2poly_names.add(_key.GetName())
        _f.Close()

for _hname in sorted(th2poly_names):
    for proc in procs:
        qcd_name, vbf_name = proc_names[proc]
        qcd_fpath = DATA_DIR / f'{qcd_name}.friend.root'
        vbf_fpath = DATA_DIR / f'{vbf_name}.friend.root'
        if not qcd_fpath.exists() or not vbf_fpath.exists():
            continue

        qcd_f = ROOT.TFile(str(qcd_fpath))
        vbf_f = ROOT.TFile(str(vbf_fpath))
        if not qcd_f or qcd_f.IsZombie() or not vbf_f or vbf_f.IsZombie():
            continue

        h_qcd = qcd_f.Get(_hname)
        h_vbf = vbf_f.Get(_hname)
        if not h_qcd or not h_vbf:
            qcd_f.Close(); vbf_f.Close()
            continue

        n_bins = h_qcd.GetNumberOfBins()
        qcd_c = np.array([h_qcd.GetBinContent(i) for i in range(1, n_bins + 1)])
        vbf_c = np.array([h_vbf.GetBinContent(i) for i in range(1, n_bins + 1)])

        qcd_tot = qcd_c.sum()
        vbf_tot = vbf_c.sum()
        if qcd_tot <= 0 or vbf_tot <= 0:
            qcd_f.Close(); vbf_f.Close()
            continue

        # probability masses (sum = 1); area cancels in ratio, S^2, KL
        qcd_prob = qcd_c / qcd_tot
        vbf_prob = vbf_c / vbf_tot

        mask2d    = qcd_prob > 0
        ratio2d   = np.where(mask2d, vbf_prob / np.where(mask2d, qcd_prob, 1.0), np.nan)

        rms2d  = float(np.sqrt(np.nanmean(ratio2d**2)))
        sep2d  = separation_power(vbf_prob, qcd_prob, 1.0)
        kld    = kl_divergence(vbf_prob, qcd_prob)

        # build ratio TH2Poly using QCD geometry as template
        h_ratio = h_qcd.Clone(f'h_ratio2D_{_hname}_{proc}')
        h_ratio.SetDirectory(0)
        h_ratio.Reset("")
        h_ratio.SetTitle('')
        for _i, _rv in enumerate(ratio2d, start=1):
            if not np.isnan(_rv):
                h_ratio.SetBinContent(_i, float(_rv))

        c2d = ROOT.TCanvas(f'cratio2D_{_hname}_{proc}', '', 700, 600)
        c2d.SetRightMargin(0.18)
        c2d.SetTopMargin(0.12)
        h_ratio.GetXaxis().SetTitle('#Delta y_{i}')
        h_ratio.GetYaxis().SetTitle('#Delta#varphi_{i}')
        h_ratio.GetZaxis().SetTitle('VBF / QCD')
        h_ratio.GetXaxis().SetNdivisions(505)
        h_ratio.GetYaxis().SetNdivisions(505)
        h_ratio.Draw('COLZ')

        lat2d = ROOT.TLatex()
        lat2d.SetNDC()
        lat2d.SetTextSize(0.030)
        lat2d.SetTextFont(42)
        stats_line = (f'RMS={rms2d:.3f}   '
                      f'#langle S^{{2}}#rangle={sep2d:.4f}   '
                      f'D_{{KL}}={kld:.4f}')
        lat2d.DrawLatex(0.13, 0.94, stats_line)

        outpath2d = str(FIGS / f'{_hname}_ratio2D_{proc_labels[proc]}.pdf')
        c2d.SaveAs(outpath2d)
        c2d.Close()
        print('  Saved', outpath2d)

        qcd_f.Close()
        vbf_f.Close()

# ── update index.php ──────────────────────────────────────────────────────────
index_path = FIGS / 'index.php'
php = '''<?php
$files = glob("*.pdf");
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
print('Updated', index_path)
