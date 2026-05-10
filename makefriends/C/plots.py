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

for root_file in sorted(DATA_DIR.glob('*slimfriend.root')):
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
for stale in ('hLeadJetSPVA_QCDHtoInv.slimfriend.pdf',
              'hLeadJetSPVA_VBFHtoInv.slimfriend.pdf',
              'hLeadJetSPVA_comparison.pdf',
              'hLeadJetSPVA_comparison.png',
              'hLeadJetSPVM_QCDHtoInv.slimfriend.pdf',
              'hLeadJetSPVM_VBFHtoInv.slimfriend.pdf',
              'jetSPVA_a-1p0b2p0_QCDHtoInv.slimfriend.pdf',
              'jetSPVA_a-1p0b2p0_VBFHtoInv.slimfriend.pdf',
              'jetSPVA_a-2p0b2p0_QCDHtoInv.slimfriend.pdf',
              'jetSPVA_a-2p0b2p0_VBFHtoInv.slimfriend.pdf',
              'jetSPVA_a1p0b-1p0_QCDHtoInv.slimfriend.pdf',
              'jetSPVA_a1p0b-1p0_VBFHtoInv.slimfriend.pdf',
              'jetSPVA_a1p0b1p0_QCDHtoInv.slimfriend.pdf',
              'jetSPVA_a1p0b1p0_VBFHtoInv.slimfriend.pdf',
              'jetSPVA_a1p0b2p0_QCDHtoInv.slimfriend.pdf',
              'jetSPVA_a1p0b2p0_VBFHtoInv.slimfriend.pdf',
              'jetSPVA_a2p0b2p0_QCDHtoInv.slimfriend.pdf',
              'jetSPVA_a2p0b2p0_VBFHtoInv.slimfriend.pdf'):
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
    fpath = DATA_DIR / f'{s["name"]}.slimfriend.root'
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
for ext in ('pdf', 'png'):
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
        fpath = DATA_DIR / f'{sample_name}.slimfriend.root'
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

    ax.hist(bin_centers, bins=bins, weights=np.where(np.isnan(ratio), 0, ratio),
            histtype='step', linewidth=3, color=proc_colors[proc],
            label=proc_labels[proc])
    ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt='none',
                color=proc_colors[proc], capsize=3, capthick=1, elinewidth=1)

ax.axhline(1.0, color='gray', linewidth=1, linestyle='--')

apply_style(ax,
            xlabel=r'$|\theta_s|$ [rad]',
            ylabel='VBF / QCD',
            title='',
            xlim=(0, np.pi),
            legend_loc='upper right')
plt.tight_layout()
for ext in ('pdf', 'png'):
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
    fpath = DATA_DIR / f'{s["name"]}.slimfriend.root'
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
for ext in ('pdf', 'png'):
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
        fpath = DATA_DIR / f'{sample_name}.slimfriend.root'
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

    ax.hist(spvm_bin_centers, bins=spvm_bins, weights=np.where(np.isnan(ratio), 0, ratio),
            histtype='step', linewidth=3, color=proc_colors[proc],
            label=proc_labels[proc])
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
for ext in ('pdf', 'png'):
    outpath = str(FIGS / f'jetSPVM_ratio.{ext}')
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
plt.close(fig)

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
