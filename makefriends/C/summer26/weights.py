#!/usr/bin/env python3
"""Produce weight-variable distributions and weights.html."""
import pathlib

from checksum import MOTHER_XS

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import uproot

DATA_DIR = pathlib.Path(__file__).parent.parent / 'friends' / 'summer26'
FIGS     = pathlib.Path(__file__).parent.parent / 'figs'   / 'summer26'
FIGS.mkdir(parents=True, exist_ok=True)

CAMPAIGN23 = pathlib.Path('/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged')
CAMPAIGN30 = pathlib.Path('/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00030/merged')
CAMPAIGN31 = pathlib.Path('/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00031/merged')

EW_SAMPLES_BY_PROC = {
    'H': [
        ('QCDHjj_mg5_pythia', CAMPAIGN31 / 'QCDHjj_mg5_pythia.root', 'QCD H  MG5+Pythia'),
        ('QCDHjj_herwig',     CAMPAIGN30 / 'QCDHjj_herwig.root',     'QCD H  Herwig'),
        ('VBFH_mg5_pythia',   CAMPAIGN23 / 'VBFH_mg5_pythia.root',   'VBF H  MG5+Pythia'),
        ('VBFH_herwig',       CAMPAIGN23 / 'VBFH_herwig.root',       'VBF H  Herwig'),
    ],
    'W': [
        ('QCDWjj_mg5_pythia', CAMPAIGN23 / 'QCDWjj_mg5_pythia.root', 'QCD W  MG5+Pythia'),
        ('QCDWjj_herwig',     CAMPAIGN23 / 'QCDWjj_herwig.root',     'QCD W  Herwig'),
        ('VBFW_mg5_pythia',   CAMPAIGN23 / 'VBFW_mg5_pythia.root',   'VBF W  MG5+Pythia'),
        ('VBFW_herwig',       CAMPAIGN23 / 'VBFW_herwig.root',       'VBF W  Herwig'),
    ],
    'Z': [
        ('QCDZjj_mg5_pythia', CAMPAIGN23 / 'QCDZjj_mg5_pythia.root', 'QCD Z  MG5+Pythia'),
        ('QCDZjj_herwig',     CAMPAIGN23 / 'QCDZjj_herwig.root',     'QCD Z  Herwig'),
        ('VBFZ_mg5_pythia',   CAMPAIGN23 / 'VBFZ_mg5_pythia.root',   'VBF Z  MG5+Pythia'),
        ('VBFZ_herwig',       CAMPAIGN23 / 'VBFZ_herwig.root',       'VBF Z  Herwig'),
    ],
}


def plot_weight_panel(weight_name, samples_by_proc, read_fn):
    for proc, samples in samples_by_proc.items():
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes = axes.flatten()
        for ax, (name, epath, label) in zip(axes, samples):
            arr = read_fn(name, epath)
            if arr is None:
                ax.set_visible(False)
                continue
            mu, std, n, wsum = float(arr.mean()), float(arr.std()), len(arr), float(arr.sum())
            xs = MOTHER_XS.get(name)
            ax.hist(arr, bins=100, histtype='step', linewidth=2, color='steelblue', density=True)
            ax.axvline(mu, color='tomato', linewidth=1.5, linestyle='--')
            ax.set_xlabel(weight_name)
            ax.set_ylabel('density')
            ax.set_yscale('log')
            ax.set_title(label, fontsize=10)
            xs_line = f'\nσ    = {xs:.5g} pb' if xs is not None else ''
            ax.text(0.97, 0.95,
                    f'N = {n:,}\nmean = {mu:.5g}\nstd  = {std:.5g}\nΣ    = {wsum:.5g}{xs_line}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=9,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
        fig.suptitle(f'{weight_name} distributions — process {proc}', fontsize=12)
        plt.tight_layout()
        outpath = str(FIGS / f'{weight_name}_{proc}.pdf')
        fig.savefig(outpath, bbox_inches='tight')
        print('Saved', outpath)
        plt.close(fig)


def read_evweight(name, epath):
    if not epath.exists():
        print(f'  Skipping {epath} (not found)')
        return None
    with uproot.open(str(epath)) as rf:
        return rf['Data']['evweight'].array(library='np').astype(float)


def read_kweight(name, epath):
    fpath = DATA_DIR / f'{name}.friend.root'
    if not fpath.exists():
        print(f'  Skipping {fpath} (not found)')
        return None
    with uproot.open(str(fpath)) as rf:
        if 'events' not in rf:
            print(f'  events tree not found in {fpath}')
            return None
        return rf['events']['kWeight'].array(library='np').astype(float)


plot_weight_panel('evweight', EW_SAMPLES_BY_PROC, read_evweight)
plot_weight_panel('kWeight',  EW_SAMPLES_BY_PROC, read_kweight)

# ── Print per-sample summary table ───────────────────────────────────────────
all_names = [name for samples in EW_SAMPLES_BY_PROC.values() for name, _, _ in samples]
seen = set()
unique_names = [n for n in all_names if not (n in seen or seen.add(n))]

print()
print('| Sample | N_evweight | evweight mean | evweight std | N_kWeight | kWeight mean | kWeight std |')
print('|--------|-----------|--------------|-------------|----------|------------|------------|')
for proc, samples in EW_SAMPLES_BY_PROC.items():
    for name, epath, label in samples:
        ew = read_evweight(name, epath)
        kw = read_kweight(name, epath)
        ev_str = f'{len(ew):,} | {ew.mean():.5g} | {ew.std():.5g}' if ew is not None else '— | — | —'
        kw_str = f'{len(kw):,} | {kw.mean():.5g} | {kw.std():.5g}' if kw is not None else '— | — | —'
        print(f'| `{name}` | {ev_str} | {kw_str} |')

# ── Generate weights.html ─────────────────────────────────────────────────────
_CSS = '''
*, *::before, *::after { box-sizing: border-box; }
html, body { height: 100%; margin: 0; overflow: hidden;
             font-family: sans-serif; background: #f0f0f0; }

.layout { display: flex; height: 100%; }

nav.toc {
  width: 210px; flex-shrink: 0;
  background: #1e2a38; color: #c8d8e8;
  overflow-y: auto; padding: 1em 0.8em;
}
nav.toc h2 { font-size: 0.85em; text-transform: uppercase;
             letter-spacing: .08em; color: #7aa0c0; margin: 1em 0 0.4em; }
nav.toc h2:first-child { margin-top: 0; }
nav.toc a {
  display: block; padding: 0.3em 0.4em; margin-bottom: 0.1em;
  color: #c8d8e8; text-decoration: none; font-size: 0.80em;
  border-radius: 3px;
}
nav.toc a:hover { background: #2e4260; color: #fff; }

main.content { flex: 1; overflow-y: auto; padding: 0.6em 1.2em; }

.triplet {
  height: calc(100vh - 2em); min-height: 400px;
  display: flex; flex-direction: column;
  background: white; border: 1px solid #d0d8e0;
  border-radius: 5px; padding: 0 0.6em 0.4em;
  margin-bottom: 0.6em;
}
.triplet-label {
  font-size: 0.75em; color: #888; text-align: center;
  padding: 0.2em 0; flex-shrink: 0;
}

.cols { flex: 1; min-height: 0; display: flex; gap: 3px; }

.cols figure {
  flex: 1; min-height: 0; margin: 0;
  display: flex; flex-direction: column;
}
.cols figcaption {
  text-align: center; font-size: 1.6em;
  flex-shrink: 0; margin-bottom: 0.15em; font-weight: bold;
}
.iframe-wrapper { flex: 1; min-height: 0; position: relative; }
.cols figure iframe { width: 100%; height: 100%; border: none; display: block; }
.iframe-overlay {
  position: absolute; inset: 0; z-index: 10; cursor: pointer;
}
.cols .missing {
  flex: 1; display: flex; align-items: center; justify-content: center;
  color: #bbb; font-size: 0.85em; background: #fafafa;
  border: 1px dashed #ddd; border-radius: 3px;
}
'''

BOSON_COLOR = {'W': '#1a6eb5', 'Z': '#2e8b57', 'H': '#b52020'}
WEIGHT_VARS = ['evweight', 'kWeight']


def _triplet_section(stem):
    lines = [f'<section class="triplet" id="{stem}">',
             f'  <div class="triplet-label">{stem}</div>',
             '  <div class="cols">']
    for b in ['W', 'Z', 'H']:
        fname = f'{stem}_{b}.pdf'
        color = BOSON_COLOR[b]
        lines.append('    <figure>')
        lines.append(f'      <figcaption style="color:{color}">{b}</figcaption>')
        lines.append('      <div class="iframe-wrapper">')
        lines.append(f'        <iframe src="{fname}#toolbar=0&amp;view=Fit" title="{stem} ({b})"></iframe>')
        lines.append(f'        <a class="iframe-overlay" href="{fname}" download title="Download {fname}"></a>')
        lines.append('      </div>')
        lines.append('    </figure>')
    lines += ['  </div>', '</section>']
    return '\n'.join(lines)


_nav_links = '\n'.join(
    f'<a href="#{w}">{w}</a>' for w in WEIGHT_VARS
)
_sections = '\n'.join(_triplet_section(w) for w in WEIGHT_VARS)

_html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>akcolor — weights</title>
<style>{_CSS}</style>
</head>
<body>
<div class="layout">
  <nav class="toc">
    <h2>Weights</h2>
    {_nav_links}
  </nav>
  <main class="content">
    {_sections}
  </main>
</div>
</body>
</html>'''

_html_path = FIGS / 'weights.html'
_html_path.write_text(_html)
print('Saved', _html_path)
