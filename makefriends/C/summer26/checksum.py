#!/usr/bin/env python3
"""
Read hSumW and the events tree from each ROOT friend file in friends/summer26/
and compute cross sections and entry counts.

xs_in    = hSumW.Integral(0, N+1)  (includes underflow and overflow bins)

hSumW is filled for every event weighted by kWeight (line 307 in makefriends.cpp),
before the dijet selection, so its integral equals sum(kWeight) over all events.

xs_out   = sum of kWeight from the events TTree (events that passed the full dijet
           selection applied during TTree filling in makefriends.cpp).

eff_w    = xs_out / xs_in  (weighted filter efficiency)
ratio    = n_entries_out / n_entries_in  (unweighted event efficiency, +/- binomial unc.)

Writes an ASCII table to figs/summer26/checksum.txt, overwriting any previous
content so the file always reflects the latest run.
"""

import ctypes
import pathlib
import numpy as np
import uproot
import ROOT

HERE         = pathlib.Path(__file__).parent
DATA_DIR     = HERE.parent / 'friends' / 'summer26'
CHECKSUM_OUT = HERE.parent / 'figs' / 'summer26' / 'checksum.txt'

# Static cross-sections (pb) from run.summer26.sh, keyed by sample name
MOTHER_XS = {
    'VBFH_herwig':       2.97189,
    'VBFH_mg5_pythia':   2.88990,
    'VBFW_herwig':       7.494,
    'VBFW_mg5_pythia':   7.20333,
    'VBFZ_herwig':       1.099,
    'VBFZ_mg5_pythia':   1.08581,
    'QCDHjj_herwig':     5.06051,
    'QCDHjj_mg5_pythia': 4.967,
    'QCDWjj_herwig':     1733.3,
    'QCDWjj_mg5_pythia': 1646.805,
    'QCDZjj_herwig':     340.2,
    'QCDZjj_mg5_pythia': 316.363,
}


rows = []
for name, xs_static in sorted(MOTHER_XS.items()):
    rootpath = DATA_DIR / f'{name}.friend.root'
    if not rootpath.exists():
        print(f'  MISSING: {rootpath}')
        rows.append((name, xs_static, None, None, None, None, None, None, None, None))
        continue
    f = ROOT.TFile.Open(str(rootpath))
    h = f.Get("hSumW")
    t = f.Get("events")
    err           = ctypes.c_double(0.0)
    xs_in         = h.IntegralAndError(0, h.GetNbinsX() + 1, err)
    xs_err        = err.value
    n_entries_in  = int(h.GetEntries())
    n_entries_out = int(t.GetEntries())
    f.Close()

    with uproot.open(str(rootpath)) as uf:
        xs_out = float(uf['events']['kWeight'].array(library='np').sum())

    eff_w       = xs_out / xs_in if xs_in > 0 else 0.0
    ratio       = n_entries_out / n_entries_in if n_entries_in > 0 else 0.0
    sigma_ratio = (float(np.sqrt(ratio * (1.0 - ratio) / n_entries_in))
                   if n_entries_in > 0 else 0.0)

    rows.append((name, xs_static, xs_in, xs_err,
                 n_entries_in, n_entries_out, xs_out, eff_w, ratio, sigma_ratio))
    print(f'  {name:<25s}  xs_static={xs_static:.5g} pb  '
          f'xs_in={xs_in:.5g} pb  '
          f'xs_out={xs_out:.5g} pb  '
          f'eff_w={eff_w:.4f}  ratio={ratio:.4f} +/- {sigma_ratio:.4f}  '
          f'n_in={n_entries_in}  n_out={n_entries_out}')

# ── Write ASCII table to checksum.txt (overwrite on every run) ────────────────
header = (
    '# summer26 cross sections and entry counts\n'
    '# xs_static  = static cross section from run.summer26.sh\n'
    '# xs_in      = hSumW.Integral(0,N+1) - sum of kWeight before dijet selection\n'
    '# xs_out     = sum of kWeight from events TTree (after full dijet selection)\n'
    '# eff_w      = xs_out / xs_in  (weighted filter efficiency)\n'
    '# ratio +/- sigma = n_entries_out / n_entries_in +/- binomial uncertainty\n'
    '#\n'
    f'{"sample":<25s}  {"xs_static [pb]":>14s}  {"xs_in [pb]":>14s}  '
    f'{"xs_out [pb]":>14s}  {"eff_w":>8s}  {"ratio +/- sigma":>16s}  '
    f'{"n_entries_in":>14s}  {"n_entries_out":>14s}\n'
)
lines = []
for (name, xs_static, xs_in, xs_err,
     n_entries_in, n_entries_out, xs_out, eff_w, ratio, sigma_ratio) in rows:
    if xs_in is None:
        lines.append(
            f'{name:<25s}  {xs_static:>14.5g}  {"MISSING":>14s}  {"---":>14s}  '
            f'{"---":>8s}  {"---":>16s}  {"---":>14s}  {"---":>14s}'
        )
    else:
        ratio_str = f'{ratio:.4f} +/- {sigma_ratio:.4f}'
        lines.append(
            f'{name:<25s}  {xs_static:>14.5g}  {xs_in:>14.5g}  {xs_out:>14.5g}  '
            f'{eff_w:>8.4f}  {ratio_str:>16s}  {n_entries_in:>14d}  {n_entries_out:>14d}'
        )

with open(CHECKSUM_OUT, 'w') as fout:
    fout.write(header + '\n'.join(lines) + '\n')
print(f'\nWrote {CHECKSUM_OUT}')
