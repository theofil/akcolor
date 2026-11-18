#!/usr/bin/env python3
"""Validate the fall26 friends (friends/fall26/<sample>.h5) against the Fall2026 package (lxplus).

    python3 fall26/validate.py                    # checks + per-stage tables (markdown), VBFZ and QCDZjj
    python3 fall26/validate.py --inputs           # also kWeight == evweight * factor, from every input chunk
    python3 fall26/validate.py --processes VBFW QCDWjj

Checks per sample:
  rows        one row per input event (sum of the manifest chunk sizes)
  evIdx       evIdx == arange(N)
  chunk order (chunk, chunkEntry) follow the manifest chunk order, entries 0..n-1 in each chunk
  sigma       sum of kWeight == sigma of cross_sections.csv, to 1e-6 (relative)
  boson       a boson in all but at most 1e-4 of the events, and bosonCharge matches the chunk's W
              component (0 for H, Z). Pythia stores particles with pT > 0.1 GeV only, so a decay
              neutrino emitted along the beam can be missing (7 Z events per million)
  inputs      (--inputs) kWeight == float32(evweight * pb_per_raw_weight), event by event
Tables per stage (all events; gate; VBF preselection: gate, mjj > 400 GeV, |dYjj| > 2.5): sigma,
events, N_eff = (sum w)^2 / sum w^2, and sigma_VBF / sigma_QCD per generator next to nov26.
"""
import argparse
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import manifest as M  # noqa: E402

STAGES = ('all events', 'gate', 'VBF presel.')
GEN_LABEL = {'herwig': 'Herwig', 'pythia-on': 'MG5+Pythia'}


def h5_path(campaign, sample):
    return os.path.join(BASE, 'friends', campaign, f'{sample}.h5')


def stages(path):
    """Event weights and stage masks of one h5 file."""
    with h5py.File(path, 'r') as f:
        kw = f['kWeight'][:].astype(np.float64)
        gate = f['gate'][:]
        vbf = gate & (f['mjj'][:] > 400) & (np.abs(f['dYjj'][:]) > 2.5)
    return kw, {'all events': np.ones(len(kw), bool), 'gate': gate, 'VBF presel.': vbf}


def n_eff(w):
    return w.sum() ** 2 / (w * w).sum() if len(w) else 0.0


def check(process, chain, m, xs, inputs):
    sample = M.sample_name(process, chain)
    rows = M.chunks(process, chain, m)
    with h5py.File(h5_path('fall26', sample), 'r') as f:
        n = len(f['evIdx'])
        ev, ch, ce = f['evIdx'][:], f['chunk'][:], f['chunkEntry'][:]
        kw32 = f['kWeight'][:]
        boson_m, q = f['bosonM'][:], f['bosonCharge'][:]
    sizes = [r['events'] for r in rows]
    sigma = xs[process, chain][0]
    rel = kw32.astype(np.float64).sum() / sigma - 1
    no_boson = int((boson_m <= 0).sum())
    res = {
        'rows': n == sum(sizes),
        'evIdx': np.array_equal(ev, np.arange(n)),
        'chunk order': (np.array_equal(ch, np.repeat([r['index'] for r in rows], sizes))
                        and np.array_equal(ce, np.concatenate([np.arange(s) for s in sizes]))),
        'sigma': abs(rel) < 1e-6,
        'boson': no_boson <= 1e-4 * n and np.array_equal(q, np.repeat([r['charge'] for r in rows], sizes)),
    }
    if inputs:
        import uproot
        w = np.concatenate([uproot.open(r['path'])['Data']['evweight'].array(library='np')
                            * r['pb_per_raw_weight'] for r in rows])
        res['inputs'] = np.array_equal(kw32, w.astype(np.float32))
    print(f'| `{sample}` | {n} | {sigma:.6g} | {rel:+.1e} | {np.median(boson_m[boson_m > 0]):.2f} | {no_boson} | '
          + ' | '.join('ok' if v else '**FAILED**' for v in res.values()) + ' |')
    return all(res.values()), list(res)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--processes', nargs='+', default=['VBFZ', 'QCDZjj'])
    ap.add_argument('--inputs', action='store_true', help='compare kWeight with every input chunk (slow)')
    args = ap.parse_args()
    m, xs = M.load_manifest(), M.cross_sections()
    if len(args.processes) != 2:
        sys.exit('give one VBF and one QCD process, e.g. --processes VBFZ QCDZjj')
    vbf_p, qcd_p = args.processes

    print('### Checks\n')
    names = ['rows', 'evIdx', 'chunk order', 'sigma', 'boson'] + (['inputs'] if args.inputs else [])
    print('| Sample | Rows | σ (pb) | Σ kWeight / σ − 1 | median bosonM | no boson | ' + ' | '.join(names) + ' |')
    print('|---|---|---|---|---|---|' + '---|' * len(names))
    good = True
    for process in args.processes:
        for chain in M.CHAINS:
            ok, _ = check(process, chain, m, xs, args.inputs)
            good &= ok

    print('\n### Stages\n')
    print('| Sample | ' + ' | '.join(f'σ {s} (pb) | N | N_eff' for s in STAGES) + ' |')
    print('|---|' + '---|---|---|' * len(STAGES))
    sig = {}
    for process in args.processes:
        for chain in M.CHAINS:
            sample = M.sample_name(process, chain)
            kw, masks = stages(h5_path('fall26', sample))
            sig['fall26', process, chain] = {s: kw[masks[s]].sum() for s in STAGES}
            print(f'| `{sample}` | ' + ' | '.join(
                f'{kw[masks[s]].sum():.4g} | {int(masks[s].sum())} | {n_eff(kw[masks[s]]):.0f}' for s in STAGES) + ' |')
            old = h5_path('nov26', sample)
            if os.path.exists(old):
                kw_o, masks_o = stages(old)
                sig['nov26', process, chain] = {s: kw_o[masks_o[s]].sum() for s in STAGES}

    print('\n### σ_VBF / σ_QCD\n')
    print('| Generator | Campaign | ' + ' | '.join(STAGES) + ' |')
    print('|---|---|' + '---|' * len(STAGES))
    for chain in M.CHAINS:
        for camp in ('fall26', 'nov26'):
            if (camp, vbf_p, chain) in sig and (camp, qcd_p, chain) in sig:
                v, q = sig[camp, vbf_p, chain], sig[camp, qcd_p, chain]
                print(f'| {GEN_LABEL[chain]} | {camp} | '
                      + ' | '.join(f'{v[s] / q[s]:.3e} (1 : {q[s] / v[s]:.0f})' for s in STAGES) + ' |')
    print('\nall checks passed' if good else '\n**SOME CHECKS FAILED**')
    return 0 if good else 1


if __name__ == '__main__':
    sys.exit(main())
