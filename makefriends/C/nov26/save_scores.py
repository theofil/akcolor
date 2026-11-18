#!/usr/bin/env python3
"""Write the out-of-fold scores of the nov26 nets into friends/nov26/<sample>.h5 (lxplus; no model).

    python3 nov26/save_scores.py --net NNjBC
    python3 nov26/save_scores.py --net all        # every net of common.NETS, one h5 write per file

Merges nov26/<net>/scores_fold{0..4}.npz (written by train_fold.py) and writes, for each sample:
  <net>    (N, 2) float32  out-of-fold score of the jet in each slot; -1 for an empty slot
  subfold  (N,)   int8     sub-fold of the event (0-9): model k scored sub-folds 2k and 2k+1
Before writing it checks that every stored jet was scored exactly once, by the model whose
inference sub-folds contain its event. Rerunning root_to_hdf5_nov26.py drops these columns.
"""
import argparse
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402


def merge(sample, n_jets, sub, folds_dir, entry_stop=None):
    """(N, 2) scores from the five fold files; raises if a stored jet is missing or doubled."""
    n = len(n_jets)
    score = np.full((n, 2), -1.0, np.float32)
    count = np.zeros((n, 2), np.int16)
    for k in range(C.N_FOLDS):
        d = np.load(os.path.join(folds_dir, f'scores_fold{k}.npz'))
        ev, slot, s = d[f'{sample}__ev'], d[f'{sample}__slot'].astype(np.int64), d[f'{sample}__score']
        if not np.isin(sub[ev], C.fold_sets(k)[2]).all():
            raise SystemExit(f'{sample}: fold {k} scored events outside its inference sub-folds')
        np.add.at(count, (ev, slot), 1)
        score[ev, slot] = s
    expected = (np.arange(2)[None, :] < n_jets[:, None]).astype(np.int16)
    rows = slice(None) if entry_stop is None else slice(0, entry_stop)
    if not np.array_equal(count[rows], expected[rows]):
        bad = int((count[rows] != expected[rows]).sum())
        raise SystemExit(f'{sample}: {bad} jet slots are not scored exactly once')
    return score


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--net', required=True, choices=list(C.NETS) + ['all'])
    ap.add_argument('--folds-dir', default=None, help='default: nov26/<net> (single net only)')
    ap.add_argument('--dry-run', action='store_true', help='check the merge, write nothing')
    ap.add_argument('--entry-stop', type=int, default=None, help='check only the first N events (smoke tests)')
    args = ap.parse_args()
    nets = list(C.NETS) if args.net == 'all' else [args.net]
    if args.folds_dir and len(nets) > 1:
        raise SystemExit('--folds-dir needs a single --net')
    for sample in C.LABEL:
        path = C.h5_path(sample)
        with h5py.File(path, 'r') as f:
            gate, n_jets, ev_idx = f['gate'][:], f['nJets'][:], f['evIdx'][:]
            old_sub = f['subfold'][:] if 'subfold' in f else None
        if not np.array_equal(ev_idx, np.arange(len(ev_idx))):
            raise SystemExit(f'{sample}: evIdx is not the row index')
        sub = C.assign_subfolds(sample, gate, n_jets)
        if old_sub is not None and not np.array_equal(old_sub, sub):
            raise SystemExit(f'{sample}: existing subfold column differs from common.assign_subfolds')
        out = {'subfold': sub}
        for net in nets:
            score = merge(sample, n_jets, sub, args.folds_dir or os.path.join(C.OUT, net), args.entry_stop)
            filled = score[score >= 0]
            print(f'{sample:18s} {net:13s} {len(filled)} jets scored, mean score {filled.mean():.3f}'
                  + ('  (dry run)' if args.dry_run else ''))
            out[net] = score
        if not args.dry_run:
            C.write_datasets(path, out)
            print(f'  wrote {", ".join(nets)} (N, 2) and subfold (N,) into {path}')


if __name__ == '__main__':
    main()
