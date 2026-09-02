#!/usr/bin/env python3
"""
Verify that the current friends/summer26/ ROOT friend files are an exact
event-level skim (subset) of the old, looser-selection files preserved at
friends/summer26/prePLB/.

Context: on 2026-09-01 the dijet preselection in makefriends.cpp was
tightened from `mjj>0` (no real mass cut, just >=2 opposite-eta jets) to
`mjj>400 & |dYjj|>2.5`. The only code change was the cut itself (see
`git diff HEAD -- makefriends.cpp`); no jet definition, kWeight formula, or
constituent-filling logic changed. So the new selection is a strict logical
subset of the old one, and both old and new files already store mjj/dYjj
per event (the old cut just didn't use them) -- meaning this check can be
done entirely offline, by re-filtering the OLD arrays and comparing
row-for-row against the NEW arrays, with no reprocessing required.

Event order is deterministic and identical between the two productions
(single-threaded event loop, no RNG, same run.summer26.sh job partitioning
and hadd merge order for both), so after masking, OLD[mask] should align
positionally with NEW for every branch.

NN score columns (NNj*, NNjj*, NNjB*, NNjjBj*, ...) are added post-hoc only
to the .h5 files by each net's save_scores.py -- they never exist in the
.friend.root files. So this ROOT-level comparison automatically excludes
them; there is nothing to explicitly skip.

Writes an ASCII table to figs/summer26/verify_prePLB_skim.txt.
"""

import pathlib
import numpy as np
import uproot

HERE     = pathlib.Path(__file__).parent
NEW_DIR  = HERE.parent / 'friends' / 'summer26'
OLD_DIR  = pathlib.Path('/afs/cern.ch/user/t/theofil/files/akcolor/summer26/friends/prePLB')
OUT_FILE = HERE.parent / 'figs' / 'summer26' / 'verify_prePLB_skim.txt'

MJJ_MIN  = 400.0
DYJJ_MIN = 2.5

SAMPLES = [
    'VBFH_herwig', 'VBFH_mg5_pythia',
    'VBFW_herwig', 'VBFW_mg5_pythia',
    'VBFZ_herwig', 'VBFZ_mg5_pythia',
    'QCDHjj_herwig', 'QCDHjj_mg5_pythia',
    'QCDWjj_herwig', 'QCDWjj_mg5_pythia',
    'QCDZjj_herwig', 'QCDZjj_mg5_pythia',
]


def compare_sample(name):
    old_path = OLD_DIR / f'{name}.friend.root'
    new_path = NEW_DIR / f'{name}.friend.root'
    if not old_path.exists() or not new_path.exists():
        print(f'  MISSING: {old_path if not old_path.exists() else new_path}')
        return dict(name=name, ok=False, note='missing file')

    with uproot.open(str(old_path)) as fo, uproot.open(str(new_path)) as fn:
        told, tnew = fo['events'], fn['events']
        n_old, n_new = told.num_entries, tnew.num_entries

        mjj_old  = told['mjj'].array(library='np')
        dYjj_old = told['dYjj'].array(library='np')
        mask = (mjj_old > MJJ_MIN) & (np.abs(dYjj_old) > DYJJ_MIN)
        n_pass = int(mask.sum())

        counts_match = (n_pass == n_new)

        old_branches = set(told.keys())
        new_branches = set(tnew.keys())
        if old_branches != new_branches:
            return dict(name=name, ok=False, n_old=n_old, n_pass=n_pass, n_new=n_new,
                         counts_match=counts_match,
                         note=f'branch schema differs: only-old={old_branches - new_branches}, '
                              f'only-new={new_branches - old_branches}')

        first_mismatch = None
        max_kweight_diff = 0.0
        if counts_match:
            for branch in sorted(old_branches):
                old_arr = told[branch].array(library='np')[mask]
                new_arr = tnew[branch].array(library='np')

                if old_arr.dtype.kind in 'fc':
                    same = np.array_equal(old_arr, new_arr)
                    if not same:
                        same = bool(np.allclose(old_arr, new_arr, rtol=1e-6, atol=1e-6,
                                                 equal_nan=True))
                else:
                    same = np.array_equal(old_arr, new_arr)

                if not same:
                    diff = np.abs(old_arr.astype(np.float64) - new_arr.astype(np.float64))
                    idx = int(np.argmax(diff.reshape(diff.shape[0], -1).max(axis=1))) \
                        if diff.ndim > 1 else int(np.argmax(diff))
                    first_mismatch = (branch, idx, float(np.max(diff)))
                    break

                del old_arr, new_arr

        kweight_old_sum = float(told['kWeight'].array(library='np')[mask].sum())
        kweight_new_sum = float(tnew['kWeight'].array(library='np').sum())
        max_kweight_diff = abs(kweight_old_sum - kweight_new_sum)

    all_match = counts_match and (first_mismatch is None)
    return dict(name=name, ok=all_match, n_old=n_old, n_pass=n_pass, n_new=n_new,
                counts_match=counts_match, first_mismatch=first_mismatch,
                kweight_old_sum=kweight_old_sum, kweight_new_sum=kweight_new_sum,
                kweight_diff=max_kweight_diff, note='')


rows = []
for name in SAMPLES:
    print(f'Checking {name} ...')
    r = compare_sample(name)
    rows.append(r)
    if r.get('note') == 'missing file':
        continue
    status = 'PASS' if r['ok'] else 'FAIL'
    print(f'  {status}  n_old={r["n_old"]}  n_old_pass_cut={r["n_pass"]}  n_new={r["n_new"]}  '
          f'counts_match={r["counts_match"]}  '
          f'kWeight_old_sum={r["kweight_old_sum"]:.6g}  kWeight_new_sum={r["kweight_new_sum"]:.6g}  '
          f'kWeight_diff={r["kweight_diff"]:.3g}'
          + (f'  first_mismatch={r["first_mismatch"]}' if r.get('first_mismatch') else ''))

# ── Write ASCII summary ────────────────────────────────────────────────────
header = (
    '# Verification that friends/summer26/<sample>.friend.root is an exact skim of\n'
    '# friends/summer26/prePLB/<sample>.friend.root under the offline cut\n'
    f'# mjj > {MJJ_MIN:g} & |dYjj| > {DYJJ_MIN:g}\n'
    '#\n'
    '# n_old          = entries in the old (prePLB) friend file\n'
    '# n_old_pass_cut = old entries surviving the offline mjj/dYjj mask\n'
    '# n_new          = entries in the current friend file\n'
    '# result          = PASS if counts match AND every branch is identical\n'
    '#                    (row-for-row, after masking) between old[mask] and new\n'
    '#\n'
    f'{"sample":<22s}  {"n_old":>10s}  {"n_old_pass_cut":>14s}  {"n_new":>10s}  '
    f'{"result":>6s}  {"first_mismatch":<40s}\n'
)
lines = []
for r in rows:
    if r.get('note') == 'missing file':
        lines.append(f'{r["name"]:<22s}  {"MISSING":>10s}')
        continue
    result = 'PASS' if r['ok'] else 'FAIL'
    mismatch_str = str(r['first_mismatch']) if r.get('first_mismatch') else '-'
    lines.append(
        f'{r["name"]:<22s}  {r["n_old"]:>10d}  {r["n_pass"]:>14d}  {r["n_new"]:>10d}  '
        f'{result:>6s}  {mismatch_str:<40s}'
    )

with open(OUT_FILE, 'w') as fout:
    fout.write(header + '\n'.join(lines) + '\n')
print(f'\nWrote {OUT_FILE}')

n_fail = sum(1 for r in rows if not r.get('ok', False))
print(f'\n{len(rows) - n_fail}/{len(rows)} samples PASS' + (f', {n_fail} FAIL' if n_fail else ''))
