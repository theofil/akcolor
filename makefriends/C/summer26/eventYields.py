#!/usr/bin/env python3
import pathlib
import numpy as np
import uproot

DATA_DIR = pathlib.Path(__file__).parent.parent / 'friends' / 'summer26'
OUT      = pathlib.Path(__file__).parent / 'eventYields.txt'

samples = [
    {'name': 'QCDHjj_mg5_pythia', 'group': 'QCD'},
    {'name': 'VBFH_mg5_pythia',   'group': 'VBF'},
]

jj_vars = [
    ('mjj',    0.0,       800.0),
    ('dYjj',  -8.0,         8.0),
    ('dPhijj', -np.pi,    np.pi),
    ('ptjj',   0.0,       300.0),
]

rows = []
for var, lo, hi in jj_vars:
    for s in samples:
        fpath = DATA_DIR / f'{s["name"]}.friend.root'
        if not fpath.exists():
            print(f'  Skipping {fpath} (not found)')
            continue
        with uproot.open(str(fpath)) as rf:
            if 'events' not in rf:
                print(f'  events tree not found in {fpath}')
                continue
            val = rf['events'][var].array(library='np').astype(float)
            kw  = rf['events']['kWeight'].array(library='np').astype(float)
        mask     = np.ones(len(val), dtype=bool)  # all events (first/last bins are overflow/underflow)
        integral = float(np.sum(kw[mask]))
        raw      = int(np.sum(mask))
        rows.append((var, s['name'], s['group'], integral, raw))
        print(f'  {var:10s}  {s["name"]:25s}  integral={integral:.6e}  raw={raw}')

header = (
    '# Event yields -- summer26\n'
    '# integral = sum(kWeight) for all events (first/last bins absorb overflow/underflow)\n'
    '#\n'
    f'{"variable":<12s}  {"sample":<25s}  {"group":<6s}  {"integral":>14s}  {"raw_entries":>12s}\n'
)
lines = []
for var, name, group, integral, raw in rows:
    lines.append(f'{var:<12s}  {name:<25s}  {group:<6s}  {integral:>14.6e}  {raw:>12d}')

OUT.write_text(header + '\n'.join(lines) + '\n')
print('Wrote', OUT)
