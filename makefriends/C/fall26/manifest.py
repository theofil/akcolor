#!/usr/bin/env python3
"""Fall2026 samples: ordered chunk list with normalization (lxplus; no ROOT needed).

    python3 fall26/manifest.py jobs VBFZ QCDZjj    # one job line per chunk, for run.fall26.sh
    python3 fall26/manifest.py summary             # datasets: chunks, events, sigma, N_eff

The Fall2026 package (README.txt in FALL2026) stores each dataset (process x chain) as ROOT
chunks. An event of a chunk gets kWeight = evweight * pb_per_raw_weight, with the factor of the
chunk's component; W+ and W- are separate components with different factors. Chunks are ordered
by (component, index), and a chunk's evOffset is the number of dataset events before it, so
evIdx = evOffset + entry is the row of the merged dataset.

Job line fields: sample path pb_per_raw_weight chunk evOffset boson charge events
"""
import csv
import json
import os
import sys

FALL2026 = '/eos/home-a/apapaefs/Projects/PullPheno/Fall2026'
GENERATOR = {'herwig': 'herwig', 'pythia-on': 'mg5_pythia'}      # chain -> sample suffix (as nov26)
BOSON = {'VBFH': 'H', 'QCDHjj': 'H', 'VBFW': 'W', 'QCDWjj': 'W', 'VBFZ': 'Z', 'QCDZjj': 'Z'}
CHAINS = tuple(GENERATOR)


def sample_name(process, chain):
    return f'{process}_{GENERATOR[chain]}'


def load_manifest(directory=FALL2026):
    """manifest.json, after checking that the package is complete and verified."""
    with open(os.path.join(directory, 'manifest.json')) as f:
        m = json.load(f)
    with open(os.path.join(directory, 'verification.json')) as f:
        v = json.load(f)
    if m.get('status') != 'complete' or v.get('passed') is not True:
        raise SystemExit(f'{directory}: manifest status {m.get("status")!r}, '
                         f'verification passed = {v.get("passed")!r}')
    return m


def charge(component):
    """+1 / -1 for the W+ / W- components, 0 otherwise."""
    return +1 if component.endswith('plus') else -1 if component.endswith('minus') else 0


def chunks(process, chain, manifest=None, directory=FALL2026):
    """Chunk records of one dataset ordered by (component, index), with the absolute 'path',
    'evOffset', 'charge', 'boson' and 'sample' added."""
    m = manifest or load_manifest(directory)
    rows = sorted((r for r in m['files'] if (r['process'], r['chain']) == (process, chain)),
                  key=lambda r: (r['component'], r['index']))
    if not rows:
        raise SystemExit(f'no chunks for {process}/{chain}')
    out, offset = [], 0
    for r in rows:
        if r['events'] % 2:          # makefriends_fall26 needs even offsets (slot randomization)
            raise SystemExit(f'{r["relative_path"]}: odd number of events ({r["events"]})')
        out.append(dict(r, path=os.path.join(directory, r['relative_path']), evOffset=offset,
                        charge=charge(r['component']), boson=BOSON[process],
                        sample=sample_name(process, chain)))
        offset += r['events']
    return out


def cross_sections(directory=FALL2026):
    """{(process, chain): (sigma_pb, integration_error_pb, events)} from cross_sections.csv."""
    with open(os.path.join(directory, 'cross_sections.csv')) as f:
        return {(r['process'], r['chain']): (float(r['cross_section_pb']),
                                             float(r['integration_error_pb']), int(r['events']))
                for r in csv.DictReader(f)}


def main(argv):
    if not argv or argv[0] not in ('jobs', 'summary'):
        sys.exit(__doc__)
    m = load_manifest()
    if argv[0] == 'jobs':
        for process in argv[1:] or BOSON:
            if process not in BOSON:
                sys.exit(f'unknown process {process}; choose from {", ".join(BOSON)}')
            for chain in CHAINS:
                for r in chunks(process, chain, m):
                    print(r['sample'], r['path'], repr(r['pb_per_raw_weight']), r['index'],
                          r['evOffset'], r['boson'], r['charge'], r['events'])
        return
    xs = cross_sections()
    print(f'{"sample":20s} {"chunks":>6s} {"events":>9s} {"sigma (pb)":>12s} {"N_eff/N":>8s}')
    for process in BOSON:
        for chain in CHAINS:
            rows = chunks(process, chain, m)
            sw = sum(r['sum_raw_weights'] for r in rows)
            sw2 = sum(r['sum_raw_weights_squared'] for r in rows)
            n = sum(r['events'] for r in rows)
            print(f'{sample_name(process, chain):20s} {len(rows):6d} {n:9d} '
                  f'{xs[process, chain][0]:12.6g} {sw * sw / sw2 / n:8.3f}')


if __name__ == '__main__':
    main(sys.argv[1:])
