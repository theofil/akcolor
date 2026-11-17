"""Convert nov26 friend trees (makefriends_nov26) to HDF5.

    python3 root_to_hdf5_nov26.py friends/nov26/VBFZ_herwig.friend.root [...]
    python3 root_to_hdf5_nov26.py friends/nov26      # every *.friend.root in the folder

Writes <name>.h5 next to each <name>.friend.root. An existing h5 is deleted first
(overwriting in place on EOS leaves .sys.v# version dirs behind). The tree is read in
chunks, so an inclusive 1M-event file never has to fit in memory at once.

Layout: 2 jet slots, NC_MAX = 80 constituents per jet, gzip level 4. Dataset names
match the summer26 h5 files (summer26 has 3 slots).
    (N, 2)      jetPt jetEta jetPhi jetM jetSPVA jetPVA jetPVM
                jetBCin jetBCout jetBCmet jetBCdphi (float32), jetNC (int32)
    (N, 2, 80)  jcsPt jcsDEta jcsDEtaRaw jcsDPhi jcsM jcsW (float32), jcsPFType jcsCharge (int8)
    (N,)        kWeight mjj dYjj dPhijj ptjj bosonM bosonY bosonEta bosonPt bosonPhi (float32),
                nJets nJetsAll leadJetIndex (int32), evIdx (int64), gate (bool)
"""
import glob
import os
import sys

import h5py
import numpy as np
import uproot

NSLOT = 2
NC_MAX = 80
CHUNK = 1024          # HDF5 chunk length along N
STEP = 128 * CHUNK    # entries read per iteration (multiple of CHUNK -> whole-chunk writes)

JET_VARS = {'jetPt': np.float32, 'jetEta': np.float32, 'jetPhi': np.float32, 'jetM': np.float32,
            'jetSPVA': np.float32, 'jetPVA': np.float32, 'jetPVM': np.float32,
            'jetBCin': np.float32, 'jetBCout': np.float32, 'jetBCmet': np.float32,
            'jetBCdphi': np.float32, 'jetNC': np.int32}
JCS_VARS = {'jcsPt': np.float32, 'jcsDEta': np.float32, 'jcsDEtaRaw': np.float32,
            'jcsDPhi': np.float32, 'jcsM': np.float32, 'jcsW': np.float32,
            'jcsPFType': np.int8, 'jcsCharge': np.int8}
EVT_VARS = {'kWeight': np.float32, 'mjj': np.float32, 'dYjj': np.float32, 'dPhijj': np.float32,
            'ptjj': np.float32, 'bosonM': np.float32, 'bosonY': np.float32,
            'bosonEta': np.float32, 'bosonPt': np.float32, 'bosonPhi': np.float32,
            'nJets': np.int32, 'nJetsAll': np.int32, 'leadJetIndex': np.int32,
            'evIdx': np.int64, 'gate': np.bool_}

SLOT_GROUPS = ((JET_VARS, (NSLOT,)), (JCS_VARS, (NSLOT, NC_MAX)))
BRANCHES = [f'{v}{s}' for group, _ in SLOT_GROUPS for v in group for s in range(NSLOT)] \
           + list(EVT_VARS)


def convert(root_path):
    out_path = root_path[:-len('.friend.root')] + '.h5'
    print(f'Reading {root_path}...')
    if os.path.exists(out_path):
        os.remove(out_path)
    with uproot.open(root_path) as f, h5py.File(out_path, 'w') as h:
        tree = f['events']
        n = tree.num_entries
        dsets = {}
        for group, tail in SLOT_GROUPS + ((EVT_VARS, ()),):
            for v, dt in group.items():
                dsets[v] = h.create_dataset(v, shape=(n,) + tail, maxshape=(None,) + tail,
                                            dtype=dt, chunks=(CHUNK,) + tail,
                                            compression='gzip', compression_opts=4)
        start = 0
        for arr in tree.iterate(BRANCHES, step_size=STEP, library='np'):
            m = len(arr['kWeight'])
            rows = slice(start, start + m)
            for group, _ in SLOT_GROUPS:
                for v, dt in group.items():
                    dsets[v][rows] = np.stack([arr[f'{v}{s}'] for s in range(NSLOT)],
                                              axis=1).astype(dt)
            for v, dt in EVT_VARS.items():
                dsets[v][rows] = arr[v].astype(dt)
            start += m
        if start != n:
            raise RuntimeError(f'{root_path}: wrote {start} of {n} entries')
    print(f'  -> {out_path}  ({n} events, {os.path.getsize(out_path) / 1e6:.1f} MB)')


def main(args):
    if not args:
        sys.exit(__doc__)
    paths = []
    for a in args:
        if os.path.isdir(a):
            paths += sorted(glob.glob(os.path.join(a, '*.friend.root')))
        elif a.endswith('.friend.root') and os.path.isfile(a):
            paths.append(a)
        else:
            sys.exit(f'not a *.friend.root file or a folder: {a}')
    if not paths:
        sys.exit('no *.friend.root files found')
    for p in paths:
        convert(p)


if __name__ == '__main__':
    main(sys.argv[1:])
