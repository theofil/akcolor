import uproot, h5py, numpy as np, glob, os, sys

SCALAR_JETS  = ['jetPt', 'jetEta', 'jetPhi', 'jetM', 'jetSPVA', 'jetPVA', 'jetPVM', 'jetNC']
CONSTIT_VARS = ['jcsPt', 'jcsDEta', 'jcsDEtaRaw', 'jcsDPhi', 'jcsM', 'jcsW']
SCALAR_EVT   = ['kWeight', 'mjj', 'dYjj', 'dPhijj', 'ptjj',
                'bosonM', 'bosonY', 'bosonEta', 'bosonPt', 'bosonPhi']

folder = sys.argv[1] if len(sys.argv) > 1 else 'friends/summer26'
pattern = os.path.join(folder, '*.friend.root')

for root_path in sorted(glob.glob(pattern)):
    out_path = root_path.replace('.friend.root', '.h5')
    print(f'Reading {root_path}...')
    with uproot.open(root_path) as f:
        arrays = f['events'].arrays(library='np')

    with h5py.File(out_path, 'w') as h:
        for var in SCALAR_JETS:
            data = np.stack([arrays[f'{var}0'], arrays[f'{var}1'], arrays[f'{var}2']], axis=1)
            h.create_dataset(var, data=data, compression='gzip', compression_opts=4)
        for var in CONSTIT_VARS:
            data = np.stack([arrays[f'{var}0'], arrays[f'{var}1'], arrays[f'{var}2']], axis=1)
            h.create_dataset(var, data=data, compression='gzip', compression_opts=4)
        for var in SCALAR_EVT:
            h.create_dataset(var, data=arrays[var].astype(np.float32),
                             compression='gzip', compression_opts=4)
        h.create_dataset('leadJetIndex', data=arrays['leadJetIndex'].astype(np.int32),
                         compression='gzip', compression_opts=4)

    print(f'  -> {out_path}  ({os.path.getsize(out_path)/1e6:.1f} MB)')
