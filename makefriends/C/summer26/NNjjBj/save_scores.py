#!/usr/bin/env python3
"""
Run NNjjBj inference (both leading jets + 3rd jet + boson scalars) on this
channel's summer26 h5 files (QCD{P}jj / VBF{P}, Herwig + MG5+Pythia) and
write the per-event sigmoid score back into each file.

Per-channel convention: each process is scored by its own Herwig-trained
net (PROCESS below) — no cross-process inference. The Herwig-trained net
is applied to both the Herwig and the MG5+Pythia files of its channel.

Written dataset:
  NNjjBj — event-level score (jets 0+1 with constituents, 3rd-jet
           4-momentum + constituents when present, boson scalars)
"""

import os
import sys
import time
import pickle
import pathlib

import numpy as np
import h5py
import torch

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))

from dataset import load_features
from inference import run_nn
from model import JetNN

PROCESS = "Z"

H5_FILES = [
    HERE / f'../../friends/summer26/QCD{PROCESS}jj_herwig.h5',
    HERE / f'../../friends/summer26/QCD{PROCESS}jj_mg5_pythia.h5',
    HERE / f'../../friends/summer26/VBF{PROCESS}_herwig.h5',
    HERE / f'../../friends/summer26/VBF{PROCESS}_mg5_pythia.h5',
]

KEY = 'NNjjBj'


def write_dataset(path, key, data, attempts=3):
    """Replace the dataset in the h5 file, retrying: EOS FUSE occasionally
    fails HDF5 metadata ops with 'bad object header version number'."""
    for attempt in range(1, attempts + 1):
        try:
            with h5py.File(path, 'a') as f:
                if key in f:
                    del f[key]
                f.create_dataset(key, data=data,
                                 compression='gzip', compression_opts=4)
            return
        except (KeyError, OSError, RuntimeError) as e:
            print(f'  write attempt {attempt}/{attempts} failed: {e}')
            if attempt == attempts:
                raise
            time.sleep(10)


def main():
    with open(HERE / f'scaler_{PROCESS}.pkl', 'rb') as fh:
        scalers = pickle.load(fh)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model  = JetNN().to(device)
    model.load_state_dict(torch.load(HERE / f'best_model_{PROCESS}.pt', map_location=device))
    model.eval()
    print(f'Model loaded, device={device}')

    for path in H5_FILES:
        path = pathlib.Path(path)
        if not path.exists():
            print(f'SKIP (not found): {path}')
            continue
        print(f'Processing {path.name} ...')
        x_j0, x_j1, x_j2, x_bos, x_c0, x_c1, x_c2, _ = load_features(path)
        scores = run_nn(model, scalers, x_j0, x_j1, x_j2, x_bos,
                        x_c0, x_c1, x_c2, device).astype(np.float32)
        print(f'  {len(scores)} events  mean={scores.mean():.4f}  max={scores.max():.4f}')

        write_dataset(path, KEY, scores)
        print(f'  Written {KEY} -> {path.name}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
