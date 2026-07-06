#!/usr/bin/env python3
"""
Run NNPrePro inference on both jets for all summer26 h5 files and write
per-event sigmoid scores back into each file.

Written datasets:
  NNPrePro_jet0 — score for the leading    jet (index 0)
  NNPrePro_jet1 — score for the subleading jet (index 1)
"""

import os
import sys
import pickle
import pathlib

import numpy as np
import h5py
import torch

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))

from dataset import JET_FEATURES, CONSTIT_FEATURES, N_CONSTIT_FEAT
from model import JetNN

PROCESS = "Z"

H5_FILES = [
    HERE / '../../friends/summer26/QCDZjj_herwig.h5',
    HERE / '../../friends/summer26/QCDZjj_mg5_pythia.h5',
    HERE / '../../friends/summer26/QCDWjj_herwig.h5',
    HERE / '../../friends/summer26/QCDWjj_mg5_pythia.h5',
    HERE / '../../friends/summer26/QCDHjj_herwig.h5',
    HERE / '../../friends/summer26/QCDHjj_mg5_pythia.h5',
    HERE / '../../friends/summer26/VBFZ_herwig.h5',
    HERE / '../../friends/summer26/VBFZ_mg5_pythia.h5',
    HERE / '../../friends/summer26/VBFW_herwig.h5',
    HERE / '../../friends/summer26/VBFW_mg5_pythia.h5',
    HERE / '../../friends/summer26/VBFH_herwig.h5',
    HERE / '../../friends/summer26/VBFH_mg5_pythia.h5',
]

BATCH = 512


def load_jet_features(h5path, jet_idx):
    """Return x_jet (N, 3) and x_jcs (N, 80, 3) for the given jet index (0 or 1)."""
    with h5py.File(h5path, 'r') as f:
        jets = {k: f[k][:, jet_idx].astype(np.float32)       for k in JET_FEATURES}
        jcs  = {k: f[k][:, jet_idx, :].astype(np.float32)    for k in CONSTIT_FEATURES}

    jets['jetEta'] = np.abs(jets['jetEta'])
    x_jet = np.stack([jets[k] for k in JET_FEATURES],     axis=1)  # (N, 3)
    x_jcs = np.stack([jcs[k]  for k in CONSTIT_FEATURES], axis=2)  # (N, 80, 3)
    x_jcs[:, :, 2] /= x_jet[:, 2:3]   # jcsPt → jcsPt / jetPt
    return x_jet, x_jcs


def scale_jcs(x_jcs, scaler):
    shape = x_jcs.shape
    return scaler.transform(x_jcs.reshape(-1, N_CONSTIT_FEAT)).reshape(shape).astype(np.float32)


def run_nn(model, scalers, x_jet_raw, x_jcs_raw, device):
    x_jet = scalers['jet'].transform(x_jet_raw).astype(np.float32)
    x_jcs = scale_jcs(x_jcs_raw, scalers['jcs'])
    mask  = (x_jcs_raw[..., 2] > 0)   # jcsPt > 0 before scaling

    x_jet_t = torch.from_numpy(x_jet)
    x_jcs_t = torch.from_numpy(x_jcs)
    mask_t  = torch.from_numpy(mask)

    logits_list = []
    with torch.no_grad():
        for start in range(0, len(x_jet), BATCH):
            end = start + BATCH
            out = model(
                x_jet_t[start:end].to(device),
                x_jcs_t[start:end].to(device),
                mask_t [start:end].to(device),
            )
            logits_list.append(out.cpu())
    logits = torch.cat(logits_list).numpy().ravel()
    return (1.0 / (1.0 + np.exp(-logits))).astype(np.float32)


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
        scores = {}
        for jet_idx in [0, 1]:
            x_jet_raw, x_jcs_raw = load_jet_features(path, jet_idx)
            scores[jet_idx] = run_nn(model, scalers, x_jet_raw, x_jcs_raw, device)
            print(f'  jet{jet_idx}: {len(scores[jet_idx])} events  '
                  f'mean={scores[jet_idx].mean():.4f}  max={scores[jet_idx].max():.4f}')

        with h5py.File(path, 'a') as f:
            for jet_idx in [0, 1]:
                key = f'NNPrePro_jet{jet_idx}'
                if key in f:
                    del f[key]
                f.create_dataset(key, data=scores[jet_idx],
                                 compression='gzip', compression_opts=4)
        print(f'  Written NNPrePro_jet0, NNPrePro_jet1 -> {path.name}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
