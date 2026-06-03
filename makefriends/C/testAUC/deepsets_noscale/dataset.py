"""
DijetDataset for DeepSetsNoScale.
Constituent features: jcsDEta, jcsDPhi, jcsPt/jetPt  (3 features).
jcsM and jcsW are excluded. jcsPt is divided by the leading-jet pT before stacking.
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

N_CONSTIT_FEAT = 3
NC_MAX         = 80

BKG_FILE = '../../friends/QCDHtoInv.h5'
SIG_FILE = '../../friends/VBFHtoInv.h5'


def load_features(h5path):
    with h5py.File(h5path, 'r') as f:
        jetPt   = f['jetPt']  [:, 0].astype(np.float32)       # (N,)  leading-jet pT
        jcsDEta = f['jcsDEta'][:, 0, :].astype(np.float32)    # (N, 80)
        jcsDPhi = f['jcsDPhi'][:, 0, :].astype(np.float32)    # (N, 80)
        jcsPt   = f['jcsPt']  [:, 0, :].astype(np.float32)    # (N, 80)  raw, for fraction
        kw      = f['kWeight'][:].astype(np.float32)

    jcsPt_frac = jcsPt / jetPt[:, np.newaxis]                  # (N, 80)  dimensionless
    x_jcs = np.stack([jcsDEta, jcsDPhi, jcsPt_frac], axis=2)  # (N, 80, 3)
    return x_jcs, kw


class DijetDataset(Dataset):
    def __init__(self, x_jcs, mask_jcs, y):
        self.x_jcs = torch.from_numpy(x_jcs)
        self.mask  = torch.from_numpy(mask_jcs)
        self.y     = torch.from_numpy(y).unsqueeze(1)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.x_jcs[idx], self.mask[idx], self.y[idx]
