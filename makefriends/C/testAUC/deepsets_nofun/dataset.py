"""
DijetDataset for DeepSetsNoFun.
Loads 4 constituent features (jcsW excluded); jet-level scalars are loaded but
never passed to the model.
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

JET_FEATURES     = ['jetEta', 'jetM', 'jetNC', 'jetPVM', 'jetPt', 'jetSPVA']
CONSTIT_FEATURES = ['jcsDEta', 'jcsDPhi', 'jcsM', 'jcsPt']   # jcsW removed

N_JET_FEAT     = len(JET_FEATURES)
N_CONSTIT_FEAT = len(CONSTIT_FEATURES)
NC_MAX         = 80

BKG_FILE = '../../friends/QCDHtoInv.h5'
SIG_FILE = '../../friends/VBFHtoInv.h5'


def load_features(h5path):
    with h5py.File(h5path, 'r') as f:
        jets = {k: f[k][:, 0].astype(np.float32)       for k in JET_FEATURES}
        jcs  = {k: f[k][:, 0, :].astype(np.float32)    for k in CONSTIT_FEATURES}
        kw   = f['kWeight'][:].astype(np.float32)

    jets['jetEta'] = np.abs(jets['jetEta'])
    x_jet = np.stack([jets[k] for k in JET_FEATURES],     axis=1)
    x_jcs = np.stack([jcs[k]  for k in CONSTIT_FEATURES], axis=2)
    return x_jet, x_jcs, kw


class DijetDataset(Dataset):
    def __init__(self, x_jcs, mask_jcs, y):
        self.x_jcs = torch.from_numpy(x_jcs)
        self.mask  = torch.from_numpy(mask_jcs)
        self.y     = torch.from_numpy(y).unsqueeze(1)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.x_jcs[idx], self.mask[idx], self.y[idx]
