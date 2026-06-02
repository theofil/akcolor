"""
KinDataset: loads QCDHtoInv.h5 + VBFHtoInv.h5 using only kinematic variables.

Event-level scalars (4):
  dPhijj, dYjj, mjj, ptjj

Per-jet features (4 × 2 jets = 8):
  [j0] jetEta, jetM, jetPhi, jetPt
  [j1] jetEta, jetM, jetPhi, jetPt

Total input dimension: 12

Labels: 0 = QCD background, 1 = VBF signal
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

EVENT_FEATURES = ['dPhijj', 'dYjj', 'mjj', 'ptjj']
JET_FEATURES   = ['jetEta', 'jetM', 'jetPhi', 'jetPt']

N_EVENT_FEAT = len(EVENT_FEATURES)          # 4
N_JET_FEAT   = len(JET_FEATURES)            # 4 per jet
N_FEAT       = N_EVENT_FEAT + 2 * N_JET_FEAT  # 12

BKG_FILE = '../friends/QCDHtoInv.h5'
SIG_FILE = '../friends/VBFHtoInv.h5'


def load_features(h5path):
    """Return x (N, 12) feature matrix and kWeight (N,)."""
    with h5py.File(h5path, 'r') as f:
        ev  = {k: f[k][:].astype(np.float32) for k in EVENT_FEATURES}   # (N,)
        jet = {k: f[k][:].astype(np.float32) for k in JET_FEATURES}     # (N, 2)
        kw  = f['kWeight'][:].astype(np.float32)


    ev_arr = np.stack([ev[k]            for k in EVENT_FEATURES], axis=1)  # (N, 4)
    jet0   = np.stack([jet[k][:, 0]    for k in JET_FEATURES],   axis=1)  # (N, 4)
    jet1   = np.stack([jet[k][:, 1]    for k in JET_FEATURES],   axis=1)  # (N, 4)
    x      = np.concatenate([ev_arr, jet0, jet1], axis=1)                  # (N, 12)
    return x, kw


class KinDataset(Dataset):
    def __init__(self, x, y, w=None):
        self.x = torch.from_numpy(x)
        self.y = torch.from_numpy(y).unsqueeze(1)
        self.w = torch.from_numpy(w) if w is not None else torch.ones(len(y))

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx], self.w[idx]
