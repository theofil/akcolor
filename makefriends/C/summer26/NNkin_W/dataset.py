"""
Dataset for NNkin: loads QCDWjj_herwig.h5 + VBFW_herwig.h5.
Uses dijet event-level kinematics and leading/subleading jet 4-vectors.

Event features (4):
  dPhijj, dYjj, mjj, ptjj

Jet features per jet (4 × 2 jets = 8):
  |jetEta|, jetM, jetPhi, jetPt  for jet 0, then jet 1

Total input: N_FEAT = 12, order matches KinNN docstring.

Labels: 0 = QCD background, 1 = VBF signal
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

EVENT_FEATURES = ["dPhijj", "dYjj", "mjj", "ptjj"]
JET_FEATURES   = ["jetEta", "jetM", "jetPhi", "jetPt"]

N_FEAT = len(EVENT_FEATURES) + 2 * len(JET_FEATURES)  # 12

BKG_FILE = "../../friends/summer26/QCDWjj_herwig.h5"
SIG_FILE = "../../friends/summer26/VBFW_herwig.h5"


def load_features(h5path):
    """Return x (N, 12) and kWeight (N,)."""
    with h5py.File(h5path, "r") as f:
        ev   = np.stack([f[k][:].astype(np.float32) for k in EVENT_FEATURES], axis=1)
        jets = np.stack([f[k][:].astype(np.float32) for k in JET_FEATURES],   axis=2)
        kw   = f["kWeight"][:].astype(np.float32)

    jets[:, :, 0] = np.abs(jets[:, :, 0])   # |jetEta|
    x = np.concatenate([ev, jets[:, 0, :], jets[:, 1, :]], axis=1)  # (N, 12)
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
