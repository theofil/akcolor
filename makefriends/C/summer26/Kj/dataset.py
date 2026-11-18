"""
Dataset for Kj: loads QCDZjj_herwig.h5 + VBFZ_herwig.h5.
Feedforward, no-constituents baseline for the Z_only_table (see summer26.md):
Kj sees only the leading jet's own four-momentum, no event-level kinematics,
no second jet, no constituents.

Jet slot 0 is randomized between the true leading and subleading jet at
production time (see makefriends.cpp / README.md "Jet slot randomisation") --
reading slot 0 here needs no extra code, the 50/50 swap is already baked in.

Jet features (4) -- the leading jet's own four-momentum:
  |jetEta|, jetM, jetPhi, jetPt  for jet 0 only

Labels: 0 = QCD background, 1 = VBF signal
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

JET_FEATURES = ["jetEta", "jetM", "jetPhi", "jetPt"]
N_FEAT = len(JET_FEATURES)  # 4

BKG_FILE = "../../friends/summer26/QCDZjj_herwig.h5"
SIG_FILE = "../../friends/summer26/VBFZ_herwig.h5"


def load_features(h5path):
    """Return x (N, 4) and kWeight (N,)."""
    with h5py.File(h5path, "r") as f:
        jets = np.stack([f[k][:, 0].astype(np.float32) for k in JET_FEATURES], axis=1)
        kw   = f["kWeight"][:].astype(np.float32)

    jets[:, 0] = np.abs(jets[:, 0])   # |jetEta|
    return jets, kw


class KjDataset(Dataset):
    def __init__(self, x, y, w=None):
        self.x = torch.from_numpy(x)
        self.y = torch.from_numpy(y).unsqueeze(1)
        self.w = torch.from_numpy(w) if w is not None else torch.ones(len(y))

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx], self.w[idx]
