"""
Dataset for KjB: loads QCDWjj_herwig.h5 + VBFW_herwig.h5.
Feedforward, no-constituents counterpart to NNjB (see summer26.md
W_only_table): same jet-level + boson scalars as NNjB, but no constituents --
isolates what NNjB's constituent branch adds on top of KjB's inputs.

Jet slot 0 is randomized between the true leading and subleading jet at
production time (see makefriends.cpp / README.md "Jet slot randomisation") --
reading slot 0 here needs no extra code, the 50/50 swap is already baked in.

Jet-level scalars (7 = 4 jet + 3 boson, leading jet = index 0, |jetEta| taken):
  [0]  |jetEta|  — absolute pseudorapidity
  [1]  jetM      — jet mass (GeV)
  [2]  jetPt     — transverse momentum (GeV)
  [3]  jetPhi    — azimuthal angle
  [4]  bosonPt   — generator boson pT (GeV)
  [5]  bosonEta  — generator boson η (signed)
  [6]  bosonPhi  — generator boson φ

Labels: 0 = QCD background, 1 = VBF signal
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

JET_FEATURES   = ['jetEta', 'jetM', 'jetPt', 'jetPhi']
BOSON_FEATURES = ['bosonPt', 'bosonEta', 'bosonPhi']
N_FEAT = len(JET_FEATURES) + len(BOSON_FEATURES)  # 7

BKG_FILE = '../../friends/summer26/QCDWjj_herwig.h5'
SIG_FILE = '../../friends/summer26/VBFW_herwig.h5'


def load_features(h5path):
    """Return x (N, 7) and kWeight (N,)."""
    with h5py.File(h5path, 'r') as f:
        jets  = {k: f[k][:, 0].astype(np.float32) for k in JET_FEATURES}
        boson = {k: f[k][:].astype(np.float32)     for k in BOSON_FEATURES}
        kw    = f['kWeight'][:].astype(np.float32)

    jets['jetEta'] = np.abs(jets['jetEta'])
    x = np.stack([jets[k] for k in JET_FEATURES]
               + [boson[k] for k in BOSON_FEATURES], axis=1)
    return x, kw


class KjBDataset(Dataset):
    def __init__(self, x, y, w=None):
        self.x = torch.from_numpy(x)
        self.y = torch.from_numpy(y).unsqueeze(1)
        self.w = torch.from_numpy(w) if w is not None else torch.ones(len(y))

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx], self.w[idx]
