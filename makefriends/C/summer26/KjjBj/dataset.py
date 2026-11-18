"""
Dataset for KjjBj: loads QCDZjj_herwig.h5 + VBFZ_herwig.h5.
Feedforward, no-constituents counterpart to NNjjBj (see summer26.md
Z_only_table): same jet-level/3rd-jet/boson scalars as NNjjBj, but no
constituents at all -- isolates what NNjjBj's constituent branches add on
top of KjjBj's inputs.

Jet features (4 per jet, both leading jets):
  [0]  jetEta — signed pseudorapidity
  [1]  jetM   — jet mass (GeV)
  [2]  jetPt  — transverse momentum (GeV)
  [3]  jetPhi — azimuthal angle

3rd-jet features (4, event-level; all zero when no 3rd jet — jetPt > 30 GeV
for real jets, so pT = 0 unambiguously encodes absence):
  [0]  jetEta — signed pseudorapidity
  [1]  jetM   — jet mass (GeV)
  [2]  jetPhi — azimuthal angle
  [3]  jetPt  — transverse momentum (GeV)

Boson features (3, event-level):
  [0]  bosonPt  — generator boson pT (GeV)
  [1]  bosonEta — generator boson η (signed)
  [2]  bosonPhi — generator boson φ

Labels: 0 = QCD background, 1 = VBF signal
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

JET_FEATURES   = ['jetEta', 'jetM', 'jetPt', 'jetPhi']
JET2_FEATURES  = ['jetEta', 'jetM', 'jetPhi', 'jetPt']
BOSON_FEATURES = ['bosonPt', 'bosonEta', 'bosonPhi']
N_FEAT = len(JET_FEATURES) * 2 + len(JET2_FEATURES) + len(BOSON_FEATURES)  # 15

BKG_FILE = '../../friends/summer26/QCDZjj_herwig.h5'
SIG_FILE = '../../friends/summer26/VBFZ_herwig.h5'


def load_features(h5path):
    """Return x (N, 15) and kWeight (N,)."""
    with h5py.File(h5path, 'r') as f:
        jets0 = {k: f[k][:, 0].astype(np.float32) for k in JET_FEATURES}
        jets1 = {k: f[k][:, 1].astype(np.float32) for k in JET_FEATURES}
        jets2 = {k: f[k][:, 2].astype(np.float32) for k in JET2_FEATURES}
        boson = {k: f[k][:].astype(np.float32)    for k in BOSON_FEATURES}
        kw    = f['kWeight'][:].astype(np.float32)

    x = np.concatenate([
        np.stack([jets0[k] for k in JET_FEATURES], axis=1),
        np.stack([jets1[k] for k in JET_FEATURES], axis=1),
        np.stack([jets2[k] for k in JET2_FEATURES], axis=1),
        np.stack([boson[k] for k in BOSON_FEATURES], axis=1),
    ], axis=1)
    return x, kw


class KjjBjDataset(Dataset):
    def __init__(self, x, y, w=None):
        self.x = torch.from_numpy(x)
        self.y = torch.from_numpy(y).unsqueeze(1)
        self.w = torch.from_numpy(w) if w is not None else torch.ones(len(y))

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx], self.w[idx]
