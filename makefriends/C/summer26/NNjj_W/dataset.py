"""
DijetDataset for NNjj: loads QCDWjj_herwig.h5 + VBFW_herwig.h5.
Raw low-level features for both the leading (index 0)
and sub-leading (index 1) jets.

Jet features (3 per jet, both jets):
  [0]  jetEta — signed pseudorapidity
  [1]  jetM   — jet mass (GeV)
  [2]  jetPt  — transverse momentum (GeV)

Constituent features per jet (NC_MAX=80, 3 per constituent):
  [0]  jcsDEtaRaw — raw Δη relative to jet axis
  [1]  jcsDPhi    — ΔΦ relative to jet axis
  [2]  jcsPt      — constituent pT / jet pT (dimensionless)  [used as validity mask: > 0]

Labels: 0 = QCD background, 1 = VBF signal
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

JET_FEATURES     = ['jetEta', 'jetM', 'jetPt']
CONSTIT_FEATURES = ['jcsDEtaRaw', 'jcsDPhi', 'jcsPt']

N_JET_FEAT     = len(JET_FEATURES)      # 3
N_CONSTIT_FEAT = len(CONSTIT_FEATURES)  # 3
NC_MAX         = 80

BKG_FILE = '../../friends/summer26/QCDWjj_herwig.h5'
SIG_FILE = '../../friends/summer26/VBFW_herwig.h5'


def load_features(h5path):
    """Return x_jet0 (N,3), x_jet1 (N,3), x_jcs0 (N,80,3), x_jcs1 (N,80,3), kWeight (N,)."""
    with h5py.File(h5path, 'r') as f:
        jets0 = {k: f[k][:, 0].astype(np.float32)       for k in JET_FEATURES}
        jets1 = {k: f[k][:, 1].astype(np.float32)       for k in JET_FEATURES}
        jcs0  = {k: f[k][:, 0, :].astype(np.float32)    for k in CONSTIT_FEATURES}
        jcs1  = {k: f[k][:, 1, :].astype(np.float32)    for k in CONSTIT_FEATURES}
        kw    = f['kWeight'][:].astype(np.float32)

    x_jet0 = np.stack([jets0[k] for k in JET_FEATURES], axis=1)   # (N, 3)
    x_jet1 = np.stack([jets1[k] for k in JET_FEATURES], axis=1)   # (N, 3)
    x_jcs0 = np.stack([jcs0[k]  for k in CONSTIT_FEATURES], axis=2)  # (N, 80, 3)
    x_jcs1 = np.stack([jcs1[k]  for k in CONSTIT_FEATURES], axis=2)  # (N, 80, 3)

    x_jcs0[:, :, 2] /= x_jet0[:, 2:3]   # jcsPt → jcsPt / jetPt  (jet 0)
    x_jcs1[:, :, 2] /= x_jet1[:, 2:3]   # jcsPt → jcsPt / jetPt  (jet 1)

    return x_jet0, x_jet1, x_jcs0, x_jcs1, kw


class DijetDataset(Dataset):
    def __init__(self, x_jet0, x_jet1, x_jcs0, x_jcs1, mask0, mask1, y, w=None):
        self.x_jet0 = torch.from_numpy(x_jet0)
        self.x_jet1 = torch.from_numpy(x_jet1)
        self.x_jcs0 = torch.from_numpy(x_jcs0)
        self.x_jcs1 = torch.from_numpy(x_jcs1)
        self.mask0  = torch.from_numpy(mask0)
        self.mask1  = torch.from_numpy(mask1)
        self.y      = torch.from_numpy(y).unsqueeze(1)
        self.w      = torch.from_numpy(w) if w is not None else torch.ones(len(y))

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return (self.x_jet0[idx], self.x_jet1[idx],
                self.x_jcs0[idx], self.x_jcs1[idx],
                self.mask0[idx],  self.mask1[idx],
                self.y[idx],      self.w[idx])
