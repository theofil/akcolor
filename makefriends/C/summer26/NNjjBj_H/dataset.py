"""
DijetDataset for NNjjBj: loads QCDHjj_herwig.h5 + VBFH_herwig.h5.
Identical inputs to NNjjB (both leading jets, raw low-level features, plus the
generator-boson 4-vector as event-level scalars) extended with the 3rd
pt-ordered jet — its 4-momentum and its constituents — when a 3rd jet exists.

Jet features (3 per jet, both leading jets):
  [0]  jetEta — signed pseudorapidity
  [1]  jetM   — jet mass (GeV)
  [2]  jetPt  — transverse momentum (GeV)

3rd-jet features (4, event-level; all zero when no 3rd jet — jetPt > 30 GeV
for real jets, so pT = 0 unambiguously encodes absence):
  [0]  jetEta — signed pseudorapidity
  [1]  jetM   — jet mass (GeV)
  [2]  jetPhi — azimuthal angle
  [3]  jetPt  — transverse momentum (GeV)

Boson features (4, event-level):
  [0]  bosonPt  — generator boson pT (GeV)
  [1]  bosonEta — generator boson η (signed)
  [2]  bosonPhi — generator boson φ
  [3]  bosonM   — generator boson mass (GeV)

Constituent features per jet, all three jets (NC_MAX=80, 3 per constituent;
jet-2 slots are all zero when no 3rd jet exists, giving an all-False mask and
a zero pooled vector):
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
JET2_FEATURES    = ['jetEta', 'jetM', 'jetPhi', 'jetPt']
BOSON_FEATURES   = ['bosonPt', 'bosonEta', 'bosonPhi', 'bosonM']
CONSTIT_FEATURES = ['jcsDEtaRaw', 'jcsDPhi', 'jcsPt']

N_JET_FEAT     = len(JET_FEATURES)      # 3
N_JET2_FEAT    = len(JET2_FEATURES)     # 4
N_BOSON_FEAT   = len(BOSON_FEATURES)    # 4
N_CONSTIT_FEAT = len(CONSTIT_FEATURES)  # 3
NC_MAX         = 80

BKG_FILE = '../../friends/summer26/QCDHjj_herwig.h5'
SIG_FILE = '../../friends/summer26/VBFH_herwig.h5'


def load_features(h5path):
    """Return x_jet0 (N,3), x_jet1 (N,3), x_jet2 (N,4), x_boson (N,4),
    x_jcs0 (N,80,3), x_jcs1 (N,80,3), x_jcs2 (N,80,3), kWeight (N,)."""
    with h5py.File(h5path, 'r') as f:
        jets0 = {k: f[k][:, 0].astype(np.float32)       for k in JET_FEATURES}
        jets1 = {k: f[k][:, 1].astype(np.float32)       for k in JET_FEATURES}
        jets2 = {k: f[k][:, 2].astype(np.float32)       for k in JET2_FEATURES}
        boson = {k: f[k][:].astype(np.float32)          for k in BOSON_FEATURES}
        jcs0  = {k: f[k][:, 0, :].astype(np.float32)    for k in CONSTIT_FEATURES}
        jcs1  = {k: f[k][:, 1, :].astype(np.float32)    for k in CONSTIT_FEATURES}
        jcs2  = {k: f[k][:, 2, :].astype(np.float32)    for k in CONSTIT_FEATURES}
        kw    = f['kWeight'][:].astype(np.float32)

    x_jet0 = np.stack([jets0[k] for k in JET_FEATURES], axis=1)   # (N, 3)
    x_jet1 = np.stack([jets1[k] for k in JET_FEATURES], axis=1)   # (N, 3)
    x_jet2 = np.stack([jets2[k] for k in JET2_FEATURES], axis=1)  # (N, 4)
    x_bos  = np.stack([boson[k] for k in BOSON_FEATURES], axis=1)  # (N, 4)
    x_jcs0 = np.stack([jcs0[k]  for k in CONSTIT_FEATURES], axis=2)  # (N, 80, 3)
    x_jcs1 = np.stack([jcs1[k]  for k in CONSTIT_FEATURES], axis=2)  # (N, 80, 3)
    x_jcs2 = np.stack([jcs2[k]  for k in CONSTIT_FEATURES], axis=2)  # (N, 80, 3)

    x_jcs0[:, :, 2] /= x_jet0[:, 2:3]   # jcsPt → jcsPt / jetPt  (jet 0)
    x_jcs1[:, :, 2] /= x_jet1[:, 2:3]   # jcsPt → jcsPt / jetPt  (jet 1)
    # jet 2: guard the division — jetPt2 (column 3) is 0 when no 3rd jet
    # exists and the constituent slots are 0 too (0/0 would give NaN)
    x_jcs2[:, :, 2] /= np.where(x_jet2[:, 3:4] > 0, x_jet2[:, 3:4], 1.0)

    return x_jet0, x_jet1, x_jet2, x_bos, x_jcs0, x_jcs1, x_jcs2, kw


class DijetDataset(Dataset):
    def __init__(self, x_jet0, x_jet1, x_jet2, x_boson, x_jcs0, x_jcs1, x_jcs2,
                 mask0, mask1, mask2, y, w=None):
        self.x_jet0  = torch.from_numpy(x_jet0)
        self.x_jet1  = torch.from_numpy(x_jet1)
        self.x_jet2  = torch.from_numpy(x_jet2)
        self.x_boson = torch.from_numpy(x_boson)
        self.x_jcs0  = torch.from_numpy(x_jcs0)
        self.x_jcs1  = torch.from_numpy(x_jcs1)
        self.x_jcs2  = torch.from_numpy(x_jcs2)
        self.mask0   = torch.from_numpy(mask0)
        self.mask1   = torch.from_numpy(mask1)
        self.mask2   = torch.from_numpy(mask2)
        self.y       = torch.from_numpy(y).unsqueeze(1)
        self.w       = torch.from_numpy(w) if w is not None else torch.ones(len(y))

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return (self.x_jet0[idx], self.x_jet1[idx], self.x_jet2[idx], self.x_boson[idx],
                self.x_jcs0[idx], self.x_jcs1[idx], self.x_jcs2[idx],
                self.mask0[idx],  self.mask1[idx],  self.mask2[idx],
                self.y[idx],      self.w[idx])
