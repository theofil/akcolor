"""
DijetDataset for NNpull: loads QCDZjj_herwig.h5 + VBFZ_herwig.h5.
Full pull-vector information kept except NC (constituent count removed).

Jet features (5, leading jet = index 0, |jetEta| taken):
  [0]  |jetEta| — absolute pseudorapidity
  [1]  jetM     — jet mass (GeV)
  [2]  jetPVM   — pull vector magnitude |t⃗|
  [3]  jetPt    — transverse momentum (GeV)
  [4]  jetSPVA  — signed pull vector angle θ_s

Constituent features per jet (NC_MAX=80, 4 per constituent):
  [0]  jcsDEta — sign-flipped Δη relative to jet axis
  [1]  jcsDPhi — ΔΦ relative to jet axis
  [2]  jcsPt   — constituent pT / jet pT (dimensionless)  [used as validity mask: > 0]
  [3]  jcsW    — pull-vector weight

Labels: 0 = QCD background, 1 = VBF signal
"""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset

JET_FEATURES     = ['jetEta', 'jetM', 'jetPVM', 'jetPt', 'jetSPVA']
CONSTIT_FEATURES = ['jcsDEta', 'jcsDPhi', 'jcsPt', 'jcsW']

N_JET_FEAT     = len(JET_FEATURES)      # 5
N_CONSTIT_FEAT = len(CONSTIT_FEATURES)  # 4
NC_MAX         = 80

BKG_FILE = '../../friends/summer26/QCDZjj_herwig.h5'
SIG_FILE = '../../friends/summer26/VBFZ_herwig.h5'


def load_features(h5path):
    """Return x_jet (N, 5), x_jcs (N, NC_MAX, 4), and kWeight (N,)."""
    with h5py.File(h5path, 'r') as f:
        jets = {k: f[k][:, 0].astype(np.float32)       for k in JET_FEATURES}
        jcs  = {k: f[k][:, 0, :].astype(np.float32)    for k in CONSTIT_FEATURES}
        kw   = f['kWeight'][:].astype(np.float32)

    jets['jetEta'] = np.abs(jets['jetEta'])
    x_jet = np.stack([jets[k] for k in JET_FEATURES],     axis=1)  # (N, 5)
    x_jcs = np.stack([jcs[k]  for k in CONSTIT_FEATURES], axis=2)  # (N, 80, 4)
    x_jcs[:, :, 2] /= x_jet[:, 3:4]   # jcsPt → jcsPt / jetPt
    return x_jet, x_jcs, kw


class DijetDataset(Dataset):
    def __init__(self, x_jet, x_jcs, mask_jcs, y, w=None):
        self.x_jet = torch.from_numpy(x_jet)
        self.x_jcs = torch.from_numpy(x_jcs)
        self.mask  = torch.from_numpy(mask_jcs)   # (N, NC_MAX), bool
        self.y     = torch.from_numpy(y).unsqueeze(1)
        self.w     = torch.from_numpy(w) if w is not None else torch.ones(len(y))

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.x_jet[idx], self.x_jcs[idx], self.mask[idx], self.y[idx], self.w[idx]
