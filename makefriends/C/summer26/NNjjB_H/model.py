import torch
import torch.nn as nn
from dataset import N_JET_FEAT, N_BOSON_FEAT, N_CONSTIT_FEAT


class JetNN(nn.Module):
    """
    DeepSets-style network for VBF vs QCD classification (NNjjB variant).
    Uses raw low-level features for both the leading and sub-leading jets,
    plus the generator-boson 4-vector as event-level scalars. The phi MLP
    is shared across both jets' constituents.

    phi: shared MLP applied to each constituent → masked sum pooling per jet
    rho: MLP on [pool0, pool1, jet0_scalars, jet1_scalars, boson_scalars] → logit

    Inputs:
      x_jet0  : (B, N_JET_FEAT=3)               — jet 0 scalars: η, m, pT (normalised)
      x_jet1  : (B, N_JET_FEAT=3)               — jet 1 scalars: η, m, pT (normalised)
      x_boson : (B, N_BOSON_FEAT=4)             — boson pT, η, φ, M (normalised)
      x_jcs0  : (B, NC_MAX, N_CONSTIT_FEAT=3)   — jet 0 constituent features (normalised)
      x_jcs1  : (B, NC_MAX, N_CONSTIT_FEAT=3)   — jet 1 constituent features (normalised)
      mask0   : (B, NC_MAX)  bool               — True for valid constituents of jet 0
      mask1   : (B, NC_MAX)  bool               — True for valid constituents of jet 1
    Output:
      (B, 1) raw logit
    """

    def __init__(self, dropout=0.3):
        super().__init__()
        self.phi = nn.Sequential(
            nn.Linear(N_CONSTIT_FEAT, 64), nn.ReLU(),
            nn.Linear(64, 64),             nn.ReLU(),
        )
        # 64 (pool jet0) + 64 (pool jet1) + 3 (jet0) + 3 (jet1) + 4 (boson) = 138
        self.rho = nn.Sequential(
            nn.Linear(64 + 64 + N_JET_FEAT + N_JET_FEAT + N_BOSON_FEAT, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 64),                                                nn.BatchNorm1d(64),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

    def _pool(self, x_jcs, mask):
        B, N, F = x_jcs.shape
        h = self.phi(x_jcs.view(B * N, F)).view(B, N, -1)   # (B, NC_MAX, 64)
        return (h * mask.unsqueeze(-1)).sum(dim=1)            # (B, 64) masked sum pool

    def forward(self, x_jet0, x_jet1, x_boson, x_jcs0, x_jcs1, mask0, mask1):
        pool0 = self._pool(x_jcs0, mask0)   # (B, 64)
        pool1 = self._pool(x_jcs1, mask1)   # (B, 64)
        cat   = torch.cat([pool0, pool1, x_jet0, x_jet1, x_boson], dim=1)   # (B, 138)
        return self.rho(cat)
