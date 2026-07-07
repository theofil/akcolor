import torch
import torch.nn as nn
from dataset import N_JET_FEAT, N_CONSTIT_FEAT


class JetNN(nn.Module):
    """
    DeepSets-style network for VBF vs QCD classification (NNj variant).
    No high-level pull-vector jet scalars (no NC, |t|, θ_s); no constituent weight.

    phi: shared MLP applied to each constituent → masked sum pooling
    rho: MLP on pooled representation concatenated with jet-level scalars → logit

    Inputs:
      x_jet : (B, N_JET_FEAT=3)               — jet scalars: |η|, m, pT (normalised)
      x_jcs : (B, NC_MAX, N_CONSTIT_FEAT=3)   — constituent features (normalised)
      mask  : (B, NC_MAX)  bool               — True for valid (non-padded) constituents
    Output:
      (B, 1) raw logit
    """

    def __init__(self, dropout=0.3):
        super().__init__()
        self.phi = nn.Sequential(
            nn.Linear(N_CONSTIT_FEAT, 64), nn.ReLU(),
            nn.Linear(64, 64),             nn.ReLU(),
        )
        self.rho = nn.Sequential(
            nn.Linear(64 + N_JET_FEAT, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 64),              nn.BatchNorm1d(64),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

    def forward(self, x_jet, x_jcs, mask):
        B, N, F = x_jcs.shape
        h = self.phi(x_jcs.view(B * N, F)).view(B, N, -1)   # (B, NC_MAX, 64)
        h = (h * mask.unsqueeze(-1)).sum(dim=1)              # (B, 64) masked sum pool
        return self.rho(torch.cat([h, x_jet], dim=1))
