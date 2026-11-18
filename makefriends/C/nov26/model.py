"""DeepSets jet tagger of the nov26 nets: summer26/NNj/model.py with a variable number of jet scalars.

phi: per constituent [jcsDEta (signed), jcsDPhi, jcsPt/jetPt] (standardised) -> 64 -> 64 (ReLU),
     masked sum over constituents.
rho: [pooled 64, jet scalars] -> 128 -> 64 -> 1 (BatchNorm, ReLU, Dropout).
The jet scalars are |eta|, m, pT, followed by the net's beam-connection inputs (common.NETS):
3 for NNj, 8 for NNjBC.
"""
import torch
import torch.nn as nn

N_CONSTIT_FEAT = 3


class JetNN(nn.Module):
    """Inputs: x_jet (B, n_jet_feat), x_jcs (B, N, 3), mask (B, N) bool. Output: (B, 1) raw logit."""

    def __init__(self, n_jet_feat, dropout=0.3):
        super().__init__()
        self.phi = nn.Sequential(
            nn.Linear(N_CONSTIT_FEAT, 64), nn.ReLU(),
            nn.Linear(64, 64), nn.ReLU(),
        )
        self.rho = nn.Sequential(
            nn.Linear(64 + n_jet_feat, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 64), nn.BatchNorm1d(64), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

    def forward(self, x_jet, x_jcs, mask):
        h = self.phi(x_jcs)                                       # (B, N, 64)
        h = (h * mask.unsqueeze(-1)).sum(dim=1)                   # masked sum pool
        return self.rho(torch.cat([h, x_jet], dim=1))
