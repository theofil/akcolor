import torch
import torch.nn as nn
from dataset import N_FEAT


class KjNN(nn.Module):
    """
    Simple MLP for VBF vs QCD (Z channel) classification using only the leading
    jet's own four-momentum -- no event-level kinematics, no second jet, no
    constituents. Used as the no-constituents baseline in the Z_only_table
    (see summer26.md) to isolate the value added by jet substructure.

    Input:  (B, N_FEAT=4) — [|eta|, m, phi, pt] of jet 0
    Output: (B, 1) raw logit
    """

    def __init__(self, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(N_FEAT, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 64),    nn.BatchNorm1d(64),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 32),     nn.BatchNorm1d(32),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(32, 1),
        )

    def forward(self, x):
        return self.net(x)
