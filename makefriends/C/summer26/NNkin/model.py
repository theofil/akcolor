import torch
import torch.nn as nn
from dataset import N_FEAT


class KinNN(nn.Module):
    """
    Simple MLP for VBF vs QCD classification from kinematic variables only.

    Input:  (B, N_FEAT=12) — [dPhijj, dYjj, mjj, ptjj, |eta|0, m0, phi0, pt0,
                                                         |eta|1, m1, phi1, pt1]
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
