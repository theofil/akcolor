import torch
import torch.nn as nn
from dataset import N_FEAT


class KjBNN(nn.Module):
    """
    Simple MLP for VBF vs QCD (W channel) classification using the leading
    jet's own four-momentum plus generator-boson scalars -- no constituents.
    Feedforward, no-constituents counterpart to NNjB; used to isolate the
    value added by NNjB's constituent branch (see the W_only_table in
    summer26.md).

    Input:  (B, N_FEAT=7) — [|eta|, m, pt, phi] of jet 0 + [bosonPt, bosonEta, bosonPhi]
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
