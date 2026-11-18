import torch
import torch.nn as nn
from dataset import N_FEAT


class KjjBjNN(nn.Module):
    """
    Simple MLP for VBF vs QCD (Z channel) classification using both leading
    jets' four-momenta, the 3rd pt-ordered jet's four-momentum (zero when
    absent), and generator-boson scalars -- no constituents. Feedforward,
    no-constituents counterpart to NNjjBj; used to isolate the value added
    by NNjjBj's constituent branches (see the Z_only_table in summer26.md).

    Input:  (B, N_FEAT=15) — jet0[eta,m,pt,phi], jet1[eta,m,pt,phi],
            jet2[eta,m,phi,pt], boson[pT,eta,phi]
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
