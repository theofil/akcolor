import torch
import torch.nn as nn

N_CONSTIT_FEAT = 5


class DeepSetsNoJet(nn.Module):
    """
    DeepSets classifier using constituent features only (no jet-level scalars).
    Identical phi/rho capacity to NNpull but x_jet is never seen.
    """

    def __init__(self, dropout=0.3):
        super().__init__()
        self.phi = nn.Sequential(
            nn.Linear(N_CONSTIT_FEAT, 64), nn.ReLU(),
            nn.Linear(64, 64),             nn.ReLU(),
        )
        self.rho = nn.Sequential(
            nn.Linear(64, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 64), nn.BatchNorm1d(64),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

    def forward(self, x_jcs, mask):
        B, N, F = x_jcs.shape
        h = self.phi(x_jcs.view(B * N, F)).view(B, N, -1)
        h = (h * mask.unsqueeze(-1)).sum(dim=1)   # masked sum pool
        return self.rho(h)
