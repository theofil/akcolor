"""NNjPF: the NNj DeepSets jet tagger with particle type and charge per constituent.

Per constituent (9 inputs): [jcsDEta (signed), jcsDPhi, jcsPt/jetPt] standardised, the PF type
one-hot (h+-, e, mu, gamma, h0; built here, so empty slots give zeros) and the charge.
phi: 9 -> 64 -> 64 (ReLU), masked sum over constituents.
rho: [pooled 64, jet scalars |eta|, m, pT (standardised)] -> 128 -> 64 -> 1 (BatchNorm, ReLU,
Dropout), exactly as summer26/NNj/model.py.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F

N_JET_FEAT, N_CONT_FEAT, N_PF_TYPES = 3, 3, 5
N_CONSTIT_IN = N_CONT_FEAT + N_PF_TYPES + 1          # 9


def constituent_inputs(x_jcs, jcs_pf, jcs_q):
    """(B, N, 9): continuous features, PF one-hot (type 0 = empty -> all zeros), charge."""
    onehot = F.one_hot(jcs_pf.long(), N_PF_TYPES + 1)[..., 1:].to(x_jcs.dtype)
    return torch.cat([x_jcs, onehot, jcs_q.to(x_jcs.dtype).unsqueeze(-1)], dim=-1)


class JetPFNN(nn.Module):
    """Inputs: x_jet (B, 3), x_jcs (B, N, 3), jcs_pf (B, N) int, jcs_q (B, N), mask (B, N) bool.
    Output: (B, 1) raw logit (VBF jet vs QCD jet)."""

    def __init__(self, dropout=0.3):
        super().__init__()
        self.phi = nn.Sequential(
            nn.Linear(N_CONSTIT_IN, 64), nn.ReLU(),
            nn.Linear(64, 64), nn.ReLU(),
        )
        self.rho = nn.Sequential(
            nn.Linear(64 + N_JET_FEAT, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 64), nn.BatchNorm1d(64), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

    def forward(self, x_jet, x_jcs, jcs_pf, jcs_q, mask):
        h = self.phi(constituent_inputs(x_jcs, jcs_pf, jcs_q))     # (B, N, 64)
        h = (h * mask.unsqueeze(-1)).sum(dim=1)                    # masked sum pool
        return self.rho(torch.cat([h, x_jet], dim=1))
