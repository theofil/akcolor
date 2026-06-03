import torch
import torch.nn as nn
from dataset import N_JET_FEAT, N_CONSTIT_FEAT

D_MODEL     = 64
N_HEADS     = 4
N_LAYERS    = 2
DIM_FF      = 128


class TransJetNN(nn.Module):
    """
    Transformer-based constituent network for VBF vs QCD classification.

    Constituents are embedded, processed by a Transformer encoder with
    padding-aware self-attention, then mean-pooled over valid constituents.
    The pooled representation is concatenated with jet-level scalars and
    passed to an MLP classifier.

    Inputs:
      x_jet : (B, N_JET_FEAT)              — jet-level scalars (normalised)
      x_jcs : (B, NC_MAX, N_CONSTIT_FEAT)  — constituent features (normalised)
      mask  : (B, NC_MAX)  bool            — True for valid (non-padded) constituents
    Output:
      (B, 1) raw logit
    """

    def __init__(self, dropout=0.1, clf_dropout=0.3):
        super().__init__()

        self.embed = nn.Sequential(
            nn.Linear(N_CONSTIT_FEAT, D_MODEL),
            nn.LayerNorm(D_MODEL),
        )

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=D_MODEL,
            nhead=N_HEADS,
            dim_feedforward=DIM_FF,
            dropout=dropout,
            batch_first=True,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=N_LAYERS)

        self.classifier = nn.Sequential(
            nn.Linear(D_MODEL + N_JET_FEAT, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(clf_dropout),
            nn.Linear(128, 64),                   nn.BatchNorm1d(64),  nn.ReLU(), nn.Dropout(clf_dropout),
            nn.Linear(64, 1),
        )

    def forward(self, x_jet, x_jcs, mask):
        # mask: True = valid constituent  →  pad_mask: True = ignore (PyTorch convention)
        pad_mask = ~mask                                    # (B, NC_MAX)

        h = self.embed(x_jcs)                              # (B, NC_MAX, D_MODEL)
        h = self.encoder(h, src_key_padding_mask=pad_mask) # (B, NC_MAX, D_MODEL)

        # masked mean-pool over valid constituents
        valid = mask.unsqueeze(-1).float()                 # (B, NC_MAX, 1)
        pooled = (h * valid).sum(dim=1) / valid.sum(dim=1).clamp(min=1)  # (B, D_MODEL)

        return self.classifier(torch.cat([pooled, x_jet], dim=1))
