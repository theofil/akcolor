import torch
import torch.nn as nn

N_CONSTIT_FEAT = 5
D_MODEL        = 64
N_HEADS        = 4
N_LAYERS       = 2
DIM_FF         = 128


class TransNoJet(nn.Module):
    """
    Transformer constituent network without jet-level scalars.
    Identical transformer + classifier capacity to NNtrans but x_jet is never seen.
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
            nn.Linear(D_MODEL, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(clf_dropout),
            nn.Linear(128, 64),      nn.BatchNorm1d(64),  nn.ReLU(), nn.Dropout(clf_dropout),
            nn.Linear(64, 1),
        )

    def forward(self, x_jcs, mask):
        pad_mask = ~mask                                     # True = ignore (PyTorch convention)
        h = self.embed(x_jcs)
        h = self.encoder(h, src_key_padding_mask=pad_mask)
        valid  = mask.unsqueeze(-1).float()
        pooled = (h * valid).sum(dim=1) / valid.sum(dim=1).clamp(min=1)
        return self.classifier(pooled)
