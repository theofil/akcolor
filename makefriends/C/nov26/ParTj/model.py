"""ParTj: Particle-Transformer-style jet tagger with pairwise interaction features.

Follows ParT (Qu, Gouskos, Kagan, arXiv:2202.03772), scaled down for ~1M training jets:
  particle tokens  MLP 9 -> 128 -> D on [jcsDEta (signed), jcsDPhi, jcsPt/jetPt] standardised,
                   PF type one-hot (5) and charge (the NNjPF inputs)
  pair features    for every pair of constituents i != j, from their 4-vectors:
                   ln kT, ln z, ln dR, ln m^2  (dR in rapidity and periodic phi,
                   kT = min(pT_i, pT_j) dR, z = min(pT) / (pT_i + pT_j), m^2 = (p_i + p_j)^2),
                   BatchNorm over valid pairs, MLP 4 -> h -> h -> H: one bias per attention
                   head, added to the attention logits of every particle block
  particle blocks  L x pre-LayerNorm multi-head self-attention (+ pair bias, padded keys
                   masked) and a GELU feed-forward, with residual connections
  class blocks     C x a learnable class token attending to [class token, particles]
  head             [LayerNorm(class token), jet scalars |eta|, m, pT] -> 64 -> 1
"""
import math

import torch
import torch.nn as nn
import torch.nn.functional as F

N_JET_FEAT, N_CONT_FEAT, N_PF_TYPES = 3, 3, 5
N_CONSTIT_IN = N_CONT_FEAT + N_PF_TYPES + 1          # 9
NEG = -1e9                                           # attention logit of a padded key


def constituent_inputs(x_jcs, jcs_pf, jcs_q):
    """(B, N, 9): continuous features, PF one-hot (type 0 = empty -> all zeros), charge."""
    onehot = F.one_hot(jcs_pf.long(), N_PF_TYPES + 1)[..., 1:].to(x_jcs.dtype)
    return torch.cat([x_jcs, onehot, jcs_q.to(x_jcs.dtype).unsqueeze(-1)], dim=-1)


def pair_features(p4, eps=1e-8):
    """(B, N, N, 4) = ln kT, ln z, ln dR, ln m^2 for every pair, from p4 = (px, py, pz, E)."""
    px, py, pz, e = p4.unbind(-1)
    pt = torch.sqrt(px * px + py * py)
    rap = 0.5 * torch.log((e + pz).clamp(min=eps) / (e - pz).clamp(min=eps))
    phi = torch.atan2(py, px)
    drap = rap.unsqueeze(2) - rap.unsqueeze(1)
    dphi = torch.remainder(phi.unsqueeze(2) - phi.unsqueeze(1) + math.pi, 2 * math.pi) - math.pi
    delta = torch.sqrt(drap * drap + dphi * dphi)
    pti, ptj = pt.unsqueeze(2), pt.unsqueeze(1)
    ptmin = torch.minimum(pti, ptj)
    lnkt = torch.log((ptmin * delta).clamp(min=eps))
    lnz = torch.log((ptmin / (pti + ptj).clamp(min=eps)).clamp(min=eps))
    lndelta = torch.log(delta.clamp(min=eps))
    s = p4.unsqueeze(2) + p4.unsqueeze(1)
    m2 = s[..., 3] ** 2 - s[..., 0] ** 2 - s[..., 1] ** 2 - s[..., 2] ** 2
    lnm2 = torch.log(m2.clamp(min=eps))
    return torch.stack([lnkt, lnz, lndelta, lnm2], dim=-1)


class PairEmbed(nn.Module):
    """Pair features -> one attention bias per head, (B, H, N, N); zero for invalid pairs."""

    def __init__(self, num_heads, hidden=32, n_feat=4):
        super().__init__()
        self.num_heads = num_heads
        self.bn = nn.BatchNorm1d(n_feat)
        self.mlp = nn.Sequential(nn.Linear(n_feat, hidden), nn.GELU(),
                                 nn.Linear(hidden, hidden), nn.GELU(),
                                 nn.Linear(hidden, num_heads))

    def forward(self, p4, mask):
        b, n, _ = p4.shape
        eye = torch.eye(n, dtype=torch.bool, device=p4.device)
        valid = mask.unsqueeze(2) & mask.unsqueeze(1) & ~eye          # no padded, no self pairs
        out = p4.new_zeros(b, n, n, self.num_heads)
        if valid.any():
            out[valid] = self.mlp(self.bn(pair_features(p4)[valid]))
        return out.permute(0, 3, 1, 2)


class ParticleBlock(nn.Module):
    def __init__(self, dim, heads, ffn_ratio, dropout):
        super().__init__()
        self.ln1, self.ln2 = nn.LayerNorm(dim), nn.LayerNorm(dim)
        self.attn = nn.MultiheadAttention(dim, heads, dropout=dropout, batch_first=True)
        self.ffn = nn.Sequential(nn.Linear(dim, ffn_ratio * dim), nn.GELU(), nn.Dropout(dropout),
                                 nn.Linear(ffn_ratio * dim, dim))
        self.drop = nn.Dropout(dropout)

    def forward(self, x, bias):                    # bias (B*H, N, N): pair bias + padded keys
        h = self.ln1(x)
        x = x + self.drop(self.attn(h, h, h, attn_mask=bias, need_weights=False)[0])
        return x + self.drop(self.ffn(self.ln2(x)))


class ClassBlock(nn.Module):
    def __init__(self, dim, heads, ffn_ratio, dropout):
        super().__init__()
        self.ln1, self.ln2 = nn.LayerNorm(dim), nn.LayerNorm(dim)
        self.attn = nn.MultiheadAttention(dim, heads, dropout=dropout, batch_first=True)
        self.ffn = nn.Sequential(nn.Linear(dim, ffn_ratio * dim), nn.GELU(), nn.Dropout(dropout),
                                 nn.Linear(ffn_ratio * dim, dim))
        self.drop = nn.Dropout(dropout)

    def forward(self, x, cls, key_pad):            # cls (B, 1, D); key_pad (B, 1+N), True = ignore
        u = self.ln1(torch.cat([cls, x], dim=1))
        a = self.attn(u[:, :1], u, u, key_padding_mask=key_pad, need_weights=False)[0]
        cls = cls + self.drop(a)
        return cls + self.drop(self.ffn(self.ln2(cls)))


class ParticleTransformerJet(nn.Module):
    """Inputs: x_jet (B, 3), x_jcs (B, N, 3), jcs_pf (B, N) int, jcs_q (B, N),
    jcs_p4 (B, N, 4) in GeV, mask (B, N) bool. Output: (B, 1) raw logit."""

    def __init__(self, embed_dim=64, num_heads=4, num_layers=3, num_cls_layers=2, ffn_ratio=4,
                 pair_hidden=32, dropout=0.1, head_dropout=0.3):
        super().__init__()
        self.num_heads = num_heads
        self.embed = nn.Sequential(nn.Linear(N_CONSTIT_IN, 128), nn.GELU(), nn.Linear(128, embed_dim))
        self.pair = PairEmbed(num_heads, pair_hidden)
        self.blocks = nn.ModuleList([ParticleBlock(embed_dim, num_heads, ffn_ratio, dropout)
                                     for _ in range(num_layers)])
        self.cls_token = nn.Parameter(torch.zeros(1, 1, embed_dim))
        nn.init.trunc_normal_(self.cls_token, std=0.02)
        self.cls_blocks = nn.ModuleList([ClassBlock(embed_dim, num_heads, ffn_ratio, dropout)
                                         for _ in range(num_cls_layers)])
        self.norm = nn.LayerNorm(embed_dim)
        self.head = nn.Sequential(nn.Linear(embed_dim + N_JET_FEAT, 64), nn.GELU(),
                                  nn.Dropout(head_dropout), nn.Linear(64, 1))

    def forward(self, x_jet, x_jcs, jcs_pf, jcs_q, jcs_p4, mask):
        b, n = mask.shape
        x = self.embed(constituent_inputs(x_jcs, jcs_pf, jcs_q))
        bias = self.pair(jcs_p4, mask).masked_fill(~mask[:, None, None, :], NEG)
        bias = bias.reshape(b * self.num_heads, n, n)
        for blk in self.blocks:
            x = blk(x, bias)
        key_pad = torch.cat([mask.new_zeros(b, 1), ~mask], dim=1)
        cls = self.cls_token.expand(b, -1, -1)
        for blk in self.cls_blocks:
            cls = blk(x, cls, key_pad)
        return self.head(torch.cat([self.norm(cls.squeeze(1)), x_jet], dim=1))
