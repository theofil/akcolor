#!/usr/bin/env python3
"""
Constituent-level gradient attribution for DeepSetsNoScale.

Computes |d(logit)/d(x_i)| for each constituent i via autograd, then plots:
  - 2D hexbin in (dEta, dPhi) coloured by mean gradient magnitude,
    separately for signal and background jets
  - d(logit)/d(z) vs z (sensitivity to pT fraction)

Output: phi_gradient_maps.pdf
"""

import os
import sys
import pathlib

import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from interpret_utils import (
    load_model, load_test_set,
    apply_style, FIG_SIZE, LABEL_SIZE, LEGEND_SIZE, TICK_SIZE,
)

SUBSET    = 5000
GRAD_BATCH = 128
SEED      = 42


def compute_gradients(model, x_scaled_sub, mask_sub, device):
    """Returns |d(logit)/d(x_i)| of shape (N, 80, 3)."""
    grads = []
    N = len(x_scaled_sub)
    for s in range(0, N, GRAD_BATCH):
        xb = torch.from_numpy(x_scaled_sub[s:s + GRAD_BATCH]).to(device)
        mb = torch.from_numpy(mask_sub[s:s + GRAD_BATCH]).to(device)
        xb.requires_grad_(True)
        model(xb, mb).sum().backward()
        grads.append(xb.grad.detach().abs().cpu().numpy())
    return np.concatenate(grads, axis=0)   # (N, 80, 3)


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    model = load_model(device)
    model.eval()

    print('Loading test set ...')
    x_raw, x_scaled, mask, y_test = load_test_set()
    print(f'  {len(y_test)} jets  (sig={int(y_test.sum())}, bkg={int((y_test==0).sum())})')

    rng     = np.random.default_rng(SEED)
    idx_sub = rng.choice(len(y_test), size=SUBSET, replace=False)
    x_raw_sub   = x_raw[idx_sub]        # (5000, 80, 3) physical units
    x_scaled_sub = x_scaled[idx_sub]    # (5000, 80, 3) scaled
    mask_sub    = mask[idx_sub]          # (5000, 80)
    y_sub       = y_test[idx_sub]

    print(f'Computing gradients on {SUBSET} jets (batch={GRAD_BATCH}) ...')
    grad_all = compute_gradients(model, x_scaled_sub, mask_sub, device)
    # grad_all: (5000, 80, 3) — columns: [|d/d(dEta)|, |d/d(dPhi)|, |d/d(z)|]

    # flatten to constituents, mask padding
    mask_bool = mask_sub.astype(bool)
    dEta_flat = x_raw_sub[..., 0][mask_bool]
    dPhi_flat = x_raw_sub[..., 1][mask_bool]
    z_flat    = x_raw_sub[..., 2][mask_bool]
    y_flat    = np.repeat(y_sub, mask_sub.sum(axis=1).astype(int))

    # spatial gradient magnitude (dEta + dPhi channels)
    gmag_flat = (grad_all[..., 0][mask_bool] + grad_all[..., 1][mask_bool])
    # z-channel gradient
    gz_flat   = grad_all[..., 2][mask_bool]

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    for col, (label, sel) in enumerate([('signal', y_flat == 1),
                                         ('background', y_flat == 0)]):
        # ── top row: spatial hexbin ───────────────────────────────────────────
        ax = axes[0, col]
        hb = ax.hexbin(dEta_flat[sel], dPhi_flat[sel],
                       C=gmag_flat[sel],
                       reduce_C_function=np.mean,
                       gridsize=50, cmap='YlOrRd', mincnt=5)
        fig.colorbar(hb, ax=ax, label=r'mean $|{\partial}\,\mathrm{logit}/{\partial}\,\Delta\eta| + |{\partial}/{\partial}\,\Delta\phi|$',
                     fraction=0.046, pad=0.04)
        ax.set_xlabel(r'$\Delta\eta$', fontsize=LABEL_SIZE)
        ax.set_ylabel(r'$\Delta\phi$', fontsize=LABEL_SIZE)
        ax.set_title(f'Spatial attribution — {label}', fontsize=LABEL_SIZE)
        ax.tick_params(labelsize=TICK_SIZE)
        ax.set_aspect('equal')

        # ── bottom row: dlogit/dz vs z ────────────────────────────────────────
        ax = axes[1, col]
        ax.hexbin(z_flat[sel], gz_flat[sel],
                  gridsize=50, cmap='Blues', mincnt=5)
        ax.set_xlabel(r'$z = p_T^{\rm constit} / p_T^{\rm jet}$', fontsize=LABEL_SIZE)
        ax.set_ylabel(r'$|{\partial}\,\mathrm{logit}/{\partial}\,z|$', fontsize=LABEL_SIZE)
        ax.set_title(f'$p_T$ sensitivity — {label}', fontsize=LABEL_SIZE)
        ax.tick_params(labelsize=TICK_SIZE)
        ax.grid(True, linestyle=':', alpha=0.4)

    plt.tight_layout()
    out = HERE / 'phi_gradient_maps.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
