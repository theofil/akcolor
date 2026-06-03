#!/usr/bin/env python3
"""
Visualise the effective per-constituent kernel of DeepSetsNoScale.

The effective contribution of constituent i to the logit is  w · phi(x_i)
where  w = E[d rho(h)/d h]  (mean Jacobian of rho w.r.t. its input, shape 64).

Outputs:
  phi_kernel_heatmap.pdf   — 2-page PDF:
      page 1: 2×2 panel  w·phi(dEta, dPhi) at four fixed z values
      page 2: radial profiles  w·phi vs r  for each z
"""

import os
import sys
import pathlib

import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from interpret_utils import (
    load_model, load_scaler, load_test_set, scale_jcs,
    compute_phi_aggregates, compute_mean_rho_jacobian,
    apply_style, FIG_SIZE, LABEL_SIZE, LEGEND_SIZE, TICK_SIZE,
)


GRID_N  = 60
ETA_LIM = 0.6
Z_VALS  = [0.05, 0.10, 0.20, 0.50]
Z_COLS  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']


def phi_kernel_grid(model, scaler, w_t, z_fixed, device):
    """Evaluate w·phi on a GRID_N×GRID_N (dEta, dPhi) grid at fixed z."""
    eta = np.linspace(-ETA_LIM, ETA_LIM, GRID_N, dtype=np.float32)
    dEta_g, dPhi_g = np.meshgrid(eta, eta)              # (N, N)
    z_col = np.full(GRID_N * GRID_N, z_fixed, dtype=np.float32)
    pts   = np.column_stack([dEta_g.ravel(), dPhi_g.ravel(), z_col])  # (N*N, 3)
    pts_s = scaler.transform(pts).astype(np.float32)
    with torch.no_grad():
        phi_out = model.phi(torch.from_numpy(pts_s).to(device))   # (N*N, 64)
        f_vals  = (phi_out * w_t).sum(-1).cpu().numpy()            # (N*N,)
    return dEta_g, dPhi_g, f_vals.reshape(GRID_N, GRID_N)


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    model  = load_model(device)
    scaler = load_scaler()

    print('Loading test set ...')
    x_raw, x_scaled, mask, y_test = load_test_set()
    print(f'  {len(y_test)} jets  (sig={int(y_test.sum())}, bkg={int((y_test==0).sum())})')

    print('Computing phi aggregates ...')
    h_test = compute_phi_aggregates(model, x_scaled, mask, device)

    print('Computing mean rho Jacobian w ...')
    w = compute_mean_rho_jacobian(model, h_test, device)
    w_t = torch.from_numpy(w).to(device)
    print(f'  |w| = {np.linalg.norm(w):.4f}')

    out_pdf = HERE / 'phi_kernel_heatmap.pdf'
    with PdfPages(out_pdf) as pdf:

        # ── page 1: 2×2 heatmaps ─────────────────────────────────────────────
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        for ax, z, col in zip(axes.ravel(), Z_VALS, Z_COLS):
            dEta_g, dPhi_g, f2d = phi_kernel_grid(model, scaler, w_t, z, device)
            vmax = np.abs(f2d).max()
            vmax = vmax if vmax > 0 else 1.0
            im = ax.pcolormesh(dEta_g, dPhi_g, f2d,
                               cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                               shading='auto')
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            ax.set_title(rf'$z = {z}$', fontsize=LABEL_SIZE)
            ax.set_xlabel(r'$\Delta\eta$', fontsize=LABEL_SIZE)
            ax.set_ylabel(r'$\Delta\phi$', fontsize=LABEL_SIZE)
            ax.tick_params(labelsize=TICK_SIZE)
            ax.set_aspect('equal')
        fig.suptitle(r'$w \cdot \phi(\Delta\eta, \Delta\phi, z)$  —  effective per-constituent kernel',
                     fontsize=14, y=1.01)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ── page 2: radial profiles ───────────────────────────────────────────
        fig, ax = plt.subplots(figsize=FIG_SIZE)
        N_BINS = 20
        r_edges = np.linspace(0, ETA_LIM * np.sqrt(2), N_BINS + 1)
        r_centres = 0.5 * (r_edges[:-1] + r_edges[1:])

        for z, col in zip(Z_VALS, Z_COLS):
            _, _, f2d = phi_kernel_grid(model, scaler, w_t, z, device)
            eta = np.linspace(-ETA_LIM, ETA_LIM, GRID_N)
            dEta_g, dPhi_g = np.meshgrid(eta, eta)
            r_flat = np.sqrt(dEta_g.ravel()**2 + dPhi_g.ravel()**2)
            f_flat = f2d.ravel()
            means, stds = [], []
            for lo, hi in zip(r_edges[:-1], r_edges[1:]):
                sel = (r_flat >= lo) & (r_flat < hi)
                if sel.sum() > 1:
                    means.append(f_flat[sel].mean())
                    stds.append(f_flat[sel].std())
                else:
                    means.append(np.nan)
                    stds.append(0.0)
            means = np.array(means)
            stds  = np.array(stds)
            ax.plot(r_centres, means, color=col, lw=2, label=rf'$z={z}$')
            ax.fill_between(r_centres, means - stds, means + stds,
                            color=col, alpha=0.2)

        ax.axhline(0, color='gray', lw=1, linestyle='--')
        apply_style(ax,
                    xlabel=r'$r = \sqrt{\Delta\eta^2 + \Delta\phi^2}$',
                    ylabel=r'$w \cdot \phi$',
                    title='Radial profile of effective kernel',
                    legend_loc='upper right')
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    print(f'Saved {out_pdf}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
