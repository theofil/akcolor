#!/usr/bin/env python3
"""
Symbolic regression on the effective per-constituent kernel of DeepSetsNoScale.

The effective contribution of constituent i is  w · phi(x_i)
where  w = E[d rho(h)/d h].  We fit a closed-form f(dEta, dPhi, z) ≈ w·phi(x_i)
using an 18-term polynomial basis (Ridge regression).  PySR is attempted first
if available; otherwise the polynomial fit is used.

Outputs:
  poly_model.npz     — coef (18,) for use in observable_roc.py
  pysr_summary.txt   — top terms, R², AUC of the resulting classical observable
  pysr_equations.csv — (only if PySR is available)
"""

import os
import sys
import pathlib

import numpy as np
import torch
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from interpret_utils import (
    load_model, load_test_set,
    compute_phi_aggregates, compute_mean_rho_jacobian,
    make_poly_basis, TERM_NAMES,
    roc_from_scores,
)

MAX_CONSTITUENTS = 200_000
SEED = 42


def collect_constituent_targets(model, w_t, x_raw, x_scaled, mask, device, batch=512):
    """
    For each valid constituent i: compute y_i = w · phi(scaled_x_i).
    Returns x_valid (M, 3) in physical units and y_valid (M,).
    """
    x_raw_list, y_list = [], []
    with torch.no_grad():
        for s in range(0, len(x_scaled), batch):
            xb  = torch.from_numpy(x_scaled[s:s + batch]).to(device)   # (B,80,3)
            mb  = mask[s:s + batch].astype(bool)                        # (B,80)
            B, N, F = xb.shape
            phi_out = model.phi(xb.view(B * N, F)).view(B, N, 64)       # (B,80,64)
            scores  = (phi_out * w_t).sum(-1).cpu().numpy()              # (B,80)
            valid   = mb.ravel()
            x_raw_list.append(x_raw[s:s + batch].reshape(-1, 3)[valid])
            y_list.append(scores.reshape(-1)[valid])
    return np.concatenate(x_raw_list), np.concatenate(y_list)


def compute_O_jet(coef, x_raw_test, mask_test):
    """Sum the polynomial observable over valid constituents for each jet."""
    N = len(x_raw_test)
    B_all  = make_poly_basis(x_raw_test.reshape(-1, 3))   # (N*80, 18)
    per_c  = (B_all @ coef).reshape(N, 80)                # (N, 80)
    per_c *= mask_test.astype(bool)                        # zero padding
    return per_c.sum(axis=1)                               # (N,)


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    model = load_model(device)

    print('Loading test set ...')
    x_raw, x_scaled, mask, y_test = load_test_set()
    print(f'  {len(y_test)} jets')

    print('Computing phi aggregates ...')
    h_test = compute_phi_aggregates(model, x_scaled, mask, device)

    print('Computing mean rho Jacobian w ...')
    w   = compute_mean_rho_jacobian(model, h_test, device)
    w_t = torch.from_numpy(w).to(device)

    print('Collecting per-constituent targets y_i = w·phi(x_i) ...')
    x_valid, y_valid = collect_constituent_targets(
        model, w_t, x_raw, x_scaled, mask, device)
    print(f'  {len(y_valid):,} valid constituents')

    # subsample if needed
    if len(y_valid) > MAX_CONSTITUENTS:
        rng = np.random.default_rng(SEED)
        idx = rng.choice(len(y_valid), MAX_CONSTITUENTS, replace=False)
        x_fit, y_fit = x_valid[idx], y_valid[idx]
        print(f'  subsampled to {MAX_CONSTITUENTS:,}')
    else:
        x_fit, y_fit = x_valid, y_valid

    # train/val split for reporting R²
    split = int(0.8 * len(x_fit))
    X_tr  = make_poly_basis(x_fit[:split])
    X_val = make_poly_basis(x_fit[split:])
    y_tr  = y_fit[:split]
    y_val = y_fit[split:]

    # ── PySR (optional) ───────────────────────────────────────────────────────
    pysr_ok = False
    try:
        from pysr import PySRRegressor
        print('\nTrying PySR ...')
        pysr_model = PySRRegressor(
            niterations=40,
            binary_operators=['+', '-', '*', '/'],
            unary_operators=['sqrt', 'square', 'abs'],
            maxsize=15,
            random_state=SEED,
            verbosity=0,
        )
        pysr_model.fit(x_fit, y_fit)
        pysr_model.equations_.to_csv(HERE / 'pysr_equations.csv', index=False)
        print('  PySR succeeded — equations saved to pysr_equations.csv')
        pysr_ok = True
    except Exception as exc:
        print(f'  PySR not available ({exc}); using polynomial Ridge fit.')

    # ── Polynomial Ridge (always run) ─────────────────────────────────────────
    print('\nFitting 18-term polynomial basis (Ridge) ...')
    ridge = Ridge(alpha=1e-4, fit_intercept=False)
    ridge.fit(X_tr, y_tr)
    coef = ridge.coef_

    r2_train = r2_score(y_tr,  ridge.predict(X_tr))
    r2_val   = r2_score(y_val, ridge.predict(X_val))
    print(f'  R² train={r2_train:.4f}  val={r2_val:.4f}')

    # sort by |coef| descending
    order = np.argsort(np.abs(coef))[::-1]
    print('\nTop 10 polynomial terms:')
    lines = ['Rank  Term        Coefficient']
    for rank, i in enumerate(order[:10], 1):
        line = f'  {rank:2d}   {TERM_NAMES[i]:12s}  {coef[i]:+.4f}'
        print(line)
        lines.append(line)

    # ── AUC of polynomial observable ─────────────────────────────────────────
    print('\nComputing AUC of polynomial observable ...')
    O_jet = compute_O_jet(coef, x_raw, mask)
    fpr, tpr, auc_poly = roc_from_scores(y_test, O_jet)
    print(f'  Polynomial observable AUC = {auc_poly:.4f}')
    print(f'  (DeepSets no-scale AUC   = 0.7810 from training)')

    # ── save ──────────────────────────────────────────────────────────────────
    np.savez(HERE / 'poly_model.npz', coef=coef)

    with open(HERE / 'pysr_summary.txt', 'w') as fh:
        fh.write('DeepSetsNoScale — polynomial kernel approximation\n')
        fh.write('=' * 50 + '\n\n')
        fh.write(f'Valid constituents used for fit: {len(y_valid):,}\n')
        fh.write(f'Subsample used:                  {len(x_fit):,}\n')
        fh.write(f'Basis:                           18-term polynomial\n')
        fh.write(f'R² (train):                      {r2_train:.4f}\n')
        fh.write(f'R² (val):                        {r2_val:.4f}\n')
        fh.write(f'AUC of O_jet = Σ f(dEta,dPhi,z):{auc_poly:.4f}\n')
        fh.write(f'DeepSets (no scale) AUC:         0.7810\n')
        fh.write(f'PySR available:                  {pysr_ok}\n\n')
        fh.write('\n'.join(lines) + '\n\n')
        fh.write('Full coefficient table:\n')
        fh.write(f'{"Term":<14}  {"Coefficient":>12}\n')
        for i in order:
            fh.write(f'{TERM_NAMES[i]:<14}  {coef[i]:+12.6f}\n')

    print(f'\nSaved poly_model.npz and pysr_summary.txt')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
