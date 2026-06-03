#!/usr/bin/env python3
"""
ROC comparison: polynomial classical observable vs. DeepSets (no scale) vs.
pull-vector variables jetPVM and |jetSPVA|.

Requires poly_model.npz produced by phi_pysr.py.

Output: observable_roc.pdf
"""

import os
import sys
import pathlib

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from interpret_utils import (
    load_test_set, load_jet_features_test,
    make_poly_basis, TERM_NAMES,
    roc_from_scores,
    DS_DIR,
    apply_style, FIG_SIZE, LABEL_SIZE, LEGEND_SIZE, TICK_SIZE,
)


def compute_O_jet(coef, x_raw_test, mask_test):
    """Σ_i f(dEta_i, dPhi_i, z_i) over valid constituents per jet."""
    N     = len(x_raw_test)
    B_all = make_poly_basis(x_raw_test.reshape(-1, 3))   # (N*80, 18)
    per_c = (B_all @ coef).reshape(N, 80)
    per_c *= mask_test.astype(bool)
    return per_c.sum(axis=1)                              # (N,)


def top_terms_label(coef, n=4):
    """Build a short label showing the dominant polynomial terms."""
    order = np.argsort(np.abs(coef))[::-1]
    parts = []
    for i in order[:n]:
        sign = '+' if coef[i] >= 0 else '-'
        parts.append(f'{sign}{abs(coef[i]):.2f}·{TERM_NAMES[i]}')
    return r'$O_{\rm jet}$: ' + '  '.join(parts[:2]) + ' ...'


def main():
    poly_path = HERE / 'poly_model.npz'
    if not poly_path.exists():
        raise FileNotFoundError(
            f'{poly_path} not found — run phi_pysr.py first.')

    coef = np.load(poly_path)['coef']

    print('Loading test set ...')
    x_raw, _, mask, y_test = load_test_set()
    print(f'  {len(y_test)} jets')

    print('Computing polynomial observable O_jet ...')
    O_jet = compute_O_jet(coef, x_raw, mask)
    fpr_poly, tpr_poly, auc_poly = roc_from_scores(y_test, O_jet)
    print(f'  O_jet AUC = {auc_poly:.4f}')

    print('Loading pre-computed DeepSets (no scale) scores ...')
    data_nn    = np.load(DS_DIR / 'scores.npz')
    nn_scores  = data_nn['nn_scores']
    y_nn       = data_nn['y_test']
    fpr_nn, tpr_nn, auc_nn = roc_from_scores(y_nn, nn_scores)
    print(f'  DeepSets AUC = {auc_nn:.4f}')

    print('Loading jet pull variables ...')
    jetPVM, jetSPVA = load_jet_features_test()
    fpr_pvm,  tpr_pvm,  auc_pvm  = roc_from_scores(y_test, jetPVM)
    fpr_spva, tpr_spva, auc_spva = roc_from_scores(y_test, np.abs(jetSPVA))
    print(f'  jetPVM  AUC = {auc_pvm:.4f}')
    print(f'  |jetSPVA| AUC = {auc_spva:.4f}')

    # ── plot ──────────────────────────────────────────────────────────────────
    obs_label = top_terms_label(coef)

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.plot(fpr_nn,   tpr_nn,   color='black',   lw=3,
            label=f'DeepSets (no scale)  AUC={auc_nn:.3f}')
    ax.plot(fpr_poly, tpr_poly, color='#1f77b4', lw=2,
            label=f'Polynomial {obs_label}  AUC={auc_poly:.3f}')
    ax.plot(fpr_pvm,  tpr_pvm,  color='#d62728', lw=2, linestyle='--',
            label=rf'$|\vec{{t}}\,|$  AUC={auc_pvm:.3f}')
    ax.plot(fpr_spva, tpr_spva, color='#2ca02c', lw=2, linestyle='--',
            label=rf'$|\theta_s|$  AUC={auc_spva:.3f}')
    ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')

    apply_style(ax,
                xlabel='Background efficiency',
                ylabel='Signal efficiency',
                title='Classical observable vs network',
                xlim=(0, 1), ylim=(0, 1),
                legend_loc='lower right')
    plt.tight_layout()
    out = HERE / 'observable_roc.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
