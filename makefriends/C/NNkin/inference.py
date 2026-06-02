#!/usr/bin/env python3
"""
Inference script for NNkin: ROC curves on the held-out test set,
overlaid with mjj and dYjj single-variable ROC curves.

Output:
  roc_nnkin.pdf  — saved next to this script
"""

import os
import sys
import pickle
import pathlib

import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc as sk_auc

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))

from dataset import load_features, BKG_FILE, SIG_FILE, EVENT_FEATURES, JET_FEATURES
from model import KinNN
from style import apply_style, FIG_SIZE


def roc_from_scores(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    area = sk_auc(fpr, tpr)
    if area < 0.5:
        fpr, tpr, area = 1 - fpr[::-1], 1 - tpr[::-1], 1 - area
    return fpr, tpr, area


def roc_from_histograms(sig_vals, bkg_vals, bins):
    sig_counts, _ = np.histogram(sig_vals, bins=bins)
    bkg_counts, _ = np.histogram(bkg_vals, bins=bins)
    sig_counts = sig_counts.astype(float)
    bkg_counts = bkg_counts.astype(float)
    if sig_counts.sum() <= 0 or bkg_counts.sum() <= 0:
        return None, None, None
    sig_cdf = np.concatenate([np.cumsum(sig_counts[::-1])[::-1], [0]]) / sig_counts.sum()
    bkg_cdf = np.concatenate([np.cumsum(bkg_counts[::-1])[::-1], [0]]) / bkg_counts.sum()
    tpr = sig_cdf[::-1]
    fpr = bkg_cdf[::-1]
    area = float(np.trapz(tpr, fpr))
    if area < 0.5:
        tpr = 1 - tpr[::-1]
        fpr = 1 - fpr[::-1]
        area = 1 - area
    return fpr, tpr, area


# event-level feature indices
_EV_IDX = {k: i for i, k in enumerate(EVENT_FEATURES)}
# dPhijj=0, dYjj=1, mjj=2, ptjj=3

ROC_VARS = [
    (r'$m_{jj}$',       np.linspace(0, 3000, 101),  _EV_IDX['mjj']),
    (r'$\Delta Y_{jj}$', np.linspace(0, 8,    81),   _EV_IDX['dYjj']),
    (r'$\Delta\phi_{jj}$', np.linspace(0, np.pi, 51), _EV_IDX['dPhijj']),
    (r'$p_T^{jj}$',     np.linspace(0, 500, 51),    _EV_IDX['ptjj']),
]

_LS = ['--', '-.', ':', (0,(3,1,1,1))]


def main():
    print('Loading data ...')
    x_bkg, _ = load_features(HERE / BKG_FILE)
    x_sig, _ = load_features(HERE / SIG_FILE)

    x_all = np.concatenate([x_bkg, x_sig], axis=0)
    y_all = np.concatenate([np.zeros(len(x_bkg), dtype=np.float32),
                            np.ones (len(x_sig), dtype=np.float32)], axis=0)

    idx_test = np.load(HERE / 'split_indices.npz')['test_idx']
    x_test_raw = x_all[idx_test]
    y_test     = y_all[idx_test]
    print(f'Test set: {len(y_test)} events  (sig={int(y_test.sum())}, bkg={int((y_test==0).sum())})')

    with open(HERE / 'scaler.pkl', 'rb') as fh:
        scaler = pickle.load(fh)
    x_test = scaler.transform(x_test_raw).astype(np.float32)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model  = KinNN().to(device)
    model.load_state_dict(torch.load(HERE / 'best_model.pt', map_location=device))
    model.eval()

    with torch.no_grad():
        logits = model(torch.from_numpy(x_test).to(device)).cpu().numpy().ravel()
    nn_scores = 1.0 / (1.0 + np.exp(-logits))

    fpr_nn, tpr_nn, auc_nn = roc_from_scores(y_test, nn_scores)
    print(f'NNkin test AUC = {auc_nn:.4f}')

    # single-variable ROCs use raw (unscaled) event features
    x_bkg_test = x_test_raw[y_test == 0]
    x_sig_test = x_test_raw[y_test == 1]

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.plot(fpr_nn, tpr_nn, linewidth=3, color='black',
            label=f'NNkin  AUC={auc_nn:.3f}')

    for i, (label, bins, col_idx) in enumerate(ROC_VARS):
        fpr_v, tpr_v, auc_v = roc_from_histograms(
            x_sig_test[:, col_idx], x_bkg_test[:, col_idx], bins)
        if fpr_v is None:
            continue
        ax.plot(fpr_v, tpr_v, linewidth=2, color='gray',
                linestyle=_LS[i % len(_LS)],
                label=f'{label}  AUC={auc_v:.3f}')

    ax.plot([0, 1], [0, 1], color='gray', linewidth=1, linestyle='--')
    apply_style(ax, xlabel='Background efficiency', ylabel='Signal efficiency',
                title='', xlim=(0, 1), ylim=(0, 1), legend_loc='lower right')
    plt.tight_layout()
    out = HERE / 'roc_nnkin.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
