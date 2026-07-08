#!/usr/bin/env python3
"""
Inference script for NNkin (summer26, Herwig Z): ROC curves comparing performance
on the Herwig Z test split (training domain) and MG5_Pythia Z (transfer test).

Both generator ROC curves plus single-variable baselines are overlaid on one plot.
The model ROC curves are also saved to roc_data_{PROCESS}.npz for the summary
plot made by summer26/plots.py.

Output:
  ../../figs/summer26/roc_nnkin.pdf
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

HERE    = pathlib.Path(__file__).parent
FIG_DIR = HERE.parent.parent / 'figs' / 'summer26'
sys.path.insert(0, str(HERE))

from dataset import load_features, BKG_FILE, SIG_FILE, EVENT_FEATURES, JET_FEATURES
from model import KinNN
from style import apply_style, FIG_SIZE

PROCESS = "Z"

BKG_FILE_MG5 = '../../friends/summer26/QCDZjj_mg5_pythia.h5'
SIG_FILE_MG5 = '../../friends/summer26/VBFZ_mg5_pythia.h5'


def roc_from_scores(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    area = sk_auc(fpr, tpr)
    if area < 0.5:
        fpr, tpr, area = 1 - fpr[::-1], 1 - tpr[::-1], 1 - area
    return fpr, tpr, area


def decimate_roc(fpr, tpr, max_points=5000):
    """Uniform-stride subsampling keeping both endpoints (for the summary npz)."""
    if len(fpr) <= max_points:
        return fpr, tpr
    idx = np.linspace(0, len(fpr) - 1, max_points).round().astype(int)
    return fpr[idx], tpr[idx]


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


_EV_IDX = {k: i for i, k in enumerate(EVENT_FEATURES)}

ROC_VARS = [
    (r'$m_{jj}$',         np.linspace(0, 3000, 101),  _EV_IDX['mjj']),
    (r'$\Delta Y_{jj}$',  np.linspace(0, 8,    81),   _EV_IDX['dYjj']),
    (r'$\Delta\phi_{jj}$', np.linspace(0, np.pi, 51), _EV_IDX['dPhijj']),
    (r'$p_T^{jj}$',       np.linspace(0, 500, 51),    _EV_IDX['ptjj']),
]

_LS = ['--', '-.', ':', (0, (3, 1, 1, 1))]


def run_nn(model, scaler, x_raw, device):
    x_sc = scaler.transform(x_raw).astype(np.float32)
    with torch.no_grad():
        logits = model(torch.from_numpy(x_sc).to(device)).cpu().numpy().ravel()
    return 1.0 / (1.0 + np.exp(-logits))


def main():
    # ── Load Herwig Z (training domain) ──────────────────────────────────────
    print('Loading Herwig Z ...')
    x_bkg_hw, _ = load_features(HERE / BKG_FILE)
    x_sig_hw, _ = load_features(HERE / SIG_FILE)
    x_all_hw = np.concatenate([x_bkg_hw, x_sig_hw], axis=0)
    y_all_hw = np.concatenate([np.zeros(len(x_bkg_hw), dtype=np.float32),
                               np.ones (len(x_sig_hw), dtype=np.float32)], axis=0)

    idx_test = np.load(HERE / f'split_indices_{PROCESS}.npz')['test_idx']
    x_test_raw = x_all_hw[idx_test]
    y_test_hw  = y_all_hw[idx_test]
    print(f'Herwig test: {len(y_test_hw)} events  (sig={int(y_test_hw.sum())}, bkg={int((y_test_hw==0).sum())})')

    # ── Load MG5_Pythia Z (transfer test, full dataset) ──────────────────────
    print('Loading MG5_Pythia Z ...')
    x_bkg_mg5, _ = load_features(HERE / BKG_FILE_MG5)
    x_sig_mg5, _ = load_features(HERE / SIG_FILE_MG5)
    x_all_mg5 = np.concatenate([x_bkg_mg5, x_sig_mg5], axis=0)
    y_all_mg5 = np.concatenate([np.zeros(len(x_bkg_mg5), dtype=np.float32),
                                np.ones (len(x_sig_mg5), dtype=np.float32)], axis=0)
    print(f'MG5_Pythia: {len(y_all_mg5)} events  (sig={int(y_all_mg5.sum())}, bkg={int((y_all_mg5==0).sum())})')

    # ── Load model + scaler ───────────────────────────────────────────────────
    with open(HERE / f'scaler_{PROCESS}.pkl', 'rb') as fh:
        scaler = pickle.load(fh)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model  = KinNN().to(device)
    model.load_state_dict(torch.load(HERE / f'best_model_{PROCESS}.pt', map_location=device))
    model.eval()

    # ── NN scores ─────────────────────────────────────────────────────────────
    scores_hw  = run_nn(model, scaler, x_test_raw, device)
    scores_mg5 = run_nn(model, scaler, x_all_mg5,  device)

    fpr_hw,  tpr_hw,  auc_hw  = roc_from_scores(y_test_hw,  scores_hw)
    fpr_mg5, tpr_mg5, auc_mg5 = roc_from_scores(y_all_mg5,  scores_mg5)
    print(f'NNkin Herwig    AUC = {auc_hw:.4f}')
    print(f'NNkin MG5+Py    AUC = {auc_mg5:.4f}')

    # ── Save ROC data for the all-NN summary plot (summer26/plots.py) ────────
    fpr_hw_s,  tpr_hw_s  = decimate_roc(fpr_hw,  tpr_hw)
    fpr_mg5_s, tpr_mg5_s = decimate_roc(fpr_mg5, tpr_mg5)
    np.savez(HERE / f'roc_data_{PROCESS}.npz',
             fpr_hw=fpr_hw_s, tpr_hw=tpr_hw_s, auc_hw=auc_hw,
             fpr_mg5=fpr_mg5_s, tpr_mg5=tpr_mg5_s, auc_mg5=auc_mg5)
    print(f'Saved {HERE / f"roc_data_{PROCESS}.npz"}')


    # ── Single-variable ROCs (Herwig test set) ────────────────────────────────
    x_bkg_test = x_test_raw[y_test_hw == 0]
    x_sig_test = x_test_raw[y_test_hw == 1]

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=FIG_SIZE)

    ax.plot(fpr_hw,  tpr_hw,  linewidth=3, color='black',
            label=f'NNkin Herwig Z  AUC={auc_hw:.3f}')
    ax.plot(fpr_mg5, tpr_mg5, linewidth=3, color='#1f77b4', linestyle='--',
            label=f'NNkin MG5+Py Z  AUC={auc_mg5:.3f}')

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
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / f'roc_nnkin_{PROCESS}.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
