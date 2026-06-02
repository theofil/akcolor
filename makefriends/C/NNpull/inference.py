#!/usr/bin/env python3
"""
Inference script: loads the best trained model and produces ROC curves on the
held-out test partition, overlaid with individual-variable ROC curves.

Output:
  roc_nn.pdf  — ROC comparison plot saved next to this script
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

from dataset import load_features, BKG_FILE, SIG_FILE, JET_FEATURES, N_CONSTIT_FEAT
from model import JetNN
from style import apply_style, FIG_SIZE


# ── ROC helpers ───────────────────────────────────────────────────────────────

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
    sig_total = sig_counts.sum()
    bkg_total = bkg_counts.sum()
    if sig_total <= 0 or bkg_total <= 0:
        return None, None, None
    sig_cdf = np.concatenate([np.cumsum(sig_counts[::-1])[::-1], [0]]) / sig_total
    bkg_cdf = np.concatenate([np.cumsum(bkg_counts[::-1])[::-1], [0]]) / bkg_total
    tpr = sig_cdf[::-1]
    fpr = bkg_cdf[::-1]
    area = float(np.trapz(tpr, fpr))
    if area < 0.5:
        tpr = 1 - tpr[::-1]
        fpr = 1 - fpr[::-1]
        area = 1 - area
    return fpr, tpr, area


# ── Individual-variable ROC configurations ────────────────────────────────────

JET_IDX = {k: i for i, k in enumerate(JET_FEATURES)}
# JET_FEATURES order: jetEta(abs), jetM, jetNC, jetPVM, jetPt, jetSPVA

ROC_VARS = [
    (r'$|\theta_s|$ lead',  np.linspace(0,   np.pi, 51), lambda b, s: (np.abs(b[:, JET_IDX['jetSPVA']]), np.abs(s[:, JET_IDX['jetSPVA']]))),
    (r'$|\vec{t}\,|$ lead', np.linspace(0,   0.06,  51), lambda b, s: (b[:, JET_IDX['jetPVM']],          s[:, JET_IDX['jetPVM']])),
    (r'$p_T$ lead',         np.linspace(30,  500,   51), lambda b, s: (b[:, JET_IDX['jetPt']],           s[:, JET_IDX['jetPt']])),
    (r'$|\eta|$ lead',      np.linspace(0,   3,     31), lambda b, s: (b[:, JET_IDX['jetEta']],          s[:, JET_IDX['jetEta']])),
    (r'$m$ lead',           np.linspace(0,   150,   51), lambda b, s: (b[:, JET_IDX['jetM']],            s[:, JET_IDX['jetM']])),
    (r'$N_c$ lead',         np.linspace(0,   100,   26), lambda b, s: (b[:, JET_IDX['jetNC']],           s[:, JET_IDX['jetNC']])),
]

_LS = ['--', '-.', ':', (0,(3,1,1,1)), (0,(5,1)), (0,(1,1))]


def scale_jcs(x_jcs, scaler):
    shape = x_jcs.shape
    return scaler.transform(x_jcs.reshape(-1, N_CONSTIT_FEAT)).reshape(shape).astype(np.float32)


def main():
    # ── Load raw features ────────────────────────────────────────────────────
    print('Loading data ...')
    x_jet_bkg, x_jcs_bkg, _ = load_features(HERE / BKG_FILE)
    x_jet_sig, x_jcs_sig, _ = load_features(HERE / SIG_FILE)

    x_jet_all = np.concatenate([x_jet_bkg, x_jet_sig], axis=0)
    x_jcs_all = np.concatenate([x_jcs_bkg, x_jcs_sig], axis=0)
    y_all = np.concatenate([np.zeros(len(x_jet_bkg), dtype=np.float32),
                            np.ones (len(x_jet_sig), dtype=np.float32)], axis=0)

    indices  = np.load(HERE / 'split_indices.npz')
    idx_test = indices['test_idx']

    x_jet_test_raw = x_jet_all[idx_test]
    x_jcs_test_raw = x_jcs_all[idx_test]
    y_test         = y_all[idx_test]
    print(f'Test set: {len(y_test)} events  (sig={int(y_test.sum())}, bkg={int((y_test==0).sum())})')

    # ── Normalise ─────────────────────────────────────────────────────────────
    with open(HERE / 'scaler.pkl', 'rb') as fh:
        scalers = pickle.load(fh)

    x_jet_test = scalers['jet'].transform(x_jet_test_raw).astype(np.float32)
    x_jcs_test = scale_jcs(x_jcs_test_raw, scalers['jcs'])
    mask_test  = (x_jcs_test_raw[..., 3] > 0)   # jcsPt > 0, from raw data

    # ── Run NN ────────────────────────────────────────────────────────────────
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model  = JetNN().to(device)
    model.load_state_dict(torch.load(HERE / 'best_model.pt', map_location=device))
    model.eval()

    x_jet_t = torch.from_numpy(x_jet_test).to(device)
    x_jcs_t = torch.from_numpy(x_jcs_test).to(device)
    mask_t  = torch.from_numpy(mask_test).to(device)

    with torch.no_grad():
        logits = model(x_jet_t, x_jcs_t, mask_t).cpu().numpy().ravel()
    nn_scores = 1.0 / (1.0 + np.exp(-logits))

    fpr_nn, tpr_nn, auc_nn = roc_from_scores(y_test, nn_scores)
    print(f'NN test AUC = {auc_nn:.4f}')

    # ── Individual-variable ROC curves ────────────────────────────────────────
    x_jet_bkg_test = x_jet_test_raw[y_test == 0]
    x_jet_sig_test = x_jet_test_raw[y_test == 1]

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.plot(fpr_nn, tpr_nn, linewidth=3, color='black',
            label=f'NN  AUC={auc_nn:.3f}')

    for i, (label, bins, val_fn) in enumerate(ROC_VARS):
        bkg_vals, sig_vals = val_fn(x_jet_bkg_test, x_jet_sig_test)
        fpr_v, tpr_v, auc_v = roc_from_histograms(sig_vals, bkg_vals, bins)
        if fpr_v is None:
            continue
        ax.plot(fpr_v, tpr_v, linewidth=2, color='gray',
                linestyle=_LS[i % len(_LS)],
                label=f'{label}  AUC={auc_v:.3f}')

    ax.plot([0, 1], [0, 1], color='gray', linewidth=1, linestyle='--')
    apply_style(ax, xlabel='Background efficiency', ylabel='Signal efficiency',
                title='', xlim=(0, 1), ylim=(0, 1), legend_loc='lower right')
    plt.tight_layout()
    out = HERE / 'roc_nn.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
