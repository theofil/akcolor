#!/usr/bin/env python3
"""
Inference for DeepSetsNoJet: ROC curve on held-out test set.

Outputs (next to this script):
  scores.npz              — y_test, nn_scores  (read by compare_roc.py)
  roc_deepsets_nojet.pdf  — ROC vs individual jet variables
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
from model import DeepSetsNoJet
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


JET_IDX = {k: i for i, k in enumerate(JET_FEATURES)}
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
    print('Loading data ...')
    x_jet_bkg, x_jcs_bkg, _ = load_features(HERE / BKG_FILE)
    x_jet_sig, x_jcs_sig, _ = load_features(HERE / SIG_FILE)

    x_jet_all = np.concatenate([x_jet_bkg, x_jet_sig], axis=0)
    x_jcs_all = np.concatenate([x_jcs_bkg, x_jcs_sig], axis=0)
    y_all = np.concatenate([np.zeros(len(x_jet_bkg), dtype=np.float32),
                            np.ones (len(x_jet_sig), dtype=np.float32)], axis=0)

    idx_test = np.load(HERE / 'split_indices.npz')['test_idx']

    x_jet_test_raw = x_jet_all[idx_test]
    x_jcs_test_raw = x_jcs_all[idx_test]
    y_test         = y_all[idx_test]
    print(f'Test set: {len(y_test)} events  (sig={int(y_test.sum())}, bkg={int((y_test==0).sum())})')

    with open(HERE / 'scaler.pkl', 'rb') as fh:
        scalers = pickle.load(fh)

    x_jcs_test = scale_jcs(x_jcs_test_raw, scalers['jcs'])
    mask_test  = (x_jcs_test_raw[..., 3] > 0)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model  = DeepSetsNoJet().to(device)
    model.load_state_dict(torch.load(HERE / 'best_model.pt', map_location=device))
    model.eval()

    x_jcs_t = torch.from_numpy(x_jcs_test)
    mask_t  = torch.from_numpy(mask_test)

    BATCH = 512
    logits_list = []
    with torch.no_grad():
        for start in range(0, len(y_test), BATCH):
            end = start + BATCH
            out = model(x_jcs_t[start:end].to(device), mask_t[start:end].to(device))
            logits_list.append(out.cpu())
    logits    = torch.cat(logits_list).numpy().ravel()
    nn_scores = 1.0 / (1.0 + np.exp(-logits))

    np.savez(HERE / 'scores.npz', y_test=y_test, nn_scores=nn_scores)

    fpr_nn, tpr_nn, auc_nn = roc_from_scores(y_test, nn_scores)
    print(f'DeepSets (no jet) test AUC = {auc_nn:.4f}')

    x_jet_bkg_test = x_jet_test_raw[y_test == 0]
    x_jet_sig_test = x_jet_test_raw[y_test == 1]

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.plot(fpr_nn, tpr_nn, linewidth=3, color='black',
            label=f'DeepSets (no jet)  AUC={auc_nn:.3f}')

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
    out = HERE / 'roc_deepsets_nojet.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
