#!/usr/bin/env python3
"""
Inference script for NNjZ (summer26, Herwig Z): ROC curves comparing performance
on the Herwig H test split (training domain) and MG5_Pythia H (transfer test).

Inputs are identical to NNj plus the generator-boson 4-vector (pT, η, φ, M).
Both generator ROC curves plus single-variable baselines are overlaid on one plot.

Output:
  ../../figs/summer26/roc_nnjz.pdf
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

from dataset import load_features, BKG_FILE, SIG_FILE, JET_FEATURES, N_CONSTIT_FEAT
from model import JetNN
from style import apply_style, FIG_SIZE

PROCESS = "H"

BKG_FILE_MG5 = '../../friends/summer26/QCDHjj_mg5_pythia.h5'
SIG_FILE_MG5 = '../../friends/summer26/VBFH_mg5_pythia.h5'


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


# x_jet column order: |jetEta|, jetM, jetPt, then boson pT, η, φ, M
JET_IDX = {k: i for i, k in enumerate(JET_FEATURES)}

ROC_VARS = [
    (r'$p_T$ lead',    np.linspace(30,  500, 51), JET_IDX['jetPt']),
    (r'$|\eta|$ lead', np.linspace(0,   3,   31), JET_IDX['jetEta']),
    (r'$m$ lead',      np.linspace(0,   150, 51), JET_IDX['jetM']),
]

_LS = ['--', '-.', ':']


def scale_jcs(x_jcs, scaler):
    shape = x_jcs.shape
    return scaler.transform(x_jcs.reshape(-1, N_CONSTIT_FEAT)).reshape(shape).astype(np.float32)


def run_nn(model, scalers, x_jet_raw, x_jcs_raw, device):
    x_jet = scalers['jet'].transform(x_jet_raw).astype(np.float32)
    x_jcs = scale_jcs(x_jcs_raw, scalers['jcs'])
    mask  = (x_jcs_raw[..., 2] > 0)   # jcsPt > 0

    x_jet_t = torch.from_numpy(x_jet)
    x_jcs_t = torch.from_numpy(x_jcs)
    mask_t  = torch.from_numpy(mask)

    BATCH = 512
    logits_list = []
    with torch.no_grad():
        for start in range(0, len(x_jet), BATCH):
            end = start + BATCH
            out = model(
                x_jet_t[start:end].to(device),
                x_jcs_t[start:end].to(device),
                mask_t [start:end].to(device),
            )
            logits_list.append(out.cpu())
    logits = torch.cat(logits_list).numpy().ravel()
    return 1.0 / (1.0 + np.exp(-logits))


def main():
    # ── Load Herwig H (training domain) ──────────────────────────────────────
    print('Loading Herwig H ...')
    x_jet_bkg_hw, x_jcs_bkg_hw, _ = load_features(HERE / BKG_FILE)
    x_jet_sig_hw, x_jcs_sig_hw, _ = load_features(HERE / SIG_FILE)
    x_jet_all_hw = np.concatenate([x_jet_bkg_hw, x_jet_sig_hw], axis=0)
    x_jcs_all_hw = np.concatenate([x_jcs_bkg_hw, x_jcs_sig_hw], axis=0)
    y_all_hw = np.concatenate([np.zeros(len(x_jet_bkg_hw), dtype=np.float32),
                               np.ones (len(x_jet_sig_hw), dtype=np.float32)], axis=0)

    idx_test = np.load(HERE / f'split_indices_{PROCESS}.npz')['test_idx']
    x_jet_test_raw = x_jet_all_hw[idx_test]
    x_jcs_test_raw = x_jcs_all_hw[idx_test]
    y_test_hw      = y_all_hw[idx_test]
    print(f'Herwig test: {len(y_test_hw)} events  (sig={int(y_test_hw.sum())}, bkg={int((y_test_hw==0).sum())})')

    # ── Load MG5_Pythia H (transfer test, full dataset) ──────────────────────
    print('Loading MG5_Pythia H ...')
    x_jet_bkg_mg5, x_jcs_bkg_mg5, _ = load_features(HERE / BKG_FILE_MG5)
    x_jet_sig_mg5, x_jcs_sig_mg5, _ = load_features(HERE / SIG_FILE_MG5)
    x_jet_all_mg5 = np.concatenate([x_jet_bkg_mg5, x_jet_sig_mg5], axis=0)
    x_jcs_all_mg5 = np.concatenate([x_jcs_bkg_mg5, x_jcs_sig_mg5], axis=0)
    y_all_mg5 = np.concatenate([np.zeros(len(x_jet_bkg_mg5), dtype=np.float32),
                                np.ones (len(x_jet_sig_mg5), dtype=np.float32)], axis=0)
    print(f'MG5_Pythia: {len(y_all_mg5)} events  (sig={int(y_all_mg5.sum())}, bkg={int((y_all_mg5==0).sum())})')

    # ── Load model + scalers ──────────────────────────────────────────────────
    with open(HERE / f'scaler_{PROCESS}.pkl', 'rb') as fh:
        scalers = pickle.load(fh)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model  = JetNN().to(device)
    model.load_state_dict(torch.load(HERE / f'best_model_{PROCESS}.pt', map_location=device))
    model.eval()

    # ── NN scores ─────────────────────────────────────────────────────────────
    scores_hw  = run_nn(model, scalers, x_jet_test_raw, x_jcs_test_raw, device)
    scores_mg5 = run_nn(model, scalers, x_jet_all_mg5,  x_jcs_all_mg5,  device)

    fpr_hw,  tpr_hw,  auc_hw  = roc_from_scores(y_test_hw,  scores_hw)
    fpr_mg5, tpr_mg5, auc_mg5 = roc_from_scores(y_all_mg5,  scores_mg5)
    print(f'NNjZ Herwig    AUC = {auc_hw:.4f}')
    print(f'NNjZ MG5+Py    AUC = {auc_mg5:.4f}')

    # ── Single-variable ROCs (Herwig test set, jet scalars only) ─────────────
    x_jet_bkg_test = x_jet_test_raw[y_test_hw == 0]
    x_jet_sig_test = x_jet_test_raw[y_test_hw == 1]

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=FIG_SIZE)

    ax.plot(fpr_hw,  tpr_hw,  linewidth=3, color='black',
            label=f'NNjZ Herwig H  AUC={auc_hw:.3f}')
    ax.plot(fpr_mg5, tpr_mg5, linewidth=3, color='#1f77b4', linestyle='--',
            label=f'NNjZ MG5+Py H  AUC={auc_mg5:.3f}')

    for i, (label, bins, col_idx) in enumerate(ROC_VARS):
        fpr_v, tpr_v, auc_v = roc_from_histograms(
            x_jet_sig_test[:, col_idx], x_jet_bkg_test[:, col_idx], bins)
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
    out = FIG_DIR / f'roc_nnjz_{PROCESS}.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
