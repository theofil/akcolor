#!/usr/bin/env python3
"""
Compare ROC curves of NNkin, NNpull, and NNtrans on their respective test sets,
overlaid with mjj and dYjj single-variable ROC curves.

Requires NNkin, NNpull, and NNtrans to have been trained first.

Output:
  compare_roc.pdf  — saved next to this script
"""

import os
import sys
import pickle
import pathlib

import numpy as np
import torch
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc as sk_auc

HERE    = pathlib.Path(__file__).parent
NNPULL  = HERE.parent / 'NNpull'
NNTRANS = HERE.parent / 'NNtrans'

# ── import NNkin ──────────────────────────────────────────────────────────────
sys.path.insert(0, str(HERE))
import dataset as kin_ds
import model   as kin_mod

# ── import NNpull (use importlib to avoid name collision) ─────────────────────
import importlib.util

def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

pull_ds  = _load_module('pull_dataset',  NNPULL  / 'dataset.py')
trans_ds = _load_module('trans_dataset', NNTRANS / 'dataset.py')

# NNpull/model.py and NNtrans/model.py do `from dataset import ...`; temporarily
# expose the right module under that name so imports resolve correctly.
_prev_dataset = sys.modules.get('dataset')

sys.modules['dataset'] = pull_ds
pull_mod = _load_module('pull_model', NNPULL / 'model.py')

sys.modules['dataset'] = trans_ds
trans_mod = _load_module('trans_model', NNTRANS / 'model.py')

if _prev_dataset is None:
    del sys.modules['dataset']
else:
    sys.modules['dataset'] = _prev_dataset

BKG_FILE = HERE / kin_ds.BKG_FILE
SIG_FILE = HERE / kin_ds.SIG_FILE


# ── ROC helpers ───────────────────────────────────────────────────────────────

def roc_from_scores(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    area = sk_auc(fpr, tpr)
    if area < 0.5:
        fpr, tpr, area = 1 - fpr[::-1], 1 - tpr[::-1], 1 - area
    return fpr, tpr, area


def roc_from_histograms(sig_vals, bkg_vals, bins):
    sig_c, _ = np.histogram(sig_vals, bins=bins)
    bkg_c, _ = np.histogram(bkg_vals, bins=bins)
    sig_c = sig_c.astype(float)
    bkg_c = bkg_c.astype(float)
    if sig_c.sum() <= 0 or bkg_c.sum() <= 0:
        return None, None, None
    sig_cdf = np.concatenate([np.cumsum(sig_c[::-1])[::-1], [0]]) / sig_c.sum()
    bkg_cdf = np.concatenate([np.cumsum(bkg_c[::-1])[::-1], [0]]) / bkg_c.sum()
    tpr  = sig_cdf[::-1]
    fpr  = bkg_cdf[::-1]
    area = float(np.trapz(tpr, fpr))
    if area < 0.5:
        tpr  = 1 - tpr[::-1]
        fpr  = 1 - fpr[::-1]
        area = 1 - area
    return fpr, tpr, area


def scale_jcs(x_jcs, scaler):
    shape = x_jcs.shape
    return scaler.transform(x_jcs.reshape(-1, pull_ds.N_CONSTIT_FEAT)).reshape(shape).astype(np.float32)


def nn_scores_kin(x_all_raw, y_all, idx_test, device):
    """Run NNkin on its test partition."""
    x_test_raw = x_all_raw[idx_test]
    y_test     = y_all[idx_test]

    with open(HERE / 'scaler.pkl', 'rb') as fh:
        scaler = pickle.load(fh)
    x_test = scaler.transform(x_test_raw).astype(np.float32)

    model = kin_mod.KinNN().to(device)
    model.load_state_dict(torch.load(HERE / 'best_model.pt', map_location=device))
    model.eval()
    with torch.no_grad():
        logits = model(torch.from_numpy(x_test).to(device)).cpu().numpy().ravel()
    scores = 1.0 / (1.0 + np.exp(-logits))
    return scores, y_test, x_test_raw


def nn_scores_pull(x_jet_all, x_jcs_all, y_all, idx_test, device):
    """Run NNpull on its test partition."""
    x_jet_raw = x_jet_all[idx_test]
    x_jcs_raw = x_jcs_all[idx_test]
    y_test    = y_all[idx_test]

    with open(NNPULL / 'scaler.pkl', 'rb') as fh:
        scalers = pickle.load(fh)
    x_jet = scalers['jet'].transform(x_jet_raw).astype(np.float32)
    x_jcs = scale_jcs(x_jcs_raw, scalers['jcs'])
    mask  = (x_jcs_raw[..., 3] > 0)   # jcsPt > 0

    model = pull_mod.JetNN().to(device)
    model.load_state_dict(torch.load(NNPULL / 'best_model.pt', map_location=device))
    model.eval()
    with torch.no_grad():
        logits = model(
            torch.from_numpy(x_jet).to(device),
            torch.from_numpy(x_jcs).to(device),
            torch.from_numpy(mask).to(device),
        ).cpu().numpy().ravel()
    scores = 1.0 / (1.0 + np.exp(-logits))
    return scores, y_test


def nn_scores_trans(x_jet_all, x_jcs_all, y_all, idx_test, device):
    """Run NNtrans on its test partition (batched to handle long sequences)."""
    x_jet_raw = x_jet_all[idx_test]
    x_jcs_raw = x_jcs_all[idx_test]
    y_test    = y_all[idx_test]

    with open(NNTRANS / 'scaler.pkl', 'rb') as fh:
        scalers = pickle.load(fh)
    x_jet = scalers['jet'].transform(x_jet_raw).astype(np.float32)
    x_jcs = scale_jcs(x_jcs_raw, scalers['jcs'])
    mask  = (x_jcs_raw[..., 3] > 0)   # jcsPt > 0

    model = trans_mod.TransJetNN().to(device)
    model.load_state_dict(torch.load(NNTRANS / 'best_model.pt', map_location=device))
    model.eval()

    x_jet_t = torch.from_numpy(x_jet)
    x_jcs_t = torch.from_numpy(x_jcs)
    mask_t  = torch.from_numpy(mask)

    BATCH = 512
    logits_list = []
    with torch.no_grad():
        for start in range(0, len(y_test), BATCH):
            end = start + BATCH
            logits_list.append(model(
                x_jet_t[start:end].to(device),
                x_jcs_t[start:end].to(device),
                mask_t [start:end].to(device),
            ).cpu())
    logits = torch.cat(logits_list).numpy().ravel()
    scores = 1.0 / (1.0 + np.exp(-logits))
    return scores, y_test


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    # ── Load raw kinematic features ───────────────────────────────────────────
    print('Loading kinematic features ...')
    x_kin_bkg, _ = kin_ds.load_features(BKG_FILE)
    x_kin_sig, _ = kin_ds.load_features(SIG_FILE)
    x_kin_all = np.concatenate([x_kin_bkg, x_kin_sig], axis=0)   # (N, 12)

    # ── Load NNpull jet+constituent features ──────────────────────────────────
    print('Loading pull features ...')
    x_jet_bkg, x_jcs_bkg, _ = pull_ds.load_features(BKG_FILE)
    x_jet_sig, x_jcs_sig, _ = pull_ds.load_features(SIG_FILE)
    x_jet_all = np.concatenate([x_jet_bkg, x_jet_sig], axis=0)
    x_jcs_all = np.concatenate([x_jcs_bkg, x_jcs_sig], axis=0)

    y_all = np.concatenate([np.zeros(len(x_kin_bkg), dtype=np.float32),
                            np.ones (len(x_kin_sig), dtype=np.float32)], axis=0)

    # test partitions (same SEED → same split)
    idx_kin   = np.load(HERE    / 'split_indices.npz')['test_idx']
    idx_pull  = np.load(NNPULL  / 'split_indices.npz')['test_idx']
    idx_trans = np.load(NNTRANS / 'split_indices.npz')['test_idx']

    # ── NNkin ─────────────────────────────────────────────────────────────────
    print('Running NNkin ...')
    sc_kin, y_test_kin, x_test_raw_kin = nn_scores_kin(x_kin_all, y_all, idx_kin, device)
    fpr_kin, tpr_kin, auc_kin = roc_from_scores(y_test_kin, sc_kin)
    print(f'NNkin  AUC = {auc_kin:.4f}')

    # ── NNpull ────────────────────────────────────────────────────────────────
    print('Running NNpull ...')
    sc_pull, y_test_pull = nn_scores_pull(x_jet_all, x_jcs_all, y_all, idx_pull, device)
    fpr_pull, tpr_pull, auc_pull = roc_from_scores(y_test_pull, sc_pull)
    print(f'NNpull AUC = {auc_pull:.4f}')

    # ── NNtrans ───────────────────────────────────────────────────────────────
    print('Running NNtrans ...')
    sc_trans, y_test_trans = nn_scores_trans(x_jet_all, x_jcs_all, y_all, idx_trans, device)
    fpr_trans, tpr_trans, auc_trans = roc_from_scores(y_test_trans, sc_trans)
    print(f'NNtrans AUC = {auc_trans:.4f}')

    # ── Single-variable ROCs (mjj, dYjj) from NNkin test set ─────────────────
    ev_idx = {k: i for i, k in enumerate(kin_ds.EVENT_FEATURES)}
    x_bkg_test = x_test_raw_kin[y_test_kin == 0]
    x_sig_test = x_test_raw_kin[y_test_kin == 1]

    fpr_mjj,  tpr_mjj,  auc_mjj  = roc_from_histograms(
        x_sig_test[:, ev_idx['mjj']],  x_bkg_test[:, ev_idx['mjj']],
        np.linspace(0, 3000, 101))
    fpr_dyj,  tpr_dyj,  auc_dyj  = roc_from_histograms(
        x_sig_test[:, ev_idx['dYjj']], x_bkg_test[:, ev_idx['dYjj']],
        np.linspace(0, 8, 81))

    print(f'mjj    AUC = {auc_mjj:.4f}')
    print(f'dYjj   AUC = {auc_dyj:.4f}')

    # ── Plot ──────────────────────────────────────────────────────────────────
    from style import apply_style, FIG_SIZE
    fig, ax = plt.subplots(figsize=FIG_SIZE)

    ax.plot(fpr_trans, tpr_trans, linewidth=3, color='#2ca02c',
            label=f'NNtrans AUC={auc_trans:.3f}')
    ax.plot(fpr_pull,  tpr_pull,  linewidth=3, color='#d62728',
            label=f'NNpull  AUC={auc_pull:.3f}')
    ax.plot(fpr_kin,   tpr_kin,   linewidth=3, color='#1f77b4',
            label=f'NNkin   AUC={auc_kin:.3f}')
    ax.plot(fpr_mjj,  tpr_mjj,  linewidth=2, color='gray', linestyle='--',
            label=rf'$m_{{jj}}$   AUC={auc_mjj:.3f}')
    ax.plot(fpr_dyj,  tpr_dyj,  linewidth=2, color='gray', linestyle='-.',
            label=rf'$\Delta Y_{{jj}}$  AUC={auc_dyj:.3f}')
    ax.plot([0, 1], [0, 1], color='lightgray', linewidth=1, linestyle='--')

    apply_style(ax, xlabel='Background efficiency', ylabel='Signal efficiency',
                title='', xlim=(0, 1), ylim=(0, 1), legend_loc='lower right')
    plt.tight_layout()
    out = HERE / 'compare_roc.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
