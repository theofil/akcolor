#!/usr/bin/env python3
"""
Compare three DeepSets ablations on the same held-out test set:

  1. DeepSets (no jet)   — 5 constituent features (jcsDEta, jcsDPhi, jcsM, jcsPt, jcsW)
  2. DeepSets (no fun)   — 4 constituent features (jcsDEta, jcsDPhi, jcsM, jcsPt)
  3. DeepSets (no scale) — 3 constituent features (jcsDEta, jcsDPhi, jcsPt/jetPt)

Output: roc_comparison_DeepSetNoScale.pdf
"""

import os
import sys
import pickle
import pathlib
import importlib.util

import numpy as np
import h5py
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc as sk_auc

HERE = pathlib.Path(__file__).parent   # testAUC/
ROOT = HERE.parent                     # makefriends/C/

BKG_H5 = ROOT / 'friends' / 'QCDHtoInv.h5'
SIG_H5 = ROOT / 'friends' / 'VBFHtoInv.h5'

LABEL_SIZE  = 20
LEGEND_SIZE = 14
TICK_SIZE   = 16
FIG_SIZE    = (6, 6)


# ── helpers ───────────────────────────────────────────────────────────────────

def load_module(name, path):
    spec   = importlib.util.spec_from_file_location(name, path)
    mod    = importlib.util.module_from_spec(spec)
    parent = str(pathlib.Path(path).parent)
    sys.path.insert(0, parent)
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.path.remove(parent)
    return mod


def load_data():
    """Load all 5 constituent features + jetPt (for fraction). Returns raw arrays."""
    CONSTIT_FEATURES = ['jcsDEta', 'jcsDPhi', 'jcsM', 'jcsPt', 'jcsW']

    def _read(h5path):
        with h5py.File(h5path, 'r') as f:
            jetPt = f['jetPt'][:, 0].astype(np.float32)           # (N,)
            jcs   = {k: f[k][:, 0, :].astype(np.float32) for k in CONSTIT_FEATURES}
        x_jcs = np.stack([jcs[k] for k in CONSTIT_FEATURES], axis=2)  # (N, 80, 5)
        return x_jcs, jetPt

    x_jcs_bkg, jetPt_bkg = _read(BKG_H5)
    x_jcs_sig, jetPt_sig = _read(SIG_H5)

    x_jcs  = np.concatenate([x_jcs_bkg,  x_jcs_sig],  axis=0)
    jetPt  = np.concatenate([jetPt_bkg,   jetPt_sig],   axis=0)
    y      = np.concatenate([np.zeros(len(x_jcs_bkg), dtype=np.float32),
                             np.ones (len(x_jcs_sig), dtype=np.float32)], axis=0)
    return x_jcs, jetPt, y


def scale_jcs(x_jcs, scaler):
    n_feat = x_jcs.shape[-1]
    shape  = x_jcs.shape
    return scaler.transform(x_jcs.reshape(-1, n_feat)).reshape(shape).astype(np.float32)


def roc_from_scores(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    area = sk_auc(fpr, tpr)
    if area < 0.5:
        fpr, tpr, area = 1 - fpr[::-1], 1 - tpr[::-1], 1 - area
    return fpr, tpr, area


def run_inference_no_jet(model, x_jcs_test, mask_test, device, batch=512):
    x_jcs_t = torch.from_numpy(x_jcs_test)
    mask_t  = torch.from_numpy(mask_test)
    logits_list = []
    model.eval()
    with torch.no_grad():
        for s in range(0, len(x_jcs_test), batch):
            e = s + batch
            out = model(x_jcs_t[s:e].to(device), mask_t[s:e].to(device))
            logits_list.append(out.cpu())
    logits = torch.cat(logits_list).numpy().ravel()
    return 1.0 / (1.0 + np.exp(-logits))


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    print('Loading data ...')
    x_jcs_all, jetPt_all, y_all = load_data()

    idx_test       = np.load(HERE / 'deepsets' / 'split_indices.npz')['test_idx']
    x_jcs_test_raw = x_jcs_all[idx_test]   # (M, 80, 5)  all 5 features
    jetPt_test     = jetPt_all[idx_test]   # (M,)
    y_test         = y_all[idx_test]
    # mask from raw jcsPt (index 3 in 5-feature array)
    mask_test      = (x_jcs_test_raw[..., 3] > 0)
    print(f'Test set: {len(y_test)} events  (sig={int(y_test.sum())}, bkg={int((y_test==0).sum())})')

    results = []

    # ── 1. DeepSets (no jet) — 5 features ────────────────────────────────────
    print('\n--- DeepSets (no jet) ---')
    ds_nj_dir = HERE / 'deepsets'
    ds_nj_mod = load_module('deepsets_nojet_model', ds_nj_dir / 'model.py')
    with open(ds_nj_dir / 'scaler.pkl', 'rb') as fh:
        sc = pickle.load(fh)
    x_jcs_s = scale_jcs(x_jcs_test_raw, sc['jcs'])          # 5-feature scaler
    model = ds_nj_mod.DeepSetsNoJet().to(device)
    model.load_state_dict(torch.load(ds_nj_dir / 'best_model.pt', map_location=device))
    scores = run_inference_no_jet(model, x_jcs_s, mask_test, device)
    fpr, tpr, auc = roc_from_scores(y_test, scores)
    print(f'  AUC = {auc:.4f}')
    results.append(('DeepSets (no jet)', fpr, tpr, auc, '#d62728', '--', 2))

    # ── 2. DeepSets (no fun) — 4 features (no jcsW) ─────────────────────────
    print('--- DeepSets (no fun) ---')
    ds_nf_dir = HERE / 'deepsets_nofun'
    ds_nf_mod = load_module('deepsets_nofun_model', ds_nf_dir / 'model.py')
    with open(ds_nf_dir / 'scaler.pkl', 'rb') as fh:
        sc = pickle.load(fh)
    x_jcs_4feat = x_jcs_test_raw[..., :4]                   # drop jcsW (column 4)
    x_jcs_s = scale_jcs(x_jcs_4feat, sc['jcs'])             # 4-feature scaler
    model = ds_nf_mod.DeepSetsNoFun().to(device)
    model.load_state_dict(torch.load(ds_nf_dir / 'best_model.pt', map_location=device))
    scores = run_inference_no_jet(model, x_jcs_s, mask_test, device)
    fpr, tpr, auc = roc_from_scores(y_test, scores)
    print(f'  AUC = {auc:.4f}')
    results.append(('DeepSets (no fun)', fpr, tpr, auc, '#9467bd', '--', 2))

    # ── 3. DeepSets (no scale) — 3 features (jcsDEta, jcsDPhi, jcsPt/jetPt) ─
    print('--- DeepSets (no scale) ---')
    ds_ns_dir = HERE / 'deepsets_noscale'
    ds_ns_mod = load_module('deepsets_noscale_model', ds_ns_dir / 'model.py')
    with open(ds_ns_dir / 'scaler.pkl', 'rb') as fh:
        sc = pickle.load(fh)
    # reconstruct 3-feature array: [jcsDEta, jcsDPhi, jcsPt/jetPt]
    jcsPt_frac = x_jcs_test_raw[..., 3] / jetPt_test[:, np.newaxis]  # (M, 80)
    x_jcs_3feat = np.stack([x_jcs_test_raw[..., 0],
                             x_jcs_test_raw[..., 1],
                             jcsPt_frac], axis=2).astype(np.float32)  # (M, 80, 3)
    x_jcs_s = scale_jcs(x_jcs_3feat, sc['jcs'])             # 3-feature scaler
    model = ds_ns_mod.DeepSetsNoScale().to(device)
    model.load_state_dict(torch.load(ds_ns_dir / 'best_model.pt', map_location=device))
    scores = run_inference_no_jet(model, x_jcs_s, mask_test, device)
    fpr, tpr, auc = roc_from_scores(y_test, scores)
    print(f'  AUC = {auc:.4f}')
    results.append(('DeepSets (no scale)', fpr, tpr, auc, '#2ca02c', '--', 2))

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for label, fpr, tpr, auc, color, ls, lw in results:
        ax.plot(fpr, tpr, linewidth=lw, color=color, linestyle=ls,
                label=f'{label}  AUC={auc:.3f}')
    ax.plot([0, 1], [0, 1], color='gray', linewidth=1, linestyle='--')

    ax.set_xlabel('Background efficiency', fontsize=LABEL_SIZE)
    ax.set_ylabel('Signal efficiency',     fontsize=LABEL_SIZE)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', labelsize=TICK_SIZE)
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.legend(loc='lower right', frameon=False, fontsize=LEGEND_SIZE)

    plt.tight_layout()
    out = HERE / 'roc_comparison_DeepSetNoScale.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSaved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
