#!/usr/bin/env python3
"""
Compare ROC curves for all four models on the same test set:
  1. NNpull    — DeepSets   + jet features   (existing, ../NNpull/)
  2. NNtrans   — Transformer + jet features  (existing, ../NNtrans/)
  3. DeepSets  — DeepSets   NO jet features  (testAUC/deepsets/)
  4. Transformer NO jet features             (testAUC/transformer/)

Requires: all four best_model.pt + scaler.pkl + split_indices.npz to exist.

Output: roc_comparison.pdf  (saved next to this script)
"""

import os
import sys
import pickle
import pathlib
import importlib.util

import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc as sk_auc

HERE    = pathlib.Path(__file__).parent   # testAUC/
ROOT    = HERE.parent                     # makefriends/C/

BKG_H5  = ROOT / 'friends' / 'QCDHtoInv.h5'
SIG_H5  = ROOT / 'friends' / 'VBFHtoInv.h5'

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
    import h5py
    JET_FEATURES     = ['jetEta', 'jetM', 'jetNC', 'jetPVM', 'jetPt', 'jetSPVA']
    CONSTIT_FEATURES = ['jcsDEta', 'jcsDPhi', 'jcsM', 'jcsPt', 'jcsW']

    def _read(h5path):
        with h5py.File(h5path, 'r') as f:
            jets = {k: f[k][:, 0].astype(np.float32)    for k in JET_FEATURES}
            jcs  = {k: f[k][:, 0, :].astype(np.float32) for k in CONSTIT_FEATURES}
        jets['jetEta'] = np.abs(jets['jetEta'])
        x_jet = np.stack([jets[k] for k in JET_FEATURES],     axis=1)
        x_jcs = np.stack([jcs[k]  for k in CONSTIT_FEATURES], axis=2)
        return x_jet, x_jcs

    x_jet_bkg, x_jcs_bkg = _read(BKG_H5)
    x_jet_sig, x_jcs_sig = _read(SIG_H5)

    x_jet = np.concatenate([x_jet_bkg, x_jet_sig], axis=0)
    x_jcs = np.concatenate([x_jcs_bkg, x_jcs_sig], axis=0)
    y     = np.concatenate([np.zeros(len(x_jet_bkg), dtype=np.float32),
                            np.ones (len(x_jet_sig), dtype=np.float32)], axis=0)
    return x_jet, x_jcs, y


def scale_jcs(x_jcs, scaler):
    N_CONSTIT_FEAT = x_jcs.shape[-1]
    shape = x_jcs.shape
    return scaler.transform(x_jcs.reshape(-1, N_CONSTIT_FEAT)).reshape(shape).astype(np.float32)


def roc_from_scores(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    area = sk_auc(fpr, tpr)
    if area < 0.5:
        fpr, tpr, area = 1 - fpr[::-1], 1 - tpr[::-1], 1 - area
    return fpr, tpr, area


def run_inference_with_jet(model, x_jet_test, x_jcs_test, mask_test, device, batch=512):
    x_jet_t = torch.from_numpy(x_jet_test)
    x_jcs_t = torch.from_numpy(x_jcs_test)
    mask_t  = torch.from_numpy(mask_test)
    logits_list = []
    model.eval()
    with torch.no_grad():
        for s in range(0, len(x_jet_test), batch):
            e = s + batch
            out = model(x_jet_t[s:e].to(device), x_jcs_t[s:e].to(device), mask_t[s:e].to(device))
            logits_list.append(out.cpu())
    logits = torch.cat(logits_list).numpy().ravel()
    return 1.0 / (1.0 + np.exp(-logits))


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

    # ── Load data once ────────────────────────────────────────────────────────
    print('Loading data ...')
    x_jet_all, x_jcs_all, y_all = load_data()

    # Use deepsets split (identical to all others: same seed + same data)
    idx_test = np.load(HERE / 'deepsets' / 'split_indices.npz')['test_idx']
    x_jet_test_raw = x_jet_all[idx_test]
    x_jcs_test_raw = x_jcs_all[idx_test]
    y_test         = y_all[idx_test]
    mask_test      = (x_jcs_test_raw[..., 3] > 0)
    print(f'Test set: {len(y_test)} events  (sig={int(y_test.sum())}, bkg={int((y_test==0).sum())})')

    results = []

    # ── 1. NNpull (DeepSets + jet) ────────────────────────────────────────────
    print('\n--- NNpull (DeepSets + jet) ---')
    nnpull_dir = ROOT / 'NNpull'
    ds_mod = load_module('nnpull_model', nnpull_dir / 'model.py')
    with open(nnpull_dir / 'scaler.pkl', 'rb') as fh:
        sc = pickle.load(fh)
    x_jet_s = sc['jet'].transform(x_jet_test_raw).astype(np.float32)
    x_jcs_s = scale_jcs(x_jcs_test_raw, sc['jcs'])
    model = ds_mod.JetNN().to(device)
    model.load_state_dict(torch.load(nnpull_dir / 'best_model.pt', map_location=device))
    scores = run_inference_with_jet(model, x_jet_s, x_jcs_s, mask_test, device)
    fpr, tpr, auc = roc_from_scores(y_test, scores)
    print(f'  AUC = {auc:.4f}')
    results.append(('NNpull (DeepSets + jet)', fpr, tpr, auc, 'black', '-', 3))

    # ── 2. NNtrans (Transformer + jet) ───────────────────────────────────────
    print('--- NNtrans (Transformer + jet) ---')
    nntrans_dir = ROOT / 'NNtrans'
    tr_mod = load_module('nntrans_model', nntrans_dir / 'model.py')
    with open(nntrans_dir / 'scaler.pkl', 'rb') as fh:
        sc = pickle.load(fh)
    x_jet_s = sc['jet'].transform(x_jet_test_raw).astype(np.float32)
    x_jcs_s = scale_jcs(x_jcs_test_raw, sc['jcs'])
    model = tr_mod.TransJetNN().to(device)
    model.load_state_dict(torch.load(nntrans_dir / 'best_model.pt', map_location=device))
    scores = run_inference_with_jet(model, x_jet_s, x_jcs_s, mask_test, device)
    fpr, tpr, auc = roc_from_scores(y_test, scores)
    print(f'  AUC = {auc:.4f}')
    results.append(('NNtrans (Trans + jet)', fpr, tpr, auc, '#1f77b4', '-', 3))

    # ── 3. DeepSets (no jet) ─────────────────────────────────────────────────
    print('--- DeepSets (no jet) ---')
    ds_nj_dir = HERE / 'deepsets'
    ds_nj_mod = load_module('deepsets_nojet_model', ds_nj_dir / 'model.py')
    with open(ds_nj_dir / 'scaler.pkl', 'rb') as fh:
        sc = pickle.load(fh)
    x_jcs_s = scale_jcs(x_jcs_test_raw, sc['jcs'])
    model = ds_nj_mod.DeepSetsNoJet().to(device)
    model.load_state_dict(torch.load(ds_nj_dir / 'best_model.pt', map_location=device))
    scores = run_inference_no_jet(model, x_jcs_s, mask_test, device)
    fpr, tpr, auc = roc_from_scores(y_test, scores)
    print(f'  AUC = {auc:.4f}')
    results.append(('DeepSets (no jet)', fpr, tpr, auc, '#d62728', '--', 2))

    # ── 4. Transformer (no jet) ───────────────────────────────────────────────
    print('--- Transformer (no jet) ---')
    tr_nj_dir = HERE / 'transformer'
    tr_nj_mod = load_module('trans_nojet_model', tr_nj_dir / 'model.py')
    with open(tr_nj_dir / 'scaler.pkl', 'rb') as fh:
        sc = pickle.load(fh)
    x_jcs_s = scale_jcs(x_jcs_test_raw, sc['jcs'])
    model = tr_nj_mod.TransNoJet().to(device)
    model.load_state_dict(torch.load(tr_nj_dir / 'best_model.pt', map_location=device))
    scores = run_inference_no_jet(model, x_jcs_s, mask_test, device)
    fpr, tpr, auc = roc_from_scores(y_test, scores)
    print(f'  AUC = {auc:.4f}')
    results.append(('Trans (no jet)', fpr, tpr, auc, '#ff7f0e', '--', 2))

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
    out = HERE / 'roc_comparison.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSaved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
