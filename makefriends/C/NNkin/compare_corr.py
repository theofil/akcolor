#!/usr/bin/env python3
"""
2D density plot of NNkin score vs NNpull score, separately for QCD and VBF,
with Pearson ρ annotated on each panel.

Requires both NNkin and NNpull to have been trained (best_model.pt, scaler.pkl,
split_indices.npz must exist in their respective directories).

Output:
  compare_corr.pdf  — saved next to this script
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

HERE   = pathlib.Path(__file__).parent
NNPULL = HERE.parent / 'NNpull'

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

pull_ds  = _load_module('pull_dataset', NNPULL / 'dataset.py')

_prev_dataset = sys.modules.get('dataset')
sys.modules['dataset'] = pull_ds
pull_mod = _load_module('pull_model', NNPULL / 'model.py')
if _prev_dataset is None:
    del sys.modules['dataset']
else:
    sys.modules['dataset'] = _prev_dataset

BKG_FILE = HERE / kin_ds.BKG_FILE
SIG_FILE = HERE / kin_ds.SIG_FILE


def scale_jcs(x_jcs, scaler):
    shape = x_jcs.shape
    return scaler.transform(x_jcs.reshape(-1, pull_ds.N_CONSTIT_FEAT)).reshape(shape).astype(np.float32)


def run_kin(x_all_raw, idx_test, device):
    x_test_raw = x_all_raw[idx_test]
    with open(HERE / 'scaler.pkl', 'rb') as fh:
        scaler = pickle.load(fh)
    x_test = scaler.transform(x_test_raw).astype(np.float32)
    model = kin_mod.KinNN().to(device)
    model.load_state_dict(torch.load(HERE / 'best_model.pt', map_location=device))
    model.eval()
    with torch.no_grad():
        logits = model(torch.from_numpy(x_test).to(device)).cpu().numpy().ravel()
    return 1.0 / (1.0 + np.exp(-logits))


def run_pull(x_jet_all, x_jcs_all, idx_test, device):
    x_jet_raw = x_jet_all[idx_test]
    x_jcs_raw = x_jcs_all[idx_test]
    with open(NNPULL / 'scaler.pkl', 'rb') as fh:
        scalers = pickle.load(fh)
    x_jet = scalers['jet'].transform(x_jet_raw).astype(np.float32)
    x_jcs = scale_jcs(x_jcs_raw, scalers['jcs'])
    mask  = (x_jcs_raw[..., 3] > 0)   # jcsPt > 0 = valid constituent
    model = pull_mod.JetNN().to(device)
    model.load_state_dict(torch.load(NNPULL / 'best_model.pt', map_location=device))
    model.eval()
    with torch.no_grad():
        logits = model(
            torch.from_numpy(x_jet).to(device),
            torch.from_numpy(x_jcs).to(device),
            torch.from_numpy(mask).to(device),
        ).cpu().numpy().ravel()
    return 1.0 / (1.0 + np.exp(-logits))


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    print('Loading kinematic features ...')
    x_kin_bkg, _ = kin_ds.load_features(BKG_FILE)
    x_kin_sig, _ = kin_ds.load_features(SIG_FILE)
    x_kin_all = np.concatenate([x_kin_bkg, x_kin_sig], axis=0)

    print('Loading pull features ...')
    x_jet_bkg, x_jcs_bkg, _ = pull_ds.load_features(BKG_FILE)
    x_jet_sig, x_jcs_sig, _ = pull_ds.load_features(SIG_FILE)
    x_jet_all = np.concatenate([x_jet_bkg, x_jet_sig], axis=0)
    x_jcs_all = np.concatenate([x_jcs_bkg, x_jcs_sig], axis=0)

    y_all = np.concatenate([np.zeros(len(x_kin_bkg), dtype=np.float32),
                            np.ones (len(x_kin_sig), dtype=np.float32)], axis=0)

    idx_kin  = np.load(HERE   / 'split_indices.npz')['test_idx']
    idx_pull = np.load(NNPULL / 'split_indices.npz')['test_idx']
    assert np.array_equal(np.sort(idx_kin), np.sort(idx_pull)), \
        'NNkin and NNpull test indices differ — were they trained with the same seed?'

    y_test = y_all[idx_kin]

    print('Running NNkin inference ...')
    sc_kin  = run_kin(x_kin_all, idx_kin, device)

    print('Running NNpull inference ...')
    sc_pull = run_pull(x_jet_all, x_jcs_all, idx_kin, device)

    mask_qcd = (y_test == 0)
    mask_vbf = (y_test == 1)

    rho_qcd = np.corrcoef(sc_kin[mask_qcd], sc_pull[mask_qcd])[0, 1]
    rho_vbf = np.corrcoef(sc_kin[mask_vbf], sc_pull[mask_vbf])[0, 1]
    print(f'Pearson rho  QCD = {rho_qcd:.4f}')
    print(f'Pearson rho  VBF = {rho_vbf:.4f}')

    from style import apply_style, LABEL_SIZE, TICK_SIZE, TITLE_SIZE, LEGEND_SIZE

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    for ax, mask, title, rho in [
        (axes[0], mask_qcd, 'QCD background', rho_qcd),
        (axes[1], mask_vbf, 'VBF signal',     rho_vbf),
    ]:
        hb = ax.hexbin(sc_kin[mask], sc_pull[mask],
                       gridsize=60, bins='log', cmap='viridis',
                       extent=[0, 1, 0, 1])
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label('counts (log scale)', fontsize=LEGEND_SIZE)
        cb.ax.tick_params(labelsize=TICK_SIZE - 2)
        ax.text(0.05, 0.93, rf'$\rho = {rho:.3f}$',
                transform=ax.transAxes, fontsize=LEGEND_SIZE,
                verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('NNkin score', fontsize=LABEL_SIZE)
        ax.set_ylabel('NNpull score', fontsize=LABEL_SIZE)
        ax.set_title(title, fontsize=TITLE_SIZE)
        ax.tick_params(axis='both', labelsize=TICK_SIZE)
        ax.grid(True, linestyle=':', alpha=0.4)

    plt.tight_layout()
    out = HERE / 'compare_corr.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
