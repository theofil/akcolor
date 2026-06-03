#!/usr/bin/env python3
"""
2D density plot of NNpull(jet 0) vs NNpull(jet 1) for mjj > 400 GeV events,
separately for QCD and VBF, with Pearson rho annotated on each panel.

Output:
  compare_corr_PullOnPull.pdf  — saved next to this script
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

# ── import NNpull modules ─────────────────────────────────────────────────────
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
    sys.modules.pop('dataset', None)
else:
    sys.modules['dataset'] = _prev_dataset

BKG_FILE = (NNPULL / pull_ds.BKG_FILE).resolve()
SIG_FILE = (NNPULL / pull_ds.SIG_FILE).resolve()

JET_FEATURES     = pull_ds.JET_FEATURES
CONSTIT_FEATURES = pull_ds.CONSTIT_FEATURES
N_CONSTIT_FEAT   = pull_ds.N_CONSTIT_FEAT

MJJ_CUT = 400.0   # GeV


def load_features_for_jet(h5path, jet_idx):
    """Return x_jet (N,6), x_jcs (N,80,5), mjj (N,) for the given jet index."""
    import h5py
    with h5py.File(h5path, 'r') as f:
        jets = {k: f[k][:, jet_idx].astype(np.float32) for k in JET_FEATURES}
        jcs  = {k: f[k][:, jet_idx, :].astype(np.float32) for k in CONSTIT_FEATURES}
        mjj  = f['mjj'][:].astype(np.float32)
    jets['jetEta'] = np.abs(jets['jetEta'])
    x_jet = np.stack([jets[k] for k in JET_FEATURES],     axis=1)   # (N, 6)
    x_jcs = np.stack([jcs[k]  for k in CONSTIT_FEATURES], axis=2)   # (N, 80, 5)
    return x_jet, x_jcs, mjj


def run_pull_inference(x_jet_raw, x_jcs_raw, scalers, model, device, batch=512):
    """Normalise features and run NNpull; return sigmoid scores."""
    x_jet = scalers['jet'].transform(x_jet_raw).astype(np.float32)
    shape = x_jcs_raw.shape
    x_jcs = scalers['jcs'].transform(
        x_jcs_raw.reshape(-1, N_CONSTIT_FEAT)
    ).reshape(shape).astype(np.float32)
    mask = (x_jcs_raw[..., 3] > 0)   # jcsPt > 0 — valid constituent

    x_jet_t = torch.from_numpy(x_jet)
    x_jcs_t = torch.from_numpy(x_jcs)
    mask_t  = torch.from_numpy(mask)

    logits_list = []
    with torch.no_grad():
        for start in range(0, len(x_jet), batch):
            end = start + batch
            out = model(
                x_jet_t[start:end].to(device),
                x_jcs_t[start:end].to(device),
                mask_t [start:end].to(device),
            )
            logits_list.append(out.cpu())
    logits = torch.cat(logits_list).numpy().ravel()
    return 1.0 / (1.0 + np.exp(-logits))


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    # ── Load features for both jets ───────────────────────────────────────────
    print('Loading jet-0 features ...')
    x_jet0_bkg, x_jcs0_bkg, mjj_bkg = load_features_for_jet(BKG_FILE, 0)
    x_jet0_sig, x_jcs0_sig, mjj_sig = load_features_for_jet(SIG_FILE, 0)

    print('Loading jet-1 features ...')
    x_jet1_bkg, x_jcs1_bkg, _       = load_features_for_jet(BKG_FILE, 1)
    x_jet1_sig, x_jcs1_sig, _       = load_features_for_jet(SIG_FILE, 1)

    x_jet0_all = np.concatenate([x_jet0_bkg, x_jet0_sig], axis=0)
    x_jcs0_all = np.concatenate([x_jcs0_bkg, x_jcs0_sig], axis=0)
    x_jet1_all = np.concatenate([x_jet1_bkg, x_jet1_sig], axis=0)
    x_jcs1_all = np.concatenate([x_jcs1_bkg, x_jcs1_sig], axis=0)
    mjj_all    = np.concatenate([mjj_bkg,     mjj_sig],    axis=0)
    y_all      = np.concatenate([
        np.zeros(len(x_jet0_bkg), dtype=np.float32),
        np.ones (len(x_jet0_sig), dtype=np.float32),
    ], axis=0)

    # ── Apply test split ──────────────────────────────────────────────────────
    idx_test = np.load(NNPULL / 'split_indices.npz')['test_idx']

    x_jet0_test = x_jet0_all[idx_test]
    x_jcs0_test = x_jcs0_all[idx_test]
    x_jet1_test = x_jet1_all[idx_test]
    x_jcs1_test = x_jcs1_all[idx_test]
    mjj_test    = mjj_all[idx_test]
    y_test      = y_all[idx_test]
    print(f'Test events before mjj cut: {len(y_test)}')

    # ── Apply mjj > 400 GeV cut ───────────────────────────────────────────────
    mjj_mask = mjj_test > MJJ_CUT
    x_jet0_test = x_jet0_test[mjj_mask]
    x_jcs0_test = x_jcs0_test[mjj_mask]
    x_jet1_test = x_jet1_test[mjj_mask]
    x_jcs1_test = x_jcs1_test[mjj_mask]
    y_test       = y_test[mjj_mask]
    print(f'Test events after  mjj>{MJJ_CUT:.0f} GeV cut: {len(y_test)}  '
          f'(sig={int(y_test.sum())}, bkg={int((y_test==0).sum())})')

    # ── Load NNpull model and scalers ─────────────────────────────────────────
    with open(NNPULL / 'scaler.pkl', 'rb') as fh:
        scalers = pickle.load(fh)

    model = pull_mod.JetNN().to(device)
    model.load_state_dict(torch.load(NNPULL / 'best_model.pt', map_location=device))
    model.eval()

    # ── Run inference ─────────────────────────────────────────────────────────
    print('Running NNpull inference on jet 0 ...')
    sc_jet0 = run_pull_inference(x_jet0_test, x_jcs0_test, scalers, model, device)

    print('Running NNpull inference on jet 1 ...')
    sc_jet1 = run_pull_inference(x_jet1_test, x_jcs1_test, scalers, model, device)

    # ── Pearson correlations ──────────────────────────────────────────────────
    mask_qcd = (y_test == 0)
    mask_vbf = (y_test == 1)

    rho_qcd = np.corrcoef(sc_jet0[mask_qcd], sc_jet1[mask_qcd])[0, 1]
    rho_vbf = np.corrcoef(sc_jet0[mask_vbf], sc_jet1[mask_vbf])[0, 1]
    print(f'Pearson rho  QCD = {rho_qcd:.4f}')
    print(f'Pearson rho  VBF = {rho_vbf:.4f}')

    # ── Plot ──────────────────────────────────────────────────────────────────
    from style import LABEL_SIZE, TICK_SIZE, TITLE_SIZE, LEGEND_SIZE

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    for ax, mask, title, rho in [
        (axes[0], mask_qcd, 'QCD background', rho_qcd),
        (axes[1], mask_vbf, 'VBF signal',     rho_vbf),
    ]:
        hb = ax.hexbin(sc_jet0[mask], sc_jet1[mask],
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
        ax.set_xlabel('NNpull score (jet 0)', fontsize=LABEL_SIZE)
        ax.set_ylabel('NNpull score (jet 1)', fontsize=LABEL_SIZE)
        ax.set_title(title, fontsize=TITLE_SIZE)
        ax.tick_params(axis='both', labelsize=TICK_SIZE)
        ax.grid(True, linestyle=':', alpha=0.4)

    fig.suptitle(rf'$m_{{jj}} > {MJJ_CUT:.0f}$ GeV', fontsize=TITLE_SIZE, y=1.01)
    plt.tight_layout()
    out = HERE / 'compare_corr_PullOnPull.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
