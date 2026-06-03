"""
Shared helpers for the DeepSetsNoScale interpretability pipeline.
All four interpret scripts (phi_visualize, phi_gradient, phi_pysr,
observable_roc) import from this module.
"""

import sys
import pathlib
import pickle

import numpy as np
import h5py
import torch
from sklearn.metrics import roc_curve, auc as sk_auc

HERE   = pathlib.Path(__file__).resolve().parent   # testAUC/
ROOT   = HERE.parent                               # makefriends/C/
DS_DIR = HERE / 'deepsets_noscale'
BKG_H5 = ROOT / 'friends' / 'QCDHtoInv.h5'
SIG_H5 = ROOT / 'friends' / 'VBFHtoInv.h5'

sys.path.insert(0, str(DS_DIR))
from model import DeepSetsNoScale                          # noqa: E402
from style import apply_style, FIG_SIZE, LABEL_SIZE, LEGEND_SIZE, TICK_SIZE  # noqa: E402


# ── model / scaler loading ────────────────────────────────────────────────────

def load_model(device):
    m = DeepSetsNoScale().to(device)
    m.load_state_dict(torch.load(DS_DIR / 'best_model.pt', map_location=device))
    m.eval()
    return m


def load_scaler():
    return pickle.load(open(DS_DIR / 'scaler.pkl', 'rb'))['jcs']


def scale_jcs(x_raw, scaler):
    shape = x_raw.shape
    return scaler.transform(x_raw.reshape(-1, 3)).reshape(shape).astype(np.float32)


# ── data loading ──────────────────────────────────────────────────────────────

def _read_3feat(path):
    """Return (N, 80, 3) array: [jcsDEta, jcsDPhi, jcsPt/jetPt]."""
    with h5py.File(path, 'r') as f:
        jetPt = f['jetPt'][:, 0].astype(np.float32)
        dEta  = f['jcsDEta'][:, 0, :].astype(np.float32)
        dPhi  = f['jcsDPhi'][:, 0, :].astype(np.float32)
        jcsPt = f['jcsPt'][:, 0, :].astype(np.float32)
    return np.stack([dEta, dPhi, jcsPt / jetPt[:, None]], axis=2)


def load_raw_data():
    """Returns x_all (N_all, 80, 3) in physical units and y_all (N_all,)."""
    x_bkg = _read_3feat(BKG_H5)
    x_sig = _read_3feat(SIG_H5)
    x_all = np.concatenate([x_bkg, x_sig])
    y_all = np.concatenate([np.zeros(len(x_bkg), np.float32),
                            np.ones(len(x_sig),  np.float32)])
    return x_all, y_all


def load_test_set():
    """Returns x_raw, x_scaled, mask, y — all for the test split."""
    x_all, y_all = load_raw_data()
    scaler = load_scaler()
    idx    = np.load(DS_DIR / 'split_indices.npz')['test_idx']
    x_raw  = x_all[idx]
    mask   = (x_raw[..., 2] > 0).astype(np.float32)   # jcsPt/jetPt > 0
    x_scal = scale_jcs(x_raw, scaler)
    return x_raw, x_scal, mask, y_all[idx]


def load_jet_features_test():
    """Returns jetPVM (N_test,) and jetSPVA (N_test,) for the test split."""
    def _read(path):
        with h5py.File(path, 'r') as f:
            pvm  = f['jetPVM'][:, 0].astype(np.float32)
            spva = f['jetSPVA'][:, 0].astype(np.float32)
        return pvm, spva

    pvm_b, spva_b = _read(BKG_H5)
    pvm_s, spva_s = _read(SIG_H5)
    pvm  = np.concatenate([pvm_b,  pvm_s])
    spva = np.concatenate([spva_b, spva_s])
    idx  = np.load(DS_DIR / 'split_indices.npz')['test_idx']
    return pvm[idx], spva[idx]


# ── network helpers ───────────────────────────────────────────────────────────

def compute_phi_aggregates(model, x_scaled, mask, device, batch=512):
    """Masked-sum-pooled phi outputs. Returns (N, 64) on CPU (no_grad)."""
    out = []
    with torch.no_grad():
        for s in range(0, len(x_scaled), batch):
            xb = torch.from_numpy(x_scaled[s:s + batch]).to(device)
            mb = torch.from_numpy(mask[s:s + batch]).to(device)
            B, N, F = xb.shape
            h = model.phi(xb.view(B * N, F)).view(B, N, 64)
            out.append((h * mb.unsqueeze(-1)).sum(1).cpu())
    return torch.cat(out)   # (N, 64)


def compute_mean_rho_jacobian(model, h_test, device, batch=512):
    """Mean d(rho(h))/d(h) over h_test. Returns np.ndarray (64,).
    model must already be in eval() mode."""
    w = torch.zeros(64, dtype=torch.float64)
    for s in range(0, len(h_test), batch):
        hb = h_test[s:s + batch].to(device).requires_grad_(True)
        model.rho(hb).sum().backward()          # gradient flows only through rho
        w += hb.grad.detach().cpu().double().sum(0)
    return (w / len(h_test)).float().numpy()    # (64,)


# ── polynomial basis ──────────────────────────────────────────────────────────

TERM_NAMES = [
    '1',
    'dEta', 'dPhi', 'z',
    'r',
    'dEta*dPhi',
    'dEta^2', 'dPhi^2', 'z^2', 'r^2',
    'z*r', 'z*dEta', 'z*dPhi',
    'r^3', 'z*r^2', 'z^2*r',
    'r^4', 'z^2*r^2',
]   # 18 terms


def make_poly_basis(X):
    """18-term polynomial basis. X: (M, 3) with columns [dEta, dPhi, z].
    Includes a constant column so use Ridge(fit_intercept=False)."""
    dEta, dPhi, z = X[:, 0], X[:, 1], X[:, 2]
    r = np.sqrt(dEta**2 + dPhi**2)
    return np.column_stack([
        np.ones_like(z),                               #  0: 1
        dEta, dPhi, z,                                 #  1-3
        r,                                              #  4: r
        dEta * dPhi,                                    #  5
        dEta**2, dPhi**2, z**2, r**2,                  #  6-9
        z * r,   z * dEta, z * dPhi,                   # 10-12  IRC-safe
        r**3, z * r**2, z**2 * r,                      # 13-15
        r**4, z**2 * r**2,                             # 16-17
    ])


# ── ROC helper ────────────────────────────────────────────────────────────────

def roc_from_scores(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    area = sk_auc(fpr, tpr)
    if area < 0.5:
        fpr, tpr, area = 1 - fpr[::-1], 1 - tpr[::-1], 1 - area
    return fpr, tpr, area
