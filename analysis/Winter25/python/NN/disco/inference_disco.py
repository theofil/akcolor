import sys
import json
import numpy as np
import uproot
import awkward as ak
import torch
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

DISCO_DIR    = Path(__file__).resolve().parent
NN_DIR       = DISCO_DIR.parent
DATA_DIR     = NN_DIR.parents[1]
FIGS         = NN_DIR / 'paper' / 'figs'
MODELS       = NN_DIR / 'models'
MODELS_DISCO = DISCO_DIR / 'models'

sys.path.insert(0, str(NN_DIR))
from model import MLP
from style import apply_style, FIG_SIZE

arch    = json.load(open(NN_DIR / 'architecture.json'))
MJJ_CUT = arch['mjj_cut']
COLS_D  = [6, 7, 8, 9, 10, 11, 12, 13]

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

MJJ_SLICES   = [(400, 550), (550, 700), (700, 900), (900, None)]
SLICE_COLORS = ['black', 'steelblue', 'darkgreen', 'crimson']


# ── data loading (for inference ROC: full 18-col matrix + SPVA) ──────────────
def load_file_full(path):
    print(f'Loading {path.name}...')
    with uproot.open(str(path)) as f:
        flat   = f['events'].arrays(['mjj', 'dYjj', 'dPhijj', 'ptjj', 'c21', 'c12'], library='np')
        jagged = f['events'].arrays(
            ['jetSPVA', 'jetPVM', 'jetPt', 'jetEta', 'jetPhi', 'jetM', 'jetPVA'], library='ak')

    def pad2(arr): return ak.pad_none(arr, 2, axis=1)
    def col(arr, i): return ak.to_numpy(ak.fill_none(pad2(arr)[:, i], np.nan))

    spva0_abs = np.abs(col(jagged['jetSPVA'], 0))

    X = np.column_stack([
        flat['mjj'], flat['dYjj'], flat['dPhijj'], flat['ptjj'],
        col(jagged['jetPVM'], 0), col(jagged['jetPVM'], 1),
        col(jagged['jetPt'],  0), col(jagged['jetPt'],  1),
        col(jagged['jetEta'], 0), col(jagged['jetEta'], 1),
        col(jagged['jetPhi'], 0), col(jagged['jetPhi'], 1),
        col(jagged['jetM'],   0), col(jagged['jetM'],   1),
        flat['c21'], flat['c12'],
        col(jagged['jetPVA'], 0), col(jagged['jetPVA'], 1),
    ])
    mask = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X).any(axis=1)
    return X[mask].astype(np.float32), spva0_abs[mask].astype(np.float32)


print('Loading data...')
X_bkg, spva0_bkg = load_file_full(DATA_DIR / 'QCDHtoInv.root')
X_sig, spva0_sig = load_file_full(DATA_DIR / 'VBFHtoInv.root')

X     = np.concatenate([X_bkg,     X_sig],     axis=0)
spva0 = np.concatenate([spva0_bkg, spva0_sig], axis=0)
y     = np.concatenate([np.zeros(len(X_bkg), dtype=np.float32),
                        np.ones (len(X_sig), dtype=np.float32)])

split   = np.load(NN_DIR / 'split_indices.npz')
idx_inf = split['idx_inf']

X_inf     = X[idx_inf]
spva0_inf = spva0[idx_inf]
y_inf     = y[idx_inf]
print(f'Inference set: {len(y_inf)} events  ({(y_inf==1).sum()} sig, {(y_inf==0).sum()} bkg)')


# ── scoring on inference split ────────────────────────────────────────────────
def score_model(model_path, norm_path):
    norm_data = np.load(norm_path)
    mu, std   = norm_data['mu'], norm_data['std']
    Xk = (X_inf[:, COLS_D] - mu) / std
    Xt = torch.from_numpy(Xk.astype(np.float32)).to(device)
    m  = MLP(8, arch['hidden_layers'], arch['dropout'], arch['batch_norm']).to(device)
    m.load_state_dict(torch.load(model_path, map_location=device, weights_only=True))
    m.eval()
    with torch.no_grad():
        return m(Xt).cpu().numpy()


print('Scoring baseline NNraw...')
scores_d     = score_model(MODELS / 'model_d.pt',            NN_DIR / 'norm_d.npz')
print('Scoring NNraw+DisCo...')
scores_disco = score_model(MODELS_DISCO / 'model_d_disco.pt', NN_DIR / 'norm_d.npz')


# ── mjj slice plots (all events from ROOT, not just inference split) ──────────
# Mirror NNpull_slices.py convention: use all data for better slice statistics.
def load_file_score(root_path, model_path, norm_path):
    norm_data = np.load(norm_path)
    mu, std   = norm_data['mu'], norm_data['std']
    with uproot.open(str(root_path)) as f:
        flat   = f['events'].arrays(['mjj', 'dYjj', 'dPhijj', 'ptjj', 'c21', 'c12'], library='np')
        jagged = f['events'].arrays(
            ['jetPVM', 'jetPt', 'jetEta', 'jetPhi', 'jetM', 'jetPVA'], library='ak')

    def pad2(arr): return ak.pad_none(arr, 2, axis=1)
    def col(arr, i): return ak.to_numpy(ak.fill_none(pad2(arr)[:, i], np.nan))

    X_full = np.column_stack([
        flat['mjj'], flat['dYjj'], flat['dPhijj'], flat['ptjj'],
        col(jagged['jetPVM'], 0), col(jagged['jetPVM'], 1),
        col(jagged['jetPt'],  0), col(jagged['jetPt'],  1),
        col(jagged['jetEta'], 0), col(jagged['jetEta'], 1),
        col(jagged['jetPhi'], 0), col(jagged['jetPhi'], 1),
        col(jagged['jetM'],   0), col(jagged['jetM'],   1),
        flat['c21'], flat['c12'],
        col(jagged['jetPVA'], 0), col(jagged['jetPVA'], 1),
    ])
    mask  = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X_full).any(axis=1)
    X_full = X_full[mask].astype(np.float32)
    mjj_all = X_full[:, 0]

    Xk = (X_full[:, COLS_D] - mu) / std
    Xt = torch.from_numpy(Xk.astype(np.float32)).to(device)
    m  = MLP(8, arch['hidden_layers'], arch['dropout'], arch['batch_norm']).to(device)
    m.load_state_dict(torch.load(model_path, map_location=device, weights_only=True))
    m.eval()
    with torch.no_grad():
        sc = m(Xt).cpu().numpy()
    return mjj_all, sc


bins        = np.linspace(0, 1, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_width   = bins[1] - bins[0]


def norm_hist(values):
    counts, _ = np.histogram(values, bins=bins)
    norm = counts.sum() * bin_width
    if norm > 0:
        return counts / norm, np.sqrt(counts) / norm
    return np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)


def slice_title(lo, hi):
    if hi is None:
        return fr'$m_{{jj}} > {lo}$ GeV'
    return fr'${lo} < m_{{jj}} < {hi}$ GeV'


def plot_mjj_slices(mjj_all, sc_all, title, outname):
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for color, (lo, hi) in zip(SLICE_COLORS, MJJ_SLICES):
        tag  = slice_title(lo, hi)
        mask = mjj_all > lo
        if hi is not None:
            mask = mask & (mjj_all < hi)
        h, herr = norm_hist(sc_all[mask])
        ax.hist(bin_centers, bins=bins, weights=h, histtype='step',
                linewidth=2.5, color=color, label=tag)
        ax.errorbar(bin_centers, h, yerr=herr, fmt='none',
                    color=color, capsize=3, capthick=1, elinewidth=1)
    apply_style(ax,
                xlabel='NN score',
                ylabel='Normalized density',
                title=title,
                xlim=(0, 1),
                legend_loc='upper left')
    ax.legend(loc='upper left', frameon=False, fontsize=11, ncol=1)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{outname}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {outname}.pdf/png')


print('Loading QCD events for slice plots...')
mjj_qcd_base,  sc_qcd_base  = load_file_score(
    DATA_DIR / 'QCDHtoInv.root', MODELS / 'model_d.pt',            NN_DIR / 'norm_d.npz')
mjj_qcd_disco, sc_qcd_disco = load_file_score(
    DATA_DIR / 'QCDHtoInv.root', MODELS_DISCO / 'model_d_disco.pt', NN_DIR / 'norm_d.npz')

print('Loading VBF events for slice plots...')
mjj_vbf_base,  sc_vbf_base  = load_file_score(
    DATA_DIR / 'VBFHtoInv.root', MODELS / 'model_d.pt',            NN_DIR / 'norm_d.npz')
mjj_vbf_disco, sc_vbf_disco = load_file_score(
    DATA_DIR / 'VBFHtoInv.root', MODELS_DISCO / 'model_d_disco.pt', NN_DIR / 'norm_d.npz')

plot_mjj_slices(mjj_qcd_base,  sc_qcd_base,
                r'NNraw (QCD $h$)',                      'disco_NNraw_slices_mjj_QCD')
plot_mjj_slices(mjj_qcd_disco, sc_qcd_disco,
                r'NNraw+DisCo (QCD $h$, $\lambda=10$)',  'disco_NNraw_disco_slices_mjj_QCD')
plot_mjj_slices(mjj_vbf_base,  sc_vbf_base,
                r'NNraw (VBF $h$)',                      'disco_NNraw_slices_mjj_VBF')
plot_mjj_slices(mjj_vbf_disco, sc_vbf_disco,
                r'NNraw+DisCo (VBF $h$, $\lambda=10$)',  'disco_NNraw_disco_slices_mjj_VBF')


# ── mjj distributions in NN score slices ──────────────────────────────────────
NN_SCORE_SLICES = [(0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0)]
SCORE_COLORS    = ['black', 'steelblue', 'darkgreen', 'crimson']

mjj_bins        = np.linspace(400, 2000, 21)
mjj_centers     = (mjj_bins[:-1] + mjj_bins[1:]) / 2
mjj_bin_width   = mjj_bins[1] - mjj_bins[0]


def norm_hist_mjj(values):
    counts, _ = np.histogram(values, bins=mjj_bins)
    norm = counts.sum() * mjj_bin_width
    if norm > 0:
        return counts / norm, np.sqrt(counts) / norm
    return np.zeros(len(mjj_bins) - 1), np.zeros(len(mjj_bins) - 1)


def score_slice_title(lo, hi):
    return fr'${lo:.2f} < \mathrm{{NN}} < {hi:.2f}$'


def plot_score_slices(mjj_all, sc_all, title, outname):
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for color, (lo, hi) in zip(SCORE_COLORS, NN_SCORE_SLICES):
        tag  = score_slice_title(lo, hi)
        mask = (sc_all >= lo) & (sc_all < hi)
        if mask.sum() < 2:
            continue
        h, herr = norm_hist_mjj(mjj_all[mask])
        ax.hist(mjj_centers, bins=mjj_bins, weights=h, histtype='step',
                linewidth=2.5, color=color, label=tag)
        ax.errorbar(mjj_centers, h, yerr=herr, fmt='none',
                    color=color, capsize=3, capthick=1, elinewidth=1)
    apply_style(ax,
                xlabel=r'$m_{jj}$ [GeV]',
                ylabel='Normalized density [1/GeV]',
                title=title,
                xlim=(400, 2000),
                legend_loc='upper right',
                log_y=True)
    ax.legend(loc='upper right', frameon=False, fontsize=11, ncol=1)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'{outname}.{ext}', bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {outname}.pdf/png')


plot_score_slices(mjj_qcd_base,  sc_qcd_base,
                  r'NNraw (QCD $h$)',                     'disco_NNraw_slices_score_QCD')
plot_score_slices(mjj_qcd_disco, sc_qcd_disco,
                  r'NNraw+DisCo (QCD $h$, $\lambda=10$)', 'disco_NNraw_disco_slices_score_QCD')
plot_score_slices(mjj_vbf_base,  sc_vbf_base,
                  r'NNraw (VBF $h$)',                     'disco_NNraw_slices_score_VBF')
plot_score_slices(mjj_vbf_disco, sc_vbf_disco,
                  r'NNraw+DisCo (VBF $h$, $\lambda=10$)', 'disco_NNraw_disco_slices_score_VBF')


# ── ROC comparison on inference split ─────────────────────────────────────────
discriminants = {
    r'NN$_\mathrm{raw}$':        scores_d,
    r'NN$_\mathrm{raw}$+DisCo':  scores_disco,
    r'$m_{jj}$':                  X_inf[:, 0],
    r'$|\theta_s[j_1]|$':         -spva0_inf,
}
colors = {
    r'NN$_\mathrm{raw}$':        'darkorchid',
    r'NN$_\mathrm{raw}$+DisCo':  'firebrick',
    r'$m_{jj}$':                  'darkorange',
    r'$|\theta_s[j_1]|$':         'purple',
}
linestyles = {
    r'NN$_\mathrm{raw}$':        'solid',
    r'NN$_\mathrm{raw}$+DisCo':  'solid',
    r'$m_{jj}$':                  'dashed',
    r'$|\theta_s[j_1]|$':         'dashed',
}

fig, ax = plt.subplots(figsize=FIG_SIZE)
for label, sc in discriminants.items():
    auc       = roc_auc_score(y_inf, sc)
    fpr, tpr, _ = roc_curve(y_inf, sc)
    ax.plot(fpr, tpr, lw=2, color=colors[label], linestyle=linestyles[label],
            label=f'{label}  AUC={auc:.3f}')

ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4)
ax.set_xlabel('False Positive Rate', fontsize=20)
ax.set_ylabel('True Positive Rate', fontsize=20)
ax.set_title(r'VBF $h$ vs QCD $h$', fontsize=22)
ax.legend(fontsize=13, loc='lower right', frameon=False)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.tick_params(axis='both', labelsize=16)
ax.grid(True, linestyle=':', alpha=0.4)
plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(FIGS / f'disco_inference_ROC_comparison.{ext}', bbox_inches='tight')
plt.close(fig)
print('Saved disco_inference_ROC_comparison.pdf/png')

print(f'\nAll figures saved to {FIGS}')
