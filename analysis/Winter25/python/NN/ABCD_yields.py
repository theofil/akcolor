import sys
import json
import numpy as np
import uproot
import awkward as ak
import torch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from model import MLP

NN_DIR   = Path(__file__).resolve().parent
PAPER    = NN_DIR / 'paper'
MODELS   = NN_DIR / 'models'
DATA_DIR = NN_DIR.parents[1]
SCRIPT   = Path(__file__).stem

BKG_FILE = DATA_DIR / 'QCDHtoInv.root'
SIG_FILE = DATA_DIR / 'VBFHtoInv.root'

arch         = json.load(open(NN_DIR / 'architecture.json'))
FEATURE_SETS = arch['feature_sets']
MJJ_CUT      = arch['mjj_cut']

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

LUMI = 1e5  # 100 fb^{-1} in pb^{-1}

# region boundaries (same as ABCD_spva.py)
NN_LOW_MAX    = 0.65
NN_HIGH_MIN   = 0.65
SPVA_LOW_MAX  = 0.9
SPVA_HIGH_MIN = 2.0


def load_file(path):
    with uproot.open(str(path)) as f:
        flat   = f['events'].arrays(['mjj', 'dYjj', 'dPhijj', 'ptjj', 'kWeight'], library='np')
        jagged = f['events'].arrays(['jetSPVA', 'jetPVM'], library='ak')
    spva  = ak.pad_none(jagged['jetSPVA'], 2, axis=1)
    pvm   = ak.pad_none(jagged['jetPVM'],  2, axis=1)
    spva0 = np.abs(ak.to_numpy(ak.fill_none(spva[:, 0], np.nan)))
    spva1 = np.abs(ak.to_numpy(ak.fill_none(spva[:, 1], np.nan)))
    pvm0  = ak.to_numpy(ak.fill_none(pvm[:, 0], np.nan))
    pvm1  = ak.to_numpy(ak.fill_none(pvm[:, 1], np.nan))
    X    = np.column_stack([
        flat['mjj'], np.abs(flat['dYjj']), np.abs(flat['dPhijj']), flat['ptjj'],
        spva0, spva1, pvm0, pvm1,
    ])
    w    = flat['kWeight']
    mask = (flat['mjj'] > MJJ_CUT) & ~np.isnan(X).any(axis=1)
    return X[mask].astype(np.float32), w[mask]


def get_scores(key, X):
    cfg  = FEATURE_SETS[key]
    cols = cfg['cols']
    norm = np.load(NN_DIR / f'norm_{key}.npz')
    mu, std = norm['mu'], norm['std']
    Xk  = (X[:, cols] - mu) / std
    Xt  = torch.from_numpy(Xk.astype(np.float32)).to(device)
    model = MLP(len(cols), arch['hidden_layers'], arch['dropout'], arch['batch_norm']).to(device)
    model.load_state_dict(torch.load(MODELS / f'model_{key}.pt',
                                     map_location=device, weights_only=True))
    model.eval()
    with torch.no_grad():
        return model(Xt).cpu().numpy()


def region_masks(sc, spva0):
    return {
        'A': (sc < NN_LOW_MAX)  & (spva0 > SPVA_HIGH_MIN),
        'B': (sc > NN_HIGH_MIN) & (spva0 > SPVA_HIGH_MIN),
        'C': (sc < NN_LOW_MAX)  & (spva0 < SPVA_LOW_MAX),
        'D': (sc > NN_HIGH_MIN) & (spva0 < SPVA_LOW_MAX),
    }


def weighted_yield(w, mask):
    wm = w[mask]
    return wm.sum(), np.sqrt((wm**2).sum())


def fmt(val, err):
    """Format value ± error with automatic decimal places."""
    if val == 0:
        return r'$0 \pm 0$'
    exp = int(np.floor(np.log10(abs(val))))
    if -2 <= exp <= 4:
        dec = max(0, 2 - exp)
        return rf'${val:.{dec}f} \pm {err:.{dec}f}$'
    man_v = val / 10**exp
    man_e = err / 10**exp
    return rf'${man_v:.2f} \pm {man_e:.2f} \times 10^{{{exp}}}$'


print('Loading data...')
X_bkg, w_bkg = load_file(BKG_FILE)
X_sig, w_sig = load_file(SIG_FILE)

print('Running NN inference...')
sc_bkg = get_scores('a', X_bkg)
sc_sig = get_scores('a', X_sig)

masks_bkg = region_masks(sc_bkg, X_bkg[:, 4])
masks_sig = region_masks(sc_sig, X_sig[:, 4])

# ── compute yields ─────────────────────────────────────────────────────────────
region_desc = {
    'A': r'A (low NN, high $|\theta_s|$)',
    'B': r'B (high NN, high $|\theta_s|$)',
    'C': r'C (low NN, low $|\theta_s|$)',
    'D': r'D (high NN, low $|\theta_s|$)',
}

print(f'\nYields for {LUMI/1e3:.0f} fb⁻¹:')
print(f'  {"Region":<6}  {"QCD h":>12}  {"VBF h":>12}  {"Total":>12}')

rows = {}
for r in ('A', 'B', 'C', 'D'):
    y_qcd, e_qcd = weighted_yield(w_bkg, masks_bkg[r])
    y_vbf, e_vbf = weighted_yield(w_sig, masks_sig[r])
    y_qcd, e_qcd = y_qcd * LUMI, e_qcd * LUMI
    y_vbf, e_vbf = y_vbf * LUMI, e_vbf * LUMI
    y_tot = y_qcd + y_vbf
    e_tot = np.sqrt(e_qcd**2 + e_vbf**2)
    rows[r] = (y_qcd, e_qcd, y_vbf, e_vbf, y_tot, e_tot)
    print(f'  {r:<6}  {y_qcd:>12.4g}  {y_vbf:>12.4g}  {y_tot:>12.4g}')

# ── build LaTeX table ──────────────────────────────────────────────────────────
tex_rows = []
for r in ('A', 'B', 'C', 'D'):
    y_qcd, e_qcd, y_vbf, e_vbf, y_tot, e_tot = rows[r]
    tex_rows.append(
        f'    {region_desc[r]} & {fmt(y_qcd, e_qcd)} & {fmt(y_vbf, e_vbf)} & {fmt(y_tot, e_tot)} \\\\'
    )

lumi_label = f'{LUMI/1e3:.0f}'
block = '\n'.join([
    '% BEGIN_ABCD_YIELDS',
    r'\appendix',
    r'\section{ABCD Region Yields}',
    r'',
    r'\begin{table}[h]',
    r'  \centering',
    rf'  \caption{{Expected yields for {lumi_label}\,fb$^{{-1}}$ of luminosity in the ABCD regions for',
    r'           QCD $h$ and VBF $h$ processes. Statistical uncertainties reflect the MC sample size.}',
    r'  \label{tab:abcd_yields}',
    r'  \begin{tabular}{lccc}',
    r'    \hline',
    r'    Region & QCD $h$ & VBF $h$ & Total \\',
    r'    \hline',
    *tex_rows,
    r'    \hline',
    r'  \end{tabular}',
    r'\end{table}',
    '% END_ABCD_YIELDS',
])

# ── inject into main.tex ───────────────────────────────────────────────────────
import re

main_tex = PAPER / 'main.tex'
tex_src  = main_tex.read_text()

if '% BEGIN_ABCD_YIELDS' in tex_src:
    tex_src = re.sub(
        r'% BEGIN_ABCD_YIELDS.*?% END_ABCD_YIELDS',
        lambda _: block,
        tex_src,
        flags=re.DOTALL,
    )
    print('\nReplaced existing appendix block in main.tex')
else:
    tex_src = tex_src.replace(r'\bibliographystyle', block + '\n\n' + r'\bibliographystyle')
    print(f'\nAppendix block injected into {main_tex}')

main_tex.write_text(tex_src)
print('Done.')
