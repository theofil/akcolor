#!/usr/bin/env python3
"""
Estimate the best achievable signal/background for every NN* variant directly
from its saved ROC curve (roc_data_{H,W,Z}.npz in each net directory), at
L = 300 fb^-1.

For each (net, channel, generator) combination, scan the net's own ROC curve
(tpr/fpr) and pick the operating point that maximizes significance = S/sqrt(S+B),
subject to a raw-MC-background statistical floor. S and B are formed from the
sample's post-selection sum(kWeight) (already in pb, see kWeight formula in
summer26.md) times L times the ROC efficiency at that point.

Branching-ratio rescale (decay-chain-aware, see summer26.md "Signal/background
estimation from ROC curves"): the samples were generated with H kept fully
stable (BSM-invisible convention -- tabulated sigma is the true inclusive
production rate), Z reconstructed as Z->nu_e nu_e-bar, and W as W->tau nu_tau
(real SM decay chains -- tabulated sigma already has that channel's real BR
baked in by the matrix element). To reinterpret these yields for a realistic
H->gammagamma / Z->ee,mumu / W->enu,munu search, each channel's S and B are
rescaled by the same factor (same boson, same decay requirement for signal and
background alike):

    H:  factor = BR(H->gammagamma)
    Z:  factor = [BR(Z->ee) + BR(Z->mumu)] / BR(Z->nu_e nu_e-bar)
    W:  factor = [BR(W->enu) + BR(W->munu)] / BR(W->taunu)

Because the same factor multiplies S and B within a channel, it cancels in the
S/B ratio and does not shift the argmax over the ROC curve -- the scan runs
once on un-rescaled yields, and the rescale is applied only to the winning
point's final S, B (see verify_rescale_invariance()).

Writes figs/summer26/SB_estimate.txt and prints a ready-to-paste markdown table.
"""

import pathlib
import re

import numpy as np
import h5py

HERE        = pathlib.Path(__file__).parent
FRIENDS_DIR = HERE.parent / 'friends' / 'summer26'
FIGS_DIR    = HERE.parent / 'figs' / 'summer26'
OUT_TXT     = FIGS_DIR / 'SB_estimate.txt'

L           = 300_000.0   # pb^-1 = 300 fb^-1
MIN_RAW_BKG = 10          # statistical-validity floor (repo precedent: old SR_optimization.txt)

NETS     = ['NNkin', 'NNj', 'NNjB', 'NNjj', 'NNjjB', 'NNjjBj']
CHANNELS = ['H', 'W', 'Z']
GENS     = ['hw', 'mg5']

CHANNEL_FILES = {
    'H': dict(bkg_hw='QCDHjj_herwig.h5',      sig_hw='VBFH_herwig.h5',
              bkg_mg5='QCDHjj_mg5_pythia.h5', sig_mg5='VBFH_mg5_pythia.h5'),
    'W': dict(bkg_hw='QCDWjj_herwig.h5',      sig_hw='VBFW_herwig.h5',
              bkg_mg5='QCDWjj_mg5_pythia.h5', sig_mg5='VBFW_mg5_pythia.h5'),
    'Z': dict(bkg_hw='QCDZjj_herwig.h5',      sig_hw='VBFZ_herwig.h5',
              bkg_mg5='QCDZjj_mg5_pythia.h5', sig_mg5='VBFZ_mg5_pythia.h5'),
}

# PDG branching ratios
BR_H_GAMMAGAMMA    = 2.27e-3                    # BR(H->gammagamma), mH = 125.09 GeV
BR_Z_EE, BR_Z_MUMU = 0.03363, 0.03366
BR_Z_NUE_PERFLAVOR = 0.2000 / 3.0               # BR(Z->nunu total) = 20.00%, / 3 flavors
BR_W_ENU, BR_W_MUNU, BR_W_TAUNU = 0.1071, 0.1063, 0.1138

RESCALE = {   # sigma_new = sigma_tabulated * RESCALE[channel]; applied to S and B alike
    'H': BR_H_GAMMAGAMMA,
    'Z': (BR_Z_EE + BR_Z_MUMU) / BR_Z_NUE_PERFLAVOR,
    'W': (BR_W_ENU + BR_W_MUNU) / BR_W_TAUNU,
}


def net_dir(net, channel):
    suffix = '' if channel == 'Z' else f'_{channel}'
    return HERE / f'{net}{suffix}'


def kweight_sum_and_n(h5path):
    """(sum(kWeight) [pb], row count) over the full h5 file -- no filtering."""
    with h5py.File(h5path, 'r') as f:
        kw = f['kWeight'][:].astype(np.float64)
    return float(kw.sum()), len(kw)


def herwig_test_split_counts(ndir, channel, n_bkg_full):
    """Raw (n_bkg_test, n_sig_test) from this net's held-out Herwig test split.

    inference.py concatenates [background, signal] before splitting
    (y_all_hw = concat([zeros(n_bkg), ones(n_sig)])), so a test index below
    n_bkg_full is background, else signal.
    """
    idx_test = np.load(ndir / f'split_indices_{channel}.npz')['test_idx']
    n_bkg_test = int((idx_test < n_bkg_full).sum())
    n_sig_test = int(len(idx_test) - n_bkg_test)
    return n_bkg_test, n_sig_test


def significance_curve(fpr, tpr, S_tot_pb, B_tot_pb):
    """Un-rescaled S(t), B(t), significance = S/sqrt(S+B) over the ROC grid."""
    S = S_tot_pb * L * tpr
    B = B_tot_pb * L * fpr
    denom = S + B
    sig = np.where(denom > 0, S / np.sqrt(np.maximum(denom, 1e-300)), -1.0)
    return S, B, sig


def find_optimum(sig, raw_bkg, min_raw_bkg=MIN_RAW_BKG):
    """argmax(sig) subject to raw_bkg >= min_raw_bkg.

    Returns (i_star, floor_bound) where floor_bound is True iff the guard
    changed the chosen point relative to the unconstrained argmax (including
    the degenerate case where no point satisfies the floor at all).
    """
    i_unconstrained = int(np.argmax(sig))
    allowed = raw_bkg >= min_raw_bkg
    if not allowed.any():
        return i_unconstrained, True
    i_star = int(np.argmax(np.where(allowed, sig, -1.0)))
    return i_star, (i_star != i_unconstrained)


def analyze_one(net, channel, gen):
    ndir  = net_dir(net, channel)
    files = CHANNEL_FILES[channel]
    roc   = np.load(ndir / f'roc_data_{channel}.npz')
    fpr, tpr, auc = roc[f'fpr_{gen}'], roc[f'tpr_{gen}'], float(roc[f'auc_{gen}'])

    if gen == 'hw':
        S_tot_pb, _sig_n     = kweight_sum_and_n(FRIENDS_DIR / files['sig_hw'])
        B_tot_pb, n_bkg_full = kweight_sum_and_n(FRIENDS_DIR / files['bkg_hw'])
        n_bkg_pop, n_sig_pop = herwig_test_split_counts(ndir, channel, n_bkg_full)
    else:
        S_tot_pb, n_sig_pop = kweight_sum_and_n(FRIENDS_DIR / files['sig_mg5'])
        B_tot_pb, n_bkg_pop = kweight_sum_and_n(FRIENDS_DIR / files['bkg_mg5'])

    S_raw, B_raw, sig = significance_curve(fpr, tpr, S_tot_pb, B_tot_pb)
    raw_bkg = fpr * n_bkg_pop
    raw_sig = tpr * n_sig_pop
    i_star, floor_bound = find_optimum(sig, raw_bkg)

    r = RESCALE[channel]
    S_final = r * S_raw[i_star]
    B_final = r * B_raw[i_star]
    sig_final = float(S_final / np.sqrt(S_final + B_final)) if (S_final + B_final) > 0 else 0.0
    # note: S/(S+B) is unaffected by the rescale factor r (cancels in the ratio)
    purity = float(S_final / (S_final + B_final)) if (S_final + B_final) > 0 else 0.0

    return dict(net=net, channel=channel, gen=gen, auc=auc,
                eff_sig=float(tpr[i_star]), eff_bkg=float(fpr[i_star]),
                S=S_final, B=B_final, purity=purity, significance=sig_final,
                raw_S=int(round(raw_sig[i_star])), raw_B=int(round(raw_bkg[i_star])),
                floor_bound=floor_bound,
                sig_unrescaled=float(sig[i_star]))


# ── Safety net: confirm CHANNEL_FILES matches what each net dir actually uses ──
_FILE_CONST_RE = re.compile(r"""^(BKG_FILE|SIG_FILE|BKG_FILE_MG5|SIG_FILE_MG5)\s*=\s*['"](.+?)['"]""")


def verify_channel_files():
    expected = {
        'BKG_FILE':     lambda c: f"../../friends/summer26/{CHANNEL_FILES[c]['bkg_hw']}",
        'SIG_FILE':     lambda c: f"../../friends/summer26/{CHANNEL_FILES[c]['sig_hw']}",
        'BKG_FILE_MG5': lambda c: f"../../friends/summer26/{CHANNEL_FILES[c]['bkg_mg5']}",
        'SIG_FILE_MG5': lambda c: f"../../friends/summer26/{CHANNEL_FILES[c]['sig_mg5']}",
    }
    for channel in CHANNELS:
        for net in NETS:
            ndir = net_dir(net, channel)
            found = {}
            for fname in ('dataset.py', 'inference.py'):
                fpath = ndir / fname
                if not fpath.exists():
                    continue
                for line in fpath.read_text().splitlines():
                    m = _FILE_CONST_RE.match(line.strip())
                    if m:
                        found[m.group(1)] = m.group(2)
            for key, fn in expected.items():
                if key in found and found[key] != fn(channel):
                    print(f'WARN: {ndir}/{{dataset,inference}}.py {key} = {found[key]!r} '
                          f'!= expected {fn(channel)!r}')


def verify_rescale_invariance():
    """Sanity check: the argmax index is unaffected by RESCALE (checked in main())."""
    for net, channel in (('NNjjBj', 'H'), ('NNkin', 'W')):
        ndir = net_dir(net, channel)
        files = CHANNEL_FILES[channel]
        roc = np.load(ndir / f'roc_data_{channel}.npz')
        fpr, tpr = roc['fpr_hw'], roc['tpr_hw']
        S_tot_pb, _ = kweight_sum_and_n(FRIENDS_DIR / files['sig_hw'])
        B_tot_pb, n_bkg_full = kweight_sum_and_n(FRIENDS_DIR / files['bkg_hw'])
        n_bkg_pop, _ = herwig_test_split_counts(ndir, channel, n_bkg_full)
        _, _, sig = significance_curve(fpr, tpr, S_tot_pb, B_tot_pb)
        raw_bkg = fpr * n_bkg_pop
        i_star, _ = find_optimum(sig, raw_bkg)
        r = RESCALE[channel]
        S_r, B_r, sig_r = significance_curve(fpr, tpr, r * S_tot_pb, r * B_tot_pb)
        i_star_r, _ = find_optimum(sig_r, raw_bkg)
        assert i_star == i_star_r, f'{net} {channel}: rescale changed argmax!'
        assert abs(sig_r[i_star_r] - np.sqrt(r) * sig[i_star]) < 1e-6 * sig_r[i_star_r]
    print('verify_rescale_invariance: OK')


def fmt_thousands(x):
    return f'{x:,.1f}'.replace(',', ' ')


def write_txt_table(rows):
    header = (
        '# summer26 S/B estimate from ROC curves\n'
        '# L = 300 fb^-1 = 300000 pb^-1\n'
        '# S, B, significance = S/sqrt(S+B) at the ROC operating point maximizing\n'
        '# significance, subject to raw B >= %d (statistical-validity floor).\n'
        '# BR rescale applied to S and B alike (see summer26.md):\n'
        '#   H factor = %.6g   Z factor = %.6g   W factor = %.6g\n'
        '# eff_sig/eff_bkg = ROC (tpr, fpr) at the chosen operating point.\n'
        '# floor_bound: 1 if the raw-B>=%d guard changed the unconstrained optimum.\n'
        '#\n' % (MIN_RAW_BKG, RESCALE['H'], RESCALE['Z'], RESCALE['W'], MIN_RAW_BKG)
    )
    cols = f'{"net":<9} {"ch":>2} {"gen":>4} {"AUC":>7} {"eff_sig":>8} {"eff_bkg":>8} ' \
           f'{"S":>10} {"B":>10} {"S/(S+B)":>8} {"S/sqrt(S+B)":>12} {"raw S":>7} {"raw B":>7} {"floor":>5}\n'
    lines = []
    for r in rows:
        lines.append(
            f'{r["net"]:<9} {r["channel"]:>2} {r["gen"]:>4} {r["auc"]:>7.4f} '
            f'{r["eff_sig"]:>8.4f} {r["eff_bkg"]:>8.4f} '
            f'{fmt_thousands(r["S"]):>10} {fmt_thousands(r["B"]):>10} '
            f'{r["purity"]:>8.4f} '
            f'{r["significance"]:>12.3f} {r["raw_S"]:>7d} {r["raw_B"]:>7d} '
            f'{int(r["floor_bound"]):>5d}'
        )
    OUT_TXT.write_text(header + cols + '\n'.join(lines) + '\n')
    print(f'Wrote {OUT_TXT}')


def print_markdown_table(by_key):
    print()
    print('| Net | Ch | AUC (HW) | ε_sig* (HW) | ε_bkg* (HW) | S (HW) | B (HW) | S/(S+B) (HW) | '
          'S/√(S+B) (HW) | raw S (HW) | raw B (HW) | AUC (MG5) | ε_sig* (MG5) | '
          'ε_bkg* (MG5) | S (MG5) | B (MG5) | S/(S+B) (MG5) | S/√(S+B) (MG5) | raw S (MG5) | raw B (MG5) |')
    print('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for channel in CHANNELS:
        for net in NETS:
            hw  = by_key[(net, channel, 'hw')]
            mg5 = by_key[(net, channel, 'mg5')]
            rawB_hw  = f'{hw["raw_B"]}{"†" if hw["floor_bound"] else ""}'
            rawB_mg5 = f'{mg5["raw_B"]}{"†" if mg5["floor_bound"] else ""}'
            print(
                f'| `{net}` | **{channel}** | {hw["auc"]:.4f} | {hw["eff_sig"]:.3f} | '
                f'{hw["eff_bkg"]:.3f} | {fmt_thousands(hw["S"])} | {fmt_thousands(hw["B"])} | '
                f'{hw["purity"]:.3f} | '
                f'**{hw["significance"]:.2f}** | {hw["raw_S"]} | {rawB_hw} | '
                f'{mg5["auc"]:.4f} | {mg5["eff_sig"]:.3f} | {mg5["eff_bkg"]:.3f} | '
                f'{fmt_thousands(mg5["S"])} | {fmt_thousands(mg5["B"])} | '
                f'{mg5["purity"]:.3f} | '
                f'**{mg5["significance"]:.2f}** | {mg5["raw_S"]} | {rawB_mg5} |'
            )
    print()


def main():
    verify_channel_files()
    verify_rescale_invariance()

    rows = []
    by_key = {}
    for channel in CHANNELS:
        for net in NETS:
            for gen in GENS:
                r = analyze_one(net, channel, gen)
                rows.append(r)
                by_key[(net, channel, gen)] = r

    write_txt_table(rows)
    print_markdown_table(by_key)


if __name__ == '__main__':
    main()
