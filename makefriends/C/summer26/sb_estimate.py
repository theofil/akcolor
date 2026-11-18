#!/usr/bin/env python3
"""
Estimate the best achievable signal/background for every NN* variant directly
from its saved ROC curve (roc_data_{H,W,Z}.npz in each net directory), at
L = 300 fb^-1. Herwig only -- MadGraph5+Pythia generator-transfer numbers are
not computed here (see summer26.md "Signal/background estimation from ROC
curves" for why).

For each (net, channel) combination, scan the net's own Herwig ROC curve
(tpr/fpr) and report two points, both selected using the background-dominated
figure of merit significance = S/sqrt(B):

  - the *initial*, unoptimized point (suffix `0`): S, B at the trivial
    eff_sig=eff_bkg=1 point (i.e. the full preselected sample, no NN cut at
    all). Depends only on the channel's total kWeight yields, so it is the
    same for every net within a channel.
  - the *final*, optimized point: the ROC operating point that maximizes
    S/sqrt(B), subject to a raw-MC-background statistical floor.

S and B are formed from the sample's post-selection sum(kWeight) (already in
pb, see kWeight formula in summer26.md) times L times the ROC efficiency at
that point.

`analyze_one(net, channel)` is also imported directly by `plots.py` to place
a star marker at the final operating point on each net's Herwig ROC curve.

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

Because the same factor multiplies S and B within a channel, it cancels in
S/B and does not shift the argmax of S/sqrt(B) over the ROC curve -- the
scan runs once on un-rescaled yields, and the rescale is applied only to the
winning point's final S, B (see verify_rescale_invariance()).

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

L               = 300_000.0   # pb^-1 = 300 fb^-1
MIN_RAW_BKG     = 10          # statistical-validity floor (repo precedent: old SR_optimization.txt)

NETS     = ['NNkin', 'NNj', 'NNjB', 'NNjj', 'NNjjB', 'NNjjBj']
CHANNELS = ['H', 'W', 'Z']

# Feedforward, no-constituents nets (see summer26.md "Z_only_table"/
# "W_only_table") -- kept out of NETS/CHANNELS above since net_dir()/
# analyze_one() are called generically over the full NETS x CHANNELS grid
# elsewhere in this file, and none of these have an _H directory (H not
# requested for this ablation family). Maps net -> list of channels it has.
EXTRA_NETS = {'Kj': ['Z', 'W'], 'KjB': ['Z', 'W'], 'KjjBj': ['Z', 'W']}

CHANNEL_FILES = {
    'H': dict(bkg='QCDHjj_herwig.h5', sig='VBFH_herwig.h5'),
    'W': dict(bkg='QCDWjj_herwig.h5', sig='VBFW_herwig.h5'),
    'Z': dict(bkg='QCDZjj_herwig.h5', sig='VBFZ_herwig.h5'),
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
    """Un-rescaled S(t), B(t), and significance = S/sqrt(B) over the ROC grid."""
    S = S_tot_pb * L * tpr
    B = B_tot_pb * L * fpr
    sig = np.where(B > 0, S / np.sqrt(np.maximum(B, 1e-300)), -1.0)
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


def analyze_one(net, channel, min_raw_bkg=MIN_RAW_BKG):
    """Initial (no-cut) and final operating points for one (net, channel),
    both BR-rescaled. The final operating point is the ROC point maximizing
    S/sqrt(B), subject to a raw-B >= min_raw_bkg statistical-validity floor.
    Used both by main() below and directly by plots.py to place the final
    point as a star on the Herwig ROC curve."""
    ndir  = net_dir(net, channel)
    files = CHANNEL_FILES[channel]
    roc   = np.load(ndir / f'roc_data_{channel}.npz')
    fpr, tpr, auc = roc['fpr_hw'], roc['tpr_hw'], float(roc['auc_hw'])

    S_tot_pb, _sig_n     = kweight_sum_and_n(FRIENDS_DIR / files['sig'])
    B_tot_pb, n_bkg_full = kweight_sum_and_n(FRIENDS_DIR / files['bkg'])
    n_bkg_pop, n_sig_pop = herwig_test_split_counts(ndir, channel, n_bkg_full)

    S_raw, B_raw, sig_sqrtB = significance_curve(fpr, tpr, S_tot_pb, B_tot_pb)
    raw_bkg = fpr * n_bkg_pop
    raw_sig = tpr * n_sig_pop
    r = RESCALE[channel]

    i_star, floor_bound = find_optimum(sig_sqrtB, raw_bkg, min_raw_bkg=min_raw_bkg)
    S_final    = r * S_raw[i_star]
    B_final    = r * B_raw[i_star]
    sig_final  = float(S_final / np.sqrt(B_final)) if B_final > 0 else 0.0
    purity     = float(S_final / (S_final + B_final)) if (S_final + B_final) > 0 else 0.0

    # initial: unoptimized baseline -- the full preselected sample, no NN cut
    # at all (eff_sig=eff_bkg=1). Same for every net within a channel.
    S0 = r * S_tot_pb * L
    B0 = r * B_tot_pb * L
    sig0 = float(S0 / np.sqrt(B0)) if B0 > 0 else 0.0

    return dict(net=net, channel=channel, auc=auc,
                S0=S0, B0=B0, raw_S0=n_sig_pop, raw_B0=n_bkg_pop,
                significance0=sig0,
                eff_sig=float(tpr[i_star]), eff_bkg=float(fpr[i_star]),
                S=S_final, B=B_final,
                raw_S=int(round(raw_sig[i_star])), raw_B=int(round(raw_bkg[i_star])),
                significance=sig_final, purity=purity,
                floor_bound=floor_bound)


# ── Safety net: confirm CHANNEL_FILES matches what each net dir actually uses ──
_FILE_CONST_RE = re.compile(r"""^(BKG_FILE|SIG_FILE)\s*=\s*['"](.+?)['"]""")


def verify_channel_files():
    expected = {
        'BKG_FILE': lambda c: f"../../friends/summer26/{CHANNEL_FILES[c]['bkg']}",
        'SIG_FILE': lambda c: f"../../friends/summer26/{CHANNEL_FILES[c]['sig']}",
    }
    combos = [(net, channel) for channel in CHANNELS for net in NETS] + \
             [(net, channel) for net, channels in EXTRA_NETS.items() for channel in channels]
    for net, channel in combos:
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
    """Sanity check: the argmax index of S/sqrt(B) is unaffected by RESCALE."""
    for net, channel in (('NNjjBj', 'H'), ('NNkin', 'W')):
        ndir = net_dir(net, channel)
        files = CHANNEL_FILES[channel]
        roc = np.load(ndir / f'roc_data_{channel}.npz')
        fpr, tpr = roc['fpr_hw'], roc['tpr_hw']
        S_tot_pb, _ = kweight_sum_and_n(FRIENDS_DIR / files['sig'])
        B_tot_pb, n_bkg_full = kweight_sum_and_n(FRIENDS_DIR / files['bkg'])
        n_bkg_pop, _ = herwig_test_split_counts(ndir, channel, n_bkg_full)
        _, _, sig = significance_curve(fpr, tpr, S_tot_pb, B_tot_pb)
        raw_bkg = fpr * n_bkg_pop
        r = RESCALE[channel]
        _, _, sig_r = significance_curve(fpr, tpr, r * S_tot_pb, r * B_tot_pb)

        i_star, _ = find_optimum(sig, raw_bkg)
        i_star_r, _ = find_optimum(sig_r, raw_bkg)
        assert i_star == i_star_r, f'{net} {channel}: rescale changed argmax!'
        assert abs(sig_r[i_star_r] - np.sqrt(r) * sig[i_star]) < 1e-6 * sig_r[i_star_r]
    print('verify_rescale_invariance: OK')


def fmt_thousands(x):
    return f'{x:,.1f}'.replace(',', ' ')


def write_txt_table(rows):
    header = (
        '# summer26 S/B estimate from ROC curves (Herwig only)\n'
        '# L = 300 fb^-1 = 300000 pb^-1\n'
        '# significance = S/sqrt(B), the only optimization objective.\n'
        '# _0: unoptimized baseline -- S, B at the full preselected sample, no NN\n'
        '#   cut at all (eff_sig=eff_bkg=1). Same for every net within a channel.\n'
        '# (no suffix): ROC operating point maximizing S/sqrt(B), subject to a\n'
        '#   raw B >= %d statistical-validity floor.\n'
        '# BR rescale applied to S and B alike (see summer26.md):\n'
        '#   H factor = %.6g   Z factor = %.6g   W factor = %.6g\n'
        '# eff_sig/eff_bkg = ROC (tpr, fpr) at the chosen operating point.\n'
        '# floor: 1 if the raw-B>=%d guard changed the unconstrained optimum.\n'
        '#\n' % (MIN_RAW_BKG, RESCALE['H'], RESCALE['Z'], RESCALE['W'], MIN_RAW_BKG)
    )
    cols = (
        f'{"net":<9} {"ch":>2} {"AUC":>7} '
        f'{"S_0":>10} {"B_0":>10} {"rawS_0":>7} {"rawB_0":>7} '
        f'{"S/sqrt(B)_0":>11} '
        f'{"eff_sig":>8} {"eff_bkg":>8} {"S":>10} {"B":>10} {"rawS":>7} {"rawB":>7} '
        f'{"S/sqrt(B)":>10} {"floor":>5}\n'
    )
    lines = []
    for r in rows:
        lines.append(
            f'{r["net"]:<9} {r["channel"]:>2} {r["auc"]:>7.4f} '
            f'{fmt_thousands(r["S0"]):>10} {fmt_thousands(r["B0"]):>10} '
            f'{r["raw_S0"]:>7d} {r["raw_B0"]:>7d} '
            f'{r["significance0"]:>11.3f} '
            f'{r["eff_sig"]:>8.4f} {r["eff_bkg"]:>8.4f} '
            f'{fmt_thousands(r["S"]):>10} {fmt_thousands(r["B"]):>10} '
            f'{r["raw_S"]:>7d} {r["raw_B"]:>7d} '
            f'{r["significance"]:>10.3f} {int(r["floor_bound"]):>5d}'
        )
    OUT_TXT.write_text(header + cols + '\n'.join(lines) + '\n')
    print(f'Wrote {OUT_TXT}')


def print_markdown_table(by_key):
    print()
    print('Subscript legend: `_0` = unoptimized baseline (full preselected sample, '
          'no NN cut, same for every net within a channel); unsuffixed columns = '
          'the ROC operating point maximizing S/√B (Herwig only) -- the only '
          'optimization objective. '
          '`†` on `raw B` marks a row where the `raw B ≥ 10` floor changed the '
          'unconstrained optimum.')
    print()
    print('| Net | Ch | S/√B₀ | S₀ | B₀ | raw S₀ | raw B₀ | S/√B | '
          'S | B | raw S | raw B |')
    print('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for channel in CHANNELS:
        for net in NETS:
            r = by_key[(net, channel)]
            rawB = f'{r["raw_B"]}{"†" if r["floor_bound"] else ""}'
            print(
                f'| `{net}` | **{channel}** | {r["significance0"]:.2f} | '
                f'{fmt_thousands(r["S0"])} | {fmt_thousands(r["B0"])} | '
                f'{r["raw_S0"]} | {r["raw_B0"]} | '
                f'**{r["significance"]:.2f}** | '
                f'{fmt_thousands(r["S"])} | '
                f'{fmt_thousands(r["B"])} | {r["raw_S"]} | {rawB} |'
            )
    print()


Z_ONLY_ROW_ORDER = ['Kj', 'NNj', 'KjB', 'NNjB', 'KjjBj', 'NNjjBj']
W_ONLY_ROW_ORDER = ['Kj', 'NNj', 'KjB', 'NNjB', 'KjjBj', 'NNjjBj']


def _print_only_table(by_key, channel, row_order, description):
    """Shared body for print_z_only_table/print_w_only_table: the `_0`
    unoptimized baseline is identical for every row within a channel (see
    analyze_one docstring) so it is stated once in prose instead of repeated
    as six columns."""
    r0 = by_key[(row_order[0], channel)]
    print()
    print(f'Unoptimized baseline (`_0`: full preselected sample, no NN cut, '
          f'`eff_sig=eff_bkg=1` -- identical for every net in the {channel} '
          f'channel): S/√B₀ = {r0["significance0"]:.2f}, S₀ = {fmt_thousands(r0["S0"])}, '
          f'B₀ = {fmt_thousands(r0["B0"])}, raw S₀ = {r0["raw_S0"]}, '
          f'raw B₀ = {r0["raw_B0"]}.')
    print()
    print(description)
    print()
    print('| Net | S/√B | S/(S+B) | eff_sig | eff_bkg | S | B | raw S | raw B |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for net in row_order:
        r = by_key[(net, channel)]
        rawB = f'{r["raw_B"]}{"†" if r["floor_bound"] else ""}'
        print(
            f'| `{net}` | **{r["significance"]:.2f}** | {r["purity"]:.3f} | '
            f'{r["eff_sig"]:.4f} | {r["eff_bkg"]:.4f} | '
            f'{fmt_thousands(r["S"])} | {fmt_thousands(r["B"])} | '
            f'{r["raw_S"]} | {rawB} |'
        )
    print()


def print_z_only_table(by_key):
    """Z_only_table: Kj/KjB/KjjBj (feedforward, no constituents) each
    inserted directly above their DeepSets counterpart (NNj/NNjB/NNjjBj),
    isolating the value added by jet constituents. Same `S/√B` objective as
    print_markdown_table's main table."""
    description = (
        '`Kj`/`KjB`/`KjjBj` (feedforward on jet four-vectors + boson '
        'scalars, no constituents) are each inserted directly above their '
        'DeepSets counterpart (`NNj`/`NNjB`/`NNjjBj`) so the added value of '
        'the constituents can be read off directly. Columns below are at '
        'the ROC operating point maximizing `S/√B` (the only optimization '
        'objective); `S/(S+B)` (purity) and `eff_sig`/`eff_bkg` are shown '
        'for context only. '
        '`raw S`/`raw B` are the raw MC event counts actually passing the '
        'cut in the signal region (Herwig test split). `†` on `raw B` marks '
        f'a row where the `raw B ≥ {MIN_RAW_BKG}` statistical-validity '
        'floor changed the unconstrained optimum.'
    )
    _print_only_table(by_key, 'Z', Z_ONLY_ROW_ORDER, description)


def print_w_only_table(by_key):
    """W_only_table: the W-channel analog of print_z_only_table. Kj/KjB/
    KjjBj (feedforward, no constituents) are each inserted directly above
    their DeepSets counterpart (NNj/NNjB/NNjjBj)."""
    description = (
        '`Kj`/`KjB`/`KjjBj` (feedforward on jet four-vectors + boson '
        'scalars, no constituents) are each inserted directly above their '
        'DeepSets counterpart (`NNj`/`NNjB`/`NNjjBj`) so the added value of '
        'the constituents can be read off directly -- the W-channel analog '
        'of the `Z_only_table` above. Columns below are at the ROC '
        'operating point maximizing `S/√B` (the only optimization '
        'objective); `S/(S+B)` (purity) and `eff_sig`/`eff_bkg` are shown '
        'for context only. `raw S`/`raw B` are '
        'the raw MC event counts actually passing the cut in the signal '
        f'region (Herwig test split). `†` on `raw B` marks a row where the '
        f'`raw B ≥ {MIN_RAW_BKG}` statistical-validity floor changed the '
        'unconstrained optimum.'
    )
    _print_only_table(by_key, 'W', W_ONLY_ROW_ORDER, description)


def main():
    verify_channel_files()
    verify_rescale_invariance()

    rows = []
    by_key = {}
    for channel in CHANNELS:
        for net in NETS:
            r = analyze_one(net, channel)
            rows.append(r)
            by_key[(net, channel)] = r

    for net, channels in EXTRA_NETS.items():
        for channel in channels:
            r = analyze_one(net, channel)
            rows.append(r)
            by_key[(net, channel)] = r

    write_txt_table(rows)
    print_markdown_table(by_key)
    print_z_only_table(by_key)
    print_w_only_table(by_key)


if __name__ == '__main__':
    main()
