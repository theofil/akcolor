#!/usr/bin/env python3
"""
Leading-jet pull-vector-angle (jetSPVA) distributions inside the S/sqrt(B)-optimal
signal region of NNj, NNjj, NNjjBj (Z channel), and NNjB (Z, W, and their sum
V = W+Z), all from the Herwig test split.

Reuses each net's own inference.py verbatim (load_features/run_nn/JetNN) to
recompute Herwig-test-split NN scores, then reuses sb_estimate.analyze_one()'s
already-established (eff_sig, eff_bkg) operating point (same S/sqrt(B)-optimal
point tabulated in summer26.md) to pick the score threshold on a freshly
computed, full-resolution ROC curve (roc_data_{channel}.npz only stores a
decimated curve, not per-event scores or thresholds, so this cannot be read
back directly).

jetSPVA is (N,3) per-jet in the friend h5 files; the "leading jet" here means
the true pT-leading jet, i.e. jetSPVA[:, leadJetIndex] -- NOT jet slot 0, which
is deliberately randomized 50/50 between the true leading/subleading jet at
fill time for NN training (see README.md "Jet slot randomisation"). This
matches the physics convention used by spva_plots.py's hLeadJetSPVA histograms.

V = W+Z: NNjB is trained and cut independently per channel (its own model, its
own S/sqrt(B)-optimal threshold), so the two channels' absolute-yield
histograms (each already at L=300 fb^-1, BR-rescaled by sb_estimate) are simply
added bin-by-bin -- the combined yield an inclusive W+Z search would see if run
with each channel's own best NNjB working point.

Output:
  figs/summer26/jetSPVA_SR_{NNj,NNjj,NNjjBj}_Z.pdf
  figs/summer26/jetSPVA_SR_NNjB_{Z,W,V}.pdf
"""
import importlib
import pathlib
import pickle
import sys

import numpy as np
import h5py
import torch
from sklearn.metrics import roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).parent
FIGS = HERE.parent / 'figs' / 'summer26'
NETS = ['NNj', 'NNjj', 'NNjjBj']   # Z-only nets, unchanged from before
V_NET = 'NNjB'
V_CHANNELS = ['Z', 'W']            # V = W+Z

sys.path.insert(0, str(HERE))
import sb_estimate  # noqa: E402  (for analyze_one's established operating point)

sys.path.insert(0, '/afs/cern.ch/user/t/theofil/work/akcolor/analysis/Winter25/python/NN')
from style import apply_style, FIG_SIZE  # noqa: E402

SPVA_BINS = np.linspace(0, np.pi, 21)
BIN_CENTERS = (SPVA_BINS[:-1] + SPVA_BINS[1:]) / 2
BIN_WIDTH   = SPVA_BINS[1] - SPVA_BINS[0]


def load_leadjet_spva(h5path):
    """(N,) leading (true pT-highest) jet's |SPVA|, per README's leadJetIndex
    convention. jetSPVA in the friend h5 files is the *signed* pull angle in
    [-pi, pi] (~50/50 split) -- not |theta_s| directly. spva_plots.py's
    hLeadJetSPVA ROOT histograms already take abs() at C++ fill time; the raw
    per-event branch here does not, so it must be taken here too (confirmed
    by a factor-~2 raw-count deficit in the signal-region histograms before
    this fix -- half of events were silently landing outside the [0, pi]
    binning). The -99 no-pull-vector sentinel (jetNC==1, negligible rate for
    the leading jet) maps to +99 under abs() and still falls outside [0, pi],
    so it is still silently and correctly excluded."""
    with h5py.File(h5path, 'r') as f:
        spva     = f['jetSPVA'][:]
        lead_idx = f['leadJetIndex'][:]
    return np.abs(spva[np.arange(len(lead_idx)), lead_idx])


def get_test_scores(net, channel):
    """Reuse net's own inference.py (load_features/run_nn/JetNN) to recompute
    Herwig-test-split scores, y_test, and idx_test (into concat[bkg, sig] rows)."""
    net_dir = sb_estimate.net_dir(net, channel)
    for m in ('dataset', 'model', 'inference'):
        sys.modules.pop(m, None)
    sys.path.insert(0, str(net_dir))
    inf = importlib.import_module('inference')

    feats_bkg = inf.load_features(inf.HERE / inf.BKG_FILE)
    feats_sig = inf.load_features(inf.HERE / inf.SIG_FILE)
    feats_all = [np.concatenate([b, s], axis=0) for b, s in zip(feats_bkg[:-1], feats_sig[:-1])]
    y_all = np.concatenate([np.zeros(len(feats_bkg[0]), dtype=np.float32),
                            np.ones(len(feats_sig[0]), dtype=np.float32)])

    idx_test = np.load(inf.HERE / f'split_indices_{channel}.npz')['test_idx']
    feats_test = [f[idx_test] for f in feats_all]
    y_test = y_all[idx_test]

    with open(inf.HERE / f'scaler_{channel}.pkl', 'rb') as fh:
        scalers = pickle.load(fh)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = inf.JetNN().to(device)
    model.load_state_dict(torch.load(inf.HERE / f'best_model_{channel}.pt', map_location=device))
    model.eval()

    scores_test = inf.run_nn(model, scalers, *feats_test, device)
    return scores_test, y_test, idx_test


def channel_histograms(net, channel):
    """Absolute-yield (events/rad at L=300 fb^-1) VBF/QCD histograms of the
    leading jet's |SPVA| inside net's S/sqrt(B)-optimal signal region for one
    (net, channel). Returns {'VBF': {...}, 'QCD': {...}}, each with
    density/err (arrays over SPVA_BINS), scale (BR-rescaled event yield,
    i.e. analyze_one()'s S or B), and n_pass (raw MC events passing the cut)."""
    files = sb_estimate.CHANNEL_FILES[channel]
    bkg_spva = load_leadjet_spva(sb_estimate.FRIENDS_DIR / files['bkg'])
    sig_spva = load_leadjet_spva(sb_estimate.FRIENDS_DIR / files['sig'])
    spva_all = np.concatenate([bkg_spva, sig_spva])  # same [bkg, sig] row order as every dataset.py

    scores_test, y_test, idx_test = get_test_scores(net, channel)
    spva_test = spva_all[idx_test]

    fpr, tpr, thr = roc_curve(y_test, scores_test)

    r = sb_estimate.analyze_one(net, channel)
    target_fpr, target_tpr = r['eff_bkg'], r['eff_sig']
    i_star = int(np.argmin((fpr - target_fpr) ** 2 + (tpr - target_tpr) ** 2))
    thr_cut = thr[i_star]

    passing = scores_test >= thr_cut
    spva_sig_pass = spva_test[passing & (y_test == 1)]
    spva_bkg_pass = spva_test[passing & (y_test == 0)]
    print(f'{net} {channel}: cut={thr_cut:.4f}  eff_sig={tpr[i_star]:.4f}  eff_bkg={fpr[i_star]:.4f}  '
          f'n_sig_pass={len(spva_sig_pass)}  n_bkg_pass={len(spva_bkg_pass)}  '
          f'S={r["S"]:.1f}  B={r["B"]:.1f}')

    # Absolute normalization at L = 300 fb^-1: each raw-MC passing-event
    # histogram is rescaled so its integral equals analyze_one()'s already
    # BR-rescaled, L=300fb^-1 S (or B) -- the exact same S/B tabulated in
    # summer26.md, not a separate re-derivation.
    out = {}
    for vals, key, scale in ((spva_sig_pass, 'VBF', r['S']), (spva_bkg_pass, 'QCD', r['B'])):
        counts, _ = np.histogram(vals, bins=SPVA_BINS)
        counts = counts.astype(float)
        n_pass = int(counts.sum())
        k = scale / (n_pass * BIN_WIDTH) if n_pass > 0 else 0.0
        out[key] = dict(density=counts * k, err=np.sqrt(counts) * k, scale=scale, n_pass=n_pass)
    return out


def combine_v(d_z, d_w):
    """V = W+Z: bin-by-bin sum of two channels' independently-cut absolute
    yields (errors added in quadrature, scales/raw counts added linearly)."""
    out = {}
    for key in ('VBF', 'QCD'):
        out[key] = dict(
            density=d_z[key]['density'] + d_w[key]['density'],
            err=np.sqrt(d_z[key]['err'] ** 2 + d_w[key]['err'] ** 2),
            scale=d_z[key]['scale'] + d_w[key]['scale'],
            n_pass=d_z[key]['n_pass'] + d_w[key]['n_pass'],
        )
    return out


def save_figure(densities, chan_label, outpath):
    """densities: {'VBF': {...}, 'QCD': {...}} as returned by
    channel_histograms()/combine_v(). Title shows only the expected event
    yields (and, in parenthesis, the raw MC counterparts) -- no efficiencies
    or other text. A black VBF+QCD line shows the total shape passing the cut."""
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    title_parts = []
    for key, color in (('VBF', 'steelblue'), ('QCD', 'tomato')):
        d = densities[key]
        if d['n_pass'] <= 0:
            continue
        title_parts.append(f'{key}: {sb_estimate.fmt_thousands(d["scale"])} ({d["n_pass"]})')
        ax.hist(BIN_CENTERS, bins=SPVA_BINS, weights=d['density'], histtype='step',
                linewidth=3, color=color, label=f'{key} {chan_label} Herwig')
        ax.errorbar(BIN_CENTERS, d['density'], yerr=d['err'], fmt='none',
                    color=color, capsize=3, capthick=1, elinewidth=1)
    if densities['VBF']['n_pass'] > 0 and densities['QCD']['n_pass'] > 0:
        # Black line = VBF + QCD, the total expected yield shape passing the
        # cut(s) -- what an inclusive counting experiment would actually see.
        total = densities['VBF']['density'] + densities['QCD']['density']
        ax.hist(BIN_CENTERS, bins=SPVA_BINS, weights=total, histtype='step',
                linewidth=2, color='black', label='VBF + QCD')
    # log_y: S and B can differ by up to ~20x in absolute scale depending on
    # net/channel (tight cuts push S and B close together, loose cuts don't)
    # -- linear would visually flatten the signal curve to near-zero.
    apply_style(ax, xlabel=r'$|\theta_s|$ [rad]', ylabel=r'Events / rad  ($L=300\,\mathrm{fb}^{-1}$)',
                title='   '.join(title_parts),
                xlim=(0, np.pi), log_y=True, legend_loc='upper right')
    ax.legend(loc='upper right', frameon=False, fontsize=10)  # long labels: shrink to avoid overlap
    plt.tight_layout()
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved', outpath)
    plt.close(fig)


def main():
    FIGS.mkdir(parents=True, exist_ok=True)

    for net in NETS:
        densities = channel_histograms(net, 'Z')
        save_figure(densities, 'Z', str(FIGS / f'jetSPVA_SR_{net}_Z.pdf'))

    per_channel = {}
    for channel in V_CHANNELS:
        per_channel[channel] = channel_histograms(V_NET, channel)
        save_figure(per_channel[channel], channel, str(FIGS / f'jetSPVA_SR_{V_NET}_{channel}.pdf'))

    v_densities = combine_v(per_channel['Z'], per_channel['W'])
    save_figure(v_densities, 'V=W+Z', str(FIGS / f'jetSPVA_SR_{V_NET}_V.pdf'))


if __name__ == '__main__':
    main()
