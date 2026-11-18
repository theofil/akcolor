#!/usr/bin/env python3
"""Training generator x test generator: is the Herwig -> MG5+Pythia loss a transfer loss? (lxplus)

    python3 nov26/train_gen_matrix.py --campaign fall26                 # NNjPF, 400 replicas
    python3 nov26/train_gen_matrix.py --campaign fall26 --net ParTj --replicas 200

Needs both h5 score columns: <net> (trained on Herwig) and <net>_py (trained on MG5+Pythia,
train_fold.py --train-gen pythia; save_scores.py --train-gen pythia). Both are out-of-fold on
all four samples, and model k of either scored the same inference sub-folds.

Prints, for the categories of roc_scores.py (gate jets, gate == 0, VBF presel.):
  1. the 2 x 2 matrix (trained on) x (tested on) of AUC and eps_S(1 %), with bootstrap sigma;
  2. the transfer loss on each test generator: the net trained on that generator minus the net
     trained on the other one, on the same test jets (paired: same events, same replicas), with
     the per-fold-model mean +- standard error next to it (includes the training randomness);
  3. the generator gap at matched training (trained and tested on Herwig minus trained and
     tested on MG5+Pythia; independent samples, sigmas added in quadrature).
If the in-generator nets match each other, the Herwig-trained net's loss on MG5+Pythia is a
transfer loss; if the transfer loss is zero, the MG5+Pythia jets are intrinsically harder to
separate. The metrics are event-weighted when the campaign has event weights (fall26).
Figure figs/<campaign>/roc_train_gen_<net>.{pdf,png}: rejection vs signal efficiency, colour =
test generator, solid = trained on the test generator, dashed = trained on the other one; lower
panels: rejection ratio transfer / in-generator.
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import roc_scores as R  # noqa: E402
from bootstrap_scores import event_index, trapezoid  # noqa: E402

GEN_LABEL = {'herwig': 'Herwig', 'pythia': 'MG5+Pythia'}
C_TEST = {'herwig': '#2a78d6', 'pythia': '#eb6834'}                    # categorical slots 1-2
OTHER = {'herwig': 'pythia', 'pythia': 'herwig'}


def prepare(sig, bkg, ms, mb, ev_s, ev_b, weighted):
    """Selected jets sorted once by descending score, with the last index of each tie group."""
    s = np.r_[sig['score'][ms], bkg['score'][mb]]
    y = np.r_[np.ones(ms.sum(), bool), np.zeros(mb.sum(), bool)]
    kw = np.r_[sig['kw'][ms], bkg['kw'][mb]].astype(np.float64) if weighted else np.ones(len(s))
    ev = np.r_[ev_s[ms], ev_b[mb]]
    o = np.argsort(-s, kind='stable')
    last = np.r_[np.nonzero(np.diff(s[o]))[0], len(s) - 1]
    return y[o], kw[o], ev[o], last


def metrics(p, m_s=None, m_b=None):
    """(eps_S at eps_B = 1 %, AUC) of the weighted ROC curve, optionally with event multiplicities."""
    y, kw, ev, last = p
    w = kw if m_s is None else kw * np.where(y, m_s[ev], m_b[ev])
    cs = np.cumsum(np.where(y, w, 0.))[last]
    cb = np.cumsum(np.where(y, 0., w))[last]
    tpr, fpr = np.r_[0., cs / cs[-1]], np.r_[0., cb / cb[-1]]
    return np.array([np.interp(R.EPS_B, fpr, tpr), trapezoid(tpr, fpr)])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--campaign', default='nov26', choices=sorted(C.CAMPAIGNS))
    ap.add_argument('--net', default='NNjPF', choices=R.NETS)
    ap.add_argument('--replicas', type=int, default=400)
    ap.add_argument('--seed', type=int, default=20261118)
    args = ap.parse_args()
    C.set_campaign(args.campaign)
    weighted = C.EVENT_WEIGHTS
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    # data[train, test] = (sig jets, bkg jets) scored by the net trained on `train`
    data = {}
    for train in C.GEN_SAMPLES:
        column = C.net_name(args.net, train)
        for test, (sig_name, bkg_name) in C.GEN_SAMPLES.items():
            sig, bkg = R.jets(sig_name, column), R.jets(bkg_name, column)
            if sig is None or bkg is None:
                sys.exit(f'{C.CAMPAIGN}: no {column} column in {sig_name}/{bkg_name}.h5 (run save_scores.py)')
            data[train, test] = (sig, bkg)

    # bootstrap: per test generator, one set of replicas shared by both training generators
    cats = list(R.PAIRED)
    central, sigma, reps_all, per_model = {}, {}, {}, {}
    for test, (sig_name, bkg_name) in C.GEN_SAMPLES.items():
        (ev_s, n_s), (ev_b, n_b) = event_index(sig_name), event_index(bkg_name)
        prep = {}
        for train in C.GEN_SAMPLES:
            sig, bkg = data[train, test]
            assert len(ev_s) == len(sig['score']) and len(ev_b) == len(bkg['score'])
            for cat in cats:
                sel = R.CATS[cat]
                prep[train, cat] = prepare(sig, bkg, sel(sig), sel(bkg), ev_s, ev_b, weighted)
                per_model[train, test, cat] = np.array([
                    metrics(prepare(sig, bkg, sel(sig) & (sig['model'] == k), sel(bkg) & (bkg['model'] == k),
                                    ev_s, ev_b, weighted)) for k in range(C.N_FOLDS)])
        reps = {key: [] for key in prep}
        for _ in range(args.replicas):
            m_s = rng.poisson(1., n_s).astype(np.float64)
            m_b = rng.poisson(1., n_b).astype(np.float64)
            for key, p in prep.items():
                reps[key].append(metrics(p, m_s, m_b))
        for (train, cat), p in prep.items():
            central[train, test, cat] = metrics(p)
            reps_all[train, test, cat] = np.array(reps[train, cat])
            sigma[train, test, cat] = reps_all[train, test, cat].std(0, ddof=1)
    print(f'({time.time() - t0:.0f} s, {args.replicas} bootstrap replicas)\n')

    wtxt = 'event-weighted' if weighted else 'unweighted'
    names = ('eps_S(1%)', 'AUC')                                         # order of metrics()
    print(f'{C.CAMPAIGN} {args.net}: {wtxt} figures of merit, trained on (rows) x tested on (columns)')
    print(f'  {"category":11s} {"metric":9s} {"trained on":11s} '
          + ' '.join(f'{"test " + GEN_LABEL[t]:>19s}' for t in C.GEN_SAMPLES))
    for cat in cats:
        for i, name in enumerate(names):
            for train in C.GEN_SAMPLES:
                print(f'  {cat:11s} {name:9s} {GEN_LABEL[train]:11s} ' + ' '.join(
                    f'{central[train, t, cat][i]:10.4f} ± {sigma[train, t, cat][i]:.4f}' for t in C.GEN_SAMPLES))

    print(f'\nTransfer loss on each test generator: (trained on it) - (trained on the other), same test jets')
    print(f'  {"category":11s} {"metric":9s} {"test":11s} {"bootstrap":>19s}   per fold model (mean ± s.e.)')
    for cat in cats:
        for i, name in enumerate(names):
            for test in C.GEN_SAMPLES:
                a, b = (test, test, cat), (OTHER[test], test, cat)
                d = central[a][i] - central[b][i]
                s = (reps_all[a][:, i] - reps_all[b][:, i]).std(ddof=1)
                pm = per_model[test, test, cat][:, i] - per_model[OTHER[test], test, cat][:, i]
                print(f'  {cat:11s} {name:9s} {GEN_LABEL[test]:11s} {d:+10.4f} ± {s:.4f}   '
                      + ' '.join(f'{x:+.4f}' for x in pm)
                      + f'   ({pm.mean():+.4f} ± {pm.std(ddof=1) / np.sqrt(len(pm)):.4f})')

    print(f'\nGenerator gap: (Herwig-trained on Herwig) - (MG5+Pythia-trained on MG5+Pythia), and for '
          f'reference the Herwig-trained net alone')
    for cat in cats:
        for i, name in enumerate(names):
            hh, pp, hp = ('herwig', 'herwig', cat), ('pythia', 'pythia', cat), ('herwig', 'pythia', cat)
            g = central[hh][i] - central[pp][i]
            sg = np.hypot(sigma[hh][i], sigma[pp][i])
            r = central[hh][i] - central[hp][i]
            sr = np.hypot(sigma[hh][i], sigma[hp][i])
            print(f'  {cat:11s} {name:9s} matched training {g:+.4f} ± {sg:.4f}    '
                  f'Herwig-trained H - P {r:+.4f} ± {sr:.4f}')

    # ── figure ───────────────────────────────────────────────────────────────
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, len(cats), sharex='col', squeeze=False, facecolor=R.SURFACE,
                             figsize=(5.2 * len(cats), 6.8), gridspec_kw={'height_ratios': [3, 1.4]})
    for col, cat in enumerate(cats):
        ax, axr = axes[0, col], axes[1, col]
        curves = {}
        for test in C.GEN_SAMPLES:
            for train in (test, OTHER[test]):
                sig, bkg = data[train, test]
                ms, mb = R.CATS[cat](sig), R.CATS[cat](bkg)
                tpr, fpr = R.roc_points(sig, bkg, ms, mb, weighted)
                curves[train, test] = (tpr, fpr)
                e, a = central[train, test, cat]
                ax.plot(tpr, 1 / fpr, color=C_TEST[test], lw=1.8, ls='-' if train == test else '--',
                        label=f'test {GEN_LABEL[test]}, trained {GEN_LABEL[train]}: AUC {a:.3f}, ε_S(1%) {e:.3f}')
            (t1, f1), (t2, f2) = curves[test, test], curves[OTHER[test], test]
            grid = np.linspace(max(t1[0], t2[0], 0.05), 1, 400)
            axr.plot(grid, R.rejection_at(t2, f2, grid) / R.rejection_at(t1, f1, grid), color=C_TEST[test],
                     lw=1.5, label=f'test {GEN_LABEL[test]}')
        ax.axhline(1 / R.EPS_B, color=R.INK2, lw=0.9, ls=':')
        ax.text(0.98, 1 / R.EPS_B * 1.12, 'ε_B = 1 %', color=R.INK2, fontsize=8, ha='right', va='bottom',
                transform=ax.get_yaxis_transform())
        ax.set_yscale('log')
        ax.set_title(R.PANEL_TITLE[cat], fontsize=10, color=R.INK)
        ax.legend(frameon=False, fontsize=7, loc='upper right')
        axr.axhline(1, color=R.INK2, lw=0.8)
        axr.legend(frameon=False, fontsize=8, loc='best')
        if col == 0:
            ax.set_ylabel('background rejection 1 / ε_B (QCD Z+jj jets)', color=R.INK)
            axr.set_ylabel('rejection\ntransfer / in-generator', color=R.INK)
        axr.set_xlabel('signal efficiency ε_S (VBF Z jets)', color=R.INK)
    for a in axes.flat:
        a.set_xlim(0, 1)
        a.grid(color=R.RULE, lw=0.4, alpha=0.7, which='both')
        a.tick_params(colors=R.INK2)
        for sp in a.spines.values():
            sp.set_color(R.RULE)
    fig.suptitle(f'{args.net} trained on Herwig or on MG5+Pythia, out-of-fold, {C.CAMPAIGN}, {wtxt}',
                 fontsize=11, color=R.INK)
    fig.tight_layout()
    os.makedirs(C.FIGS, exist_ok=True)
    out = os.path.join(C.FIGS, f'roc_train_gen_{args.net}')
    for ext in ('pdf', 'png'):
        fig.savefig(f'{out}.{ext}', dpi=150, facecolor=R.SURFACE)
    print(f'\nwrote {out}.pdf and .png')


if __name__ == '__main__':
    main()
