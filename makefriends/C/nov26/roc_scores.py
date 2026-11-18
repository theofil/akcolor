#!/usr/bin/env python3
"""Out-of-fold performance of the jet nets, from the scores in the h5 files (lxplus; no model).

    python3 nov26/roc_scores.py                     # nov26: every net with scores -> figs/nov26/roc_nets.{pdf,png}
    python3 nov26/roc_scores.py --campaign fall26   # fall26 (event-weighted)     -> figs/fall26/roc_nets.{pdf,png}
    python3 nov26/roc_scores.py --from-npz          # read <campaign>/<net>/scores_fold*.npz instead of the h5 columns

Two figures of merit, for VBF jets (signal) vs QCD jets (background):
  AUC            area under the ROC curve
  eps_S(1%)      signal efficiency at 1 % background efficiency (the ROC curve interpolated
                 at eps_B = 0.01), i.e. the VBF jets kept when 99 % of the QCD jets are rejected
For each net and generator (Herwig: the training generator; MG5+Pythia: transfer test) both are
printed for
  gate jets      roles A + B, both stored jets of gate events
  role A, role B the 1st (slot 0) and 2nd (slot 1) jet of gate events, separately; the slots are
                 randomised and both are out-of-fold, so the two must agree within statistics
  gate == 0      role C
  VBF presel.    gate jets of events with mjj > 400 GeV and |dYjj| > 2.5 (summer26 selection)
with the spread of the five fold models (Herwig gate jets, each model on its own inference
events). With both nets it also prints the paired difference ParTj - NNjPF per fold model: the
nets scored the same jets, so the sample fluctuation cancels; mean +- standard error of the five.
nov26 metrics are unweighted, as in summer26. fall26 metrics and ROC curves are weighted with the
event weights kWeight; the unweighted values and the effective numbers of jets are printed too.
The figure shows, for the gate jets, the gate == 0 jets and the VBF-preselected gate jets, the
background rejection 1/eps_B vs the signal efficiency down to 10 surviving background jets, with
a guide at eps_B = 1 %; when both nets have scores, lower panels show the ratio of their
rejections at equal signal efficiency.
"""
import argparse
import os
import sys

import h5py
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import save_scores  # noqa: E402

NETS = ('NNjPF', 'ParTj')
GENERATORS = {'Herwig': C.TRAIN_SAMPLES, 'MG5+Pythia': C.TRANSFER_SAMPLES}
C_NET = {'NNjPF': '#2a78d6', 'ParTj': '#eb6834'}                    # categorical slots 1-2
LS = {'Herwig': '-', 'MG5+Pythia': '--'}
INK, INK2, RULE, SURFACE = '#0b0b0b', '#52514e', '#c3c2b7', '#fcfcfb'
CATS = {'gate jets': lambda j: j['role'] != C.ROLE_C, 'role A': lambda j: j['role'] == C.ROLE_A,
        'role B': lambda j: j['role'] == C.ROLE_B,
        'gate == 0': lambda j: j['role'] == C.ROLE_C,
        'VBF presel.': lambda j: (j['role'] != C.ROLE_C) & j['vbf']}
PAIRED = ('gate jets', 'gate == 0', 'VBF presel.')
PANEL_TITLE = {'gate jets': 'gate jets (roles A + B)', 'gate == 0': 'gate == 0 jets (role C)',
               'VBF presel.': 'gate jets, VBF-preselected\n(mjj > 400 GeV, |ΔY_jj| > 2.5)'}
EPS_B = 0.01                                                         # eps_S is quoted at eps_B = 1 %


def rejection_at(tpr, fpr, grid):
    """Background rejection 1/fpr of a ROC curve, interpolated at the signal efficiencies grid."""
    t, i = np.unique(tpr, return_index=True)
    return 1 / np.interp(grid, t, fpr[i])


def jets(sample, net, from_npz=False):
    """Per-jet score, event weight, role, VBF-preselection flag and model (fold) index of one sample."""
    with h5py.File(C.h5_path(sample), 'r') as f:
        n_jets, gate, kw = f['nJets'][:], f['gate'][:], f['kWeight'][:]
        vbf = gate & (f['mjj'][:] > 400) & (np.abs(f['dYjj'][:]) > 2.5)
        if from_npz:
            folds_dir = os.path.join(C.OUT, net)
            if not all(os.path.exists(os.path.join(folds_dir, f'scores_fold{k}.npz')) for k in range(C.N_FOLDS)):
                return None
            sub = C.assign_subfolds(sample, gate, n_jets)
            score = save_scores.merge(sample, n_jets, sub, folds_dir)
        elif net in f:
            score, sub = f[net][:], f['subfold'][:]
        else:
            return None
    model = sub // 2                                                  # model k scored sub-folds 2k, 2k+1
    out = {k: [] for k in ('score', 'kw', 'role', 'vbf', 'model')}
    for s in (0, 1):
        m = s < n_jets
        out['score'].append(score[m, s])
        out['kw'].append(kw[m])
        out['role'].append(np.where(gate[m], C.ROLE_A if s == 0 else C.ROLE_B, C.ROLE_C))
        out['vbf'].append(vbf[m])
        out['model'].append(model[m])
    return {k: np.concatenate(v) for k, v in out.items()}


def labelled(sig, bkg, ms, mb, weighted):
    """Labels, scores and (if weighted) event weights of the selected signal and background jets."""
    y = np.r_[np.ones(ms.sum()), np.zeros(mb.sum())]
    s = np.r_[sig['score'][ms], bkg['score'][mb]]
    w = np.r_[sig['kw'][ms], bkg['kw'][mb]].astype(np.float64) if weighted else None
    return y, s, w


def auc_of(sig, bkg, ms, mb, weighted):
    y, s, w = labelled(sig, bkg, ms, mb, weighted)
    return roc_auc_score(y, s, sample_weight=w)


def eff_at(sig, bkg, ms, mb, weighted, eps_b=EPS_B):
    """Signal efficiency at background efficiency eps_b (ROC curve interpolated linearly)."""
    y, s, w = labelled(sig, bkg, ms, mb, weighted)
    fpr, tpr, _ = roc_curve(y, s, sample_weight=w)
    return float(np.interp(eps_b, fpr, tpr))


METRICS = {'AUC': auc_of, 'eps_S(1%)': eff_at}


def roc_points(sig, bkg, ms, mb, weighted):
    """(tpr, fpr) of the ROC curve, down to 10 surviving background jets (counted, not weighted)."""
    y, s, w = labelled(sig, bkg, ms, mb, weighted)
    fpr, tpr, thr = roc_curve(y, s, sample_weight=w)
    sb = np.sort(bkg['score'][mb])
    keep = len(sb) - np.searchsorted(sb, thr, side='left') >= 10
    return tpr[keep], fpr[keep]


def per_model(sig, bkg, cat, weighted, metric=auc_of):
    """Figure of merit of category `cat` for each fold model, on that model's own inference jets."""
    return [metric(sig, bkg, CATS[cat](sig) & (sig['model'] == k), CATS[cat](bkg) & (bkg['model'] == k),
                   weighted) for k in range(C.N_FOLDS)]


def print_table(data, weighted, metric=auc_of, title=None):
    if title:
        print(title)
    print(f'{"net":6s} {"generator":10s} ' + ' '.join(f'{c:>12s}' for c in CATS) + '   fold-model spread')
    for (net, gen), (sig, bkg) in data.items():
        row = [metric(sig, bkg, sel(sig), sel(bkg), weighted) for sel in CATS.values()]
        spread = ''
        if gen == 'Herwig':
            per = per_model(sig, bkg, 'gate jets', weighted, metric)
            spread = f'   {" ".join(f"{a:.4f}" for a in per)}  (rms {np.std(per):.4f})'
        print(f'{net:6s} {gen:10s} ' + ' '.join(f'{a:12.4f}' for a in row) + spread)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--campaign', default='nov26', choices=sorted(C.CAMPAIGNS))
    ap.add_argument('--from-npz', action='store_true', help='read the fold npz files, not the h5 columns')
    args = ap.parse_args()
    C.set_campaign(args.campaign)
    weighted = C.EVENT_WEIGHTS
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    data = {}
    for net in NETS:
        for gen, (sig_name, bkg_name) in GENERATORS.items():
            sig, bkg = jets(sig_name, net, args.from_npz), jets(bkg_name, net, args.from_npz)
            if sig is not None and bkg is not None:
                data[net, gen] = (sig, bkg)
    if not data:
        sys.exit(f'no scores found for {", ".join(NETS)} ({C.CAMPAIGN})')

    # ── tables ───────────────────────────────────────────────────────────────
    if weighted:
        print_table(data, True, auc_of, f'{C.CAMPAIGN}: event-weighted AUC')
        print_table(data, True, eff_at, f'\n{C.CAMPAIGN}: event-weighted signal efficiency at {EPS_B:.0%} '
                                         'background efficiency')
        print_table(data, False, auc_of, f'\n{C.CAMPAIGN}: unweighted AUC (for comparison)')
        print_table(data, False, eff_at, f'\n{C.CAMPAIGN}: unweighted signal efficiency at {EPS_B:.0%} '
                                          'background efficiency (for comparison)')
        print('\neffective numbers of jets N_eff = (sum w)^2 / sum w^2 (jets in parentheses)')
        for gen in GENERATORS:
            nets = [n for n, g in data if g == gen]
            if not nets:
                continue
            sig, bkg = data[nets[0], gen]
            print(f'  {gen:10s} ' + '  '.join(
                f'{c}: VBF {C.n_eff(sig["kw"][sel(sig)]):.0f} ({int(sel(sig).sum())}), '
                f'QCD {C.n_eff(bkg["kw"][sel(bkg)]):.0f} ({int(sel(bkg).sum())})'
                for c, sel in CATS.items() if c in PAIRED))
    else:
        print_table(data, False, auc_of)
        print_table(data, False, eff_at, f'\nsignal efficiency at {EPS_B:.0%} background efficiency')

    for name, metric in METRICS.items():
        for gen in GENERATORS:
            if all((net, gen) in data for net in NETS):
                for cat in PAIRED:
                    d = np.subtract(per_model(*data['ParTj', gen], cat, weighted, metric),
                                    per_model(*data['NNjPF', gen], cat, weighted, metric))
                    print(f'{name:9s} {gen:10s} {cat:11s} ParTj - NNjPF per model: '
                          + ' '.join(f'{x:+.4f}' for x in d)
                          + f'   mean {d.mean():+.4f} +- {d.std(ddof=1) / np.sqrt(len(d)):.4f}')

    # ── figure: one column per category ──────────────────────────────────────
    nets_in = [n for n in NETS if any((n, g) in data for g in GENERATORS)]
    both = all((n, g) in data for n in NETS for g in GENERATORS)
    nrow = 2 if both else 1
    fig, axes = plt.subplots(nrow, len(PAIRED), sharex='col', squeeze=False, facecolor=SURFACE,
                             figsize=(5.2 * len(PAIRED), 6.8 if both else 5.0),
                             gridspec_kw={'height_ratios': [3, 1.4]} if both else None)
    for col, cat in enumerate(PAIRED):
        ax = axes[0, col]
        curves = {}
        for (net, gen), (sig, bkg) in data.items():
            ms, mb = CATS[cat](sig), CATS[cat](bkg)
            tpr, fpr = roc_points(sig, bkg, ms, mb, weighted)
            curves[net, gen] = (tpr, fpr)
            ax.plot(tpr, 1 / fpr, color=C_NET[net], lw=1.8, ls=LS[gen],
                    label=f'{net}, {gen}: AUC {auc_of(sig, bkg, ms, mb, weighted):.3f}, '
                          f'ε_S(1%) {eff_at(sig, bkg, ms, mb, weighted):.3f}')
        ax.axhline(1 / EPS_B, color=INK2, lw=0.9, ls=':')
        ax.text(0.98, 1 / EPS_B * 1.12, 'ε_B = 1 %', color=INK2, fontsize=8, ha='right', va='bottom',
                transform=ax.get_yaxis_transform())
        ax.set_yscale('log')
        ax.set_title(PANEL_TITLE[cat], fontsize=10, color=INK)
        ax.legend(frameon=False, fontsize=7.5, loc='upper right')
        if col == 0:
            ax.set_ylabel('background rejection 1 / ε_B (QCD Z+jj jets)', color=INK)
        if both:
            axr = axes[1, col]
            axr.axhline(1, color=INK2, lw=0.8)
            for gen in GENERATORS:
                (t1, f1), (t2, f2) = curves['NNjPF', gen], curves['ParTj', gen]
                grid = np.linspace(max(t1[0], t2[0], 0.05), 1, 400)
                axr.plot(grid, rejection_at(t2, f2, grid) / rejection_at(t1, f1, grid), color=INK, lw=1.5,
                         ls=LS[gen], label=gen)
            if col == 0:
                axr.set_ylabel('rejection\nParTj / NNjPF', color=INK)
            axr.legend(frameon=False, fontsize=8, loc='upper right')
    for a in axes.flat:
        a.set_xlim(0, 1)
        a.grid(color=RULE, lw=0.4, alpha=0.7, which='both')
        a.tick_params(colors=INK2)
        for sp in a.spines.values():
            sp.set_color(RULE)
    for col in range(len(PAIRED)):
        axes[-1, col].set_xlabel('signal efficiency ε_S (VBF Z jets)', color=INK)
    fig.suptitle(f'Out-of-fold ROC curves, {C.CAMPAIGN} ({", ".join(nets_in)})'
                 + (', event-weighted' if weighted else ''), fontsize=11, color=INK)
    fig.tight_layout()
    os.makedirs(C.FIGS, exist_ok=True)
    out = os.path.join(C.FIGS, 'roc_nets')
    for ext in ('pdf', 'png'):
        fig.savefig(f'{out}.{ext}', dpi=150, facecolor=SURFACE)
    print(f'wrote {out}.pdf and .png')


if __name__ == '__main__':
    main()
