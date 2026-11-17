#!/usr/bin/env python3
"""Out-of-fold performance of the nov26 jet nets, from the scores in the h5 files (lxplus; no model).

    python3 nov26/roc_scores.py            # every net with scores -> figs/nov26/roc_nets.{pdf,png}
    python3 nov26/roc_scores.py --from-npz # read <net>/scores_fold*.npz instead of the h5 columns

For each net and generator (Herwig: the training generator; MG5+Pythia: transfer test) it prints
the AUC of VBF jets vs QCD jets for
  gate jets      roles A + B, both stored jets of gate events
  role A, role B the 1st (slot 0) and 2nd (slot 1) jet of gate events, separately; the slots are
                 randomised and both are out-of-fold, so the two AUCs must agree within statistics
  gate == 0      role C
  VBF presel.    gate jets of events with mjj > 400 GeV and |dYjj| > 2.5 (summer26 selection)
and the spread of the five fold models (Herwig gate jets, each model on its own inference
events). AUCs are unweighted, as in summer26. The figure shows the background rejection vs the
signal efficiency of the VBF-preselected gate jets, down to 10 surviving background jets; when
both nets have scores, a lower panel shows the ratio of their rejections at equal signal
efficiency (the nets scored the same jets, so the ratio is free of sample-to-sample noise).
"""
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


def rejection_at(tpr, fpr, grid):
    """Background rejection 1/fpr of a ROC curve, interpolated at the signal efficiencies grid."""
    t, i = np.unique(tpr, return_index=True)
    return 1 / np.interp(grid, t, fpr[i])


def jets(sample, net, from_npz=False):
    """Per-jet score, role, VBF-preselection flag and model (fold) index of one sample."""
    with h5py.File(C.h5_path(sample), 'r') as f:
        n_jets, gate = f['nJets'][:], f['gate'][:]
        vbf = gate & (f['mjj'][:] > 400) & (np.abs(f['dYjj'][:]) > 2.5)
        if from_npz:
            folds_dir = os.path.join(HERE, net)
            if not all(os.path.exists(os.path.join(folds_dir, f'scores_fold{k}.npz')) for k in range(C.N_FOLDS)):
                return None
            sub = C.assign_subfolds(sample, gate, n_jets)
            score = save_scores.merge(sample, n_jets, sub, folds_dir)
        elif net in f:
            score, sub = f[net][:], f['subfold'][:]
        else:
            return None
    model = sub // 2                                                  # model k scored sub-folds 2k, 2k+1
    out = {k: [] for k in ('score', 'role', 'vbf', 'model')}
    for s in (0, 1):
        m = s < n_jets
        out['score'].append(score[m, s])
        out['role'].append(np.where(gate[m], C.ROLE_A if s == 0 else C.ROLE_B, C.ROLE_C))
        out['vbf'].append(vbf[m])
        out['model'].append(model[m])
    return {k: np.concatenate(v) for k, v in out.items()}


def main():
    from_npz = '--from-npz' in sys.argv[1:]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    cats = {'gate jets': lambda j: j['role'] != C.ROLE_C, 'role A': lambda j: j['role'] == C.ROLE_A,
            'role B': lambda j: j['role'] == C.ROLE_B,
            'gate == 0': lambda j: j['role'] == C.ROLE_C,
            'VBF presel.': lambda j: (j['role'] != C.ROLE_C) & j['vbf']}
    curves = {}
    print(f'{"net":6s} {"generator":10s} ' + ' '.join(f'{c:>12s}' for c in cats) + '   fold-model spread')
    for net in NETS:
        for gen, (sig_name, bkg_name) in GENERATORS.items():
            sig, bkg = jets(sig_name, net, from_npz), jets(bkg_name, net, from_npz)
            if sig is None or bkg is None:
                continue
            row = []
            for name, sel in cats.items():
                ms, mb = sel(sig), sel(bkg)
                y = np.r_[np.ones(ms.sum()), np.zeros(mb.sum())]
                s = np.r_[sig['score'][ms], bkg['score'][mb]]
                row.append(roc_auc_score(y, s))
                if name == 'VBF presel.':
                    fpr, tpr, _ = roc_curve(y, s)
                    keep = fpr * mb.sum() >= 10                    # at least 10 background jets left
                    curves[net, gen] = (tpr[keep], fpr[keep], row[-1])
            spread = ''
            if gen == 'Herwig':
                per = []
                for k in range(C.N_FOLDS):
                    ms, mb = cats['gate jets'](sig) & (sig['model'] == k), cats['gate jets'](bkg) & (bkg['model'] == k)
                    per.append(roc_auc_score(np.r_[np.ones(ms.sum()), np.zeros(mb.sum())],
                                             np.r_[sig['score'][ms], bkg['score'][mb]]))
                spread = f'   {" ".join(f"{a:.4f}" for a in per)}  (rms {np.std(per):.4f})'
            print(f'{net:6s} {gen:10s} ' + ' '.join(f'{a:12.4f}' for a in row) + spread)
    both = all((n, g) in curves for n in NETS for g in GENERATORS)
    fig, axes = plt.subplots(2 if both else 1, 1, sharex=True, squeeze=False, facecolor=SURFACE,
                             figsize=(6.4, 6.8 if both else 5.0),
                             gridspec_kw={'height_ratios': [3, 1.4]} if both else None)
    ax = axes[0, 0]
    for (net, gen), (tpr, fpr, auc) in curves.items():
        ax.plot(tpr, 1 / fpr, color=C_NET[net], lw=1.8, ls=LS[gen], label=f'{net}, {gen}: AUC {auc:.3f}')
    ax.set_yscale('log')
    ax.set_ylabel('background rejection 1 / ε (QCD Z+jj jets)', color=INK)
    ax.set_title('Out-of-fold ROC, gate jets of VBF-preselected events\n(mjj > 400 GeV, |ΔY_jj| > 2.5)',
                 fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8.5)
    if both:                                     # ParTj / NNjPF rejection at equal signal efficiency
        axr = axes[1, 0]
        axr.axhline(1, color=INK2, lw=0.8)
        for gen in GENERATORS:
            (t1, f1, _), (t2, f2, _) = curves['NNjPF', gen], curves['ParTj', gen]
            grid = np.linspace(max(t1[0], t2[0], 0.05), 1, 400)
            axr.plot(grid, rejection_at(t2, f2, grid) / rejection_at(t1, f1, grid), color=INK, lw=1.5,
                     ls=LS[gen], label=gen)
        axr.set_ylabel('rejection\nParTj / NNjPF', color=INK)
        axr.legend(frameon=False, fontsize=8.5)
    for a in axes[:, 0]:
        a.set_xlim(0, 1)
        a.grid(color=RULE, lw=0.4, alpha=0.7, which='both')
        a.tick_params(colors=INK2)
        for sp in a.spines.values():
            sp.set_color(RULE)
    axes[-1, 0].set_xlabel('signal efficiency (VBF Z jets)', color=INK)
    fig.tight_layout()
    out = os.path.join(C.BASE, 'figs', 'nov26', 'roc_nets')
    for ext in ('pdf', 'png'):
        fig.savefig(f'{out}.{ext}', dpi=150, facecolor=SURFACE)
    print(f'wrote {out}.pdf and .png')


if __name__ == '__main__':
    main()
