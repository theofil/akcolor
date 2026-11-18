#!/usr/bin/env python3
"""Statistical uncertainties of the out-of-fold figures of merit, by bootstrap (lxplus; no model).

    python3 nov26/bootstrap_scores.py --campaign fall26          # 400 replicas, ~6 min
    python3 nov26/bootstrap_scores.py --campaign nov26 --replicas 200

Reads the score columns of friends/<campaign>/*.h5, like roc_scores.py, and prints:
  1. the background rejection 1/eps_B of NNjPF at fixed signal efficiencies, and the ratio
     ParTj / NNjPF, for the gate jets and the VBF-preselected gate jets (the ROC curves of
     figs/<campaign>/roc_nets.pdf), with the number of QCD jets above the NNjPF cut at eps_S = 5 %;
  2. for every category of roc_scores.py, and for each fold model's own gate jets, the AUC and
     eps_S(1 %) with their bootstrap standard deviation;
  3. the bootstrap standard deviation of the differences role A - role B and ParTj - NNjPF.

Bootstrap: each replica gives every event a multiplicity drawn from a Poisson distribution of
mean 1, independently for the signal and the background sample. The multiplicity multiplies the
event's weight (kWeight for fall26, 1 for nov26). Both jets of an event get the same multiplicity,
because they share one event weight, and both nets see the same replicas, so the differences
are resampled too. The spread of the replicas is the statistical uncertainty for the trained
models as they are; it does not include the randomness of the training (the spread of the five
fold models does).
"""
import argparse
import os
import sys
import time

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import roc_scores as R  # noqa: E402

GRID = np.array([0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90])     # signal efficiencies of the rejection table
trapezoid = getattr(np, 'trapezoid', None) or np.trapz


def event_index(sample):
    """Row (event) index of every jet, in the jet order of roc_scores.jets."""
    with h5py.File(C.h5_path(sample), 'r') as f:
        n_jets = f['nJets'][:]
    return np.concatenate([np.nonzero(s < n_jets)[0] for s in (0, 1)]), len(n_jets)


def rejection_table(d, cat, weighted):
    rej = {}
    for net in R.NETS:
        sig, bkg = d[net]
        tpr, fpr = R.roc_points(sig, bkg, R.CATS[cat](sig), R.CATS[cat](bkg), weighted)
        rej[net] = R.rejection_at(tpr, fpr, GRID)
    sig, bkg = d['NNjPF']
    ms, mb = R.CATS[cat](sig), R.CATS[cat](bkg)
    ss = sig['score'][ms]
    ws = sig['kw'][ms].astype(np.float64) if weighted else np.ones(ms.sum())
    o = np.argsort(-ss)
    cut5 = ss[o][np.searchsorted(np.cumsum(ws[o]) / ws.sum(), 0.05)]
    return rej, int((bkg['score'][mb] >= cut5).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--campaign', default='nov26', choices=sorted(C.CAMPAIGNS))
    ap.add_argument('--replicas', type=int, default=400)
    ap.add_argument('--seed', type=int, default=20261118)
    args = ap.parse_args()
    C.set_campaign(args.campaign)
    weighted = C.EVENT_WEIGHTS
    rng = np.random.default_rng(args.seed)

    cats = dict(R.CATS)
    for k in range(C.N_FOLDS):
        cats[f'gate jets, model {k}'] = (lambda k: lambda j: (j['role'] != C.ROLE_C) & (j['model'] == k))(k)

    for gen, (sig_name, bkg_name) in R.GENERATORS.items():
        t0 = time.time()
        d = {net: (R.jets(sig_name, net), R.jets(bkg_name, net)) for net in R.NETS}
        if any(s is None or b is None for s, b in d.values()):
            sys.exit(f'{C.CAMPAIGN}: both nets need scores in the h5 files (run save_scores.py)')
        ev_s, n_s = event_index(sig_name)
        ev_b, n_b = event_index(bkg_name)
        assert len(ev_s) == len(d['NNjPF'][0]['score']) and len(ev_b) == len(d['NNjPF'][1]['score'])

        # 1. rejection at fixed signal efficiency
        print(f'{C.CAMPAIGN}, {gen}: background rejection 1/eps_B at fixed signal efficiency'
              + (' (event-weighted)' if weighted else ''))
        print(f'  {"eps_S":28s}' + ''.join(f'{g:8.2f}' for g in GRID))
        for cat in ('gate jets', 'VBF presel.'):
            rej, nb5 = rejection_table(d, cat, weighted)
            print(f'  {cat:11s} NNjPF rejection ' + ''.join(f'{x:8.1f}' for x in rej['NNjPF']))
            print(f'  {cat:11s} ParTj / NNjPF   ' + ''.join(f'{x:8.3f}' for x in rej['ParTj'] / rej['NNjPF'])
                  + f'   ({nb5} QCD jets above the NNjPF cut at eps_S = 5 %)')

        # 2. bootstrap: per (category, net) the jets sorted once by descending score
        prep = {}
        for cat, sel in cats.items():
            for net in R.NETS:
                sig, bkg = d[net]
                ms, mb = sel(sig), sel(bkg)
                s = np.r_[sig['score'][ms], bkg['score'][mb]]
                y = np.r_[np.ones(ms.sum(), bool), np.zeros(mb.sum(), bool)]
                kw = (np.r_[sig['kw'][ms], bkg['kw'][mb]].astype(np.float64) if weighted
                      else np.ones(len(s)))
                ev = np.r_[ev_s[ms], ev_b[mb]]
                o = np.argsort(-s, kind='stable')
                last = np.r_[np.nonzero(np.diff(s[o]))[0], len(s) - 1]     # last jet of each tie group
                prep[cat, net] = (y[o], kw[o], ev[o], last)

        def metrics(cat, net, m_s=None, m_b=None):
            """(eps_S at eps_B = 1 %, AUC) of the weighted ROC curve, as sklearn's roc_curve gives it."""
            y, kw, ev, last = prep[cat, net]
            w = kw if m_s is None else kw * np.where(y, m_s[ev], m_b[ev])
            cs = np.cumsum(np.where(y, w, 0.))[last]
            cb = np.cumsum(np.where(y, 0., w))[last]
            tpr, fpr = np.r_[0., cs / cs[-1]], np.r_[0., cb / cb[-1]]
            return np.interp(R.EPS_B, fpr, tpr), trapezoid(tpr, fpr)

        central = {key: metrics(*key) for key in prep}
        reps = {key: [] for key in prep}
        for _ in range(args.replicas):
            m_s = rng.poisson(1., n_s).astype(np.float64)
            m_b = rng.poisson(1., n_b).astype(np.float64)
            for key in prep:
                reps[key].append(metrics(*key, m_s, m_b))
        reps = {k: np.array(v) for k, v in reps.items()}                 # (replicas, 2): eps_S, AUC

        print(f'\n{C.CAMPAIGN}, {gen}: figures of merit and bootstrap sigma ({args.replicas} replicas)')
        print(f'  {"category":22s} {"net":6s} {"eps_S(1%)":>10s} {"sigma":>7s} {"AUC":>8s} {"sigma":>7s}')
        for (cat, net), (e, a) in central.items():
            se, sa = reps[cat, net].std(0, ddof=1)
            print(f'  {cat:22s} {net:6s} {e:10.4f} {se:7.4f} {a:8.4f} {sa:7.4f}')
        print(f'\n{C.CAMPAIGN}, {gen}: differences and bootstrap sigma')
        diffs = [(f'{net} role A - role B', ('role A', net), ('role B', net)) for net in R.NETS]
        diffs += [(f'ParTj - NNjPF, {cat}', (cat, 'ParTj'), (cat, 'NNjPF')) for cat in R.PAIRED]
        for name, k1, k2 in diffs:
            c = np.subtract(central[k1], central[k2])
            s = (reps[k1] - reps[k2]).std(0, ddof=1)
            print(f'  {name:28s} eps_S(1%) {c[0]:+.4f} +- {s[0]:.4f}   AUC {c[1]:+.4f} +- {s[1]:.4f}')
        print(f'  ({time.time() - t0:.0f} s)\n')


if __name__ == '__main__':
    main()
