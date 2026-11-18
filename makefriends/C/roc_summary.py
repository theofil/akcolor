#!/usr/bin/env python3
"""ROC summary of the single-jet nets across campaigns (lxplus; reads stored scores, runs no model).

    python3 roc_summary.py                  # tables for roc.summary.md, 400 bootstrap replicas (~10 min)
    python3 roc_summary.py --replicas 50    # quicker, rougher uncertainties

Nets: summer26 NNj; nov26 NNjPF and ParTj; fall26 NNjPF and ParTj, and NNjPF_py (NNjPF trained on
MG5+Pythia). ParTj - NNjPF is also printed on the same jets (paired: same events, same replicas).
Figures of merit, VBF Z jets (signal) against QCD Z+jj jets (background):
  AUC
  eps_S(1%)    signal efficiency at 1 % background efficiency
  eps_B(50%)   background efficiency at 50 % signal efficiency
each event-weighted (kWeight) and unweighted, with a bootstrap sigma: every event gets a Poisson(1)
multiplicity, shared by both of its jets, independently for the signal and the background sample.

Categories:
  VBF presel.  both stored jets of the events passing the summer26 selection: two jets with pT > 30 GeV,
               |eta| < 3, eta1 * eta2 < 0 (the nov26 gate), mjj > 400 GeV, |dYjj| > 2.5. The summer26
               friends hold only these events, so this is the one category all three campaigns share.
  gate jets    both stored jets of gate events (nov26 and fall26 only)
Scores:
  summer26 NNj   h5 columns NNj_jet0 / NNj_jet1 of one Herwig-trained model. Herwig: only the events of
                 the 15 % test split (summer26/NNj/split_indices_Z.npz, indices into [QCD, VBF]); the
                 model never saw either of their jets. MG5+Pythia: all events.
  nov26, fall26  out-of-fold h5 columns (nov26/save_scores.py), read with nov26/roc_scores.jets.
Paired check (also alone: --paired-only): nov26 NNjPF and ParTj - summer26 NNj on the same events. summer26 and
nov26 read the same generator files; events are matched by their boson variables.
Output: the tables of roc.summary.md (written by hand from this printout).
"""
import argparse
import os
import sys
import time

import h5py
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(BASE, 'nov26'))
import common as C  # noqa: E402
import roc_scores as R  # noqa: E402
from bootstrap_scores import event_index, trapezoid  # noqa: E402

GENS = ('herwig', 'pythia')
GEN_LABEL = {'herwig': 'Herwig', 'pythia': 'MG5+Pythia'}
SAMPLES = {'herwig': ('VBFZ_herwig', 'QCDZjj_herwig'), 'pythia': ('VBFZ_mg5_pythia', 'QCDZjj_mg5_pythia')}
METRICS = ('AUC', 'eps_S(1%)', 'eps_B(50%)')
# (campaign, h5 column, label); summer26 is read by summer26_jets
ROWS = [('summer26', 'NNj', 'NNj'), ('nov26', 'NNjPF', 'NNjPF'), ('nov26', 'ParTj', 'ParTj'),
        ('fall26', 'NNjPF', 'NNjPF'), ('fall26', 'ParTj', 'ParTj'),
        ('fall26', 'NNjPF_py', 'NNjPF, trained on MG5+Pythia')]
CATS = ('VBF presel.', 'gate jets')


def summer26_jets(sample):
    """Both jets of the summer26 NNj evaluation events of `sample`: score, kWeight, event row."""
    path = os.path.join(BASE, 'friends', 'summer26', f'{sample}.h5')
    with h5py.File(path, 'r') as f:
        s0, s1, kw = f['NNj_jet0'][:], f['NNj_jet1'][:], f['kWeight'][:].astype(np.float64)
    rows = np.arange(len(kw))
    if sample.endswith('_herwig'):                                    # Herwig: the test split only
        test = np.load(os.path.join(BASE, 'summer26', 'NNj', 'split_indices_Z.npz'))['test_idx']
        with h5py.File(os.path.join(BASE, 'friends', 'summer26', 'QCDZjj_herwig.h5'), 'r') as f:
            n_bkg = f['kWeight'].shape[0]
        rows = np.sort(test[test < n_bkg] if sample.startswith('QCD') else test[test >= n_bkg] - n_bkg)
    return dict(score=np.r_[s0[rows], s1[rows]], kw=np.r_[kw[rows], kw[rows]], ev=np.r_[rows, rows],
                slot=np.r_[np.zeros(len(rows), np.int8), np.ones(len(rows), np.int8)]), len(kw)


def campaign_jets(campaign, column, sample, cat):
    """Jets of category `cat` from the out-of-fold column `column` of a nov26/fall26 sample."""
    C.set_campaign(campaign)
    j = R.jets(sample, column)
    if j is None:
        sys.exit(f'{campaign}: no column {column} in {sample}.h5')
    ev, n = event_index(sample)
    sel = (j['role'] != C.ROLE_C) & (j['vbf'] if cat == 'VBF presel.' else True)
    return dict(score=j['score'][sel], kw=j['kw'][sel].astype(np.float64), ev=ev[sel], model=j['model'][sel]), n


def prepare(sig, bkg, weighted):
    """Jets sorted once by descending score, with the last index of each tie group."""
    s = np.r_[sig['score'], bkg['score']]
    y = np.r_[np.ones(len(sig['score']), bool), np.zeros(len(bkg['score']), bool)]
    w = np.r_[sig['kw'], bkg['kw']] if weighted else np.ones(len(s))
    ev = np.r_[sig['ev'], bkg['ev']]
    o = np.argsort(-s, kind='stable')
    last = np.r_[np.nonzero(np.diff(s[o]))[0], len(s) - 1]
    return y[o], w[o], ev[o], last


def metrics(p, m_s=None, m_b=None):
    """(AUC, eps_S at eps_B = 1 %, eps_B at eps_S = 50 %) of the (weighted) ROC curve."""
    y, w, ev, last = p
    if m_s is not None:                       # signal and background rows index different files
        m = np.empty(len(y))
        m[y], m[~y] = m_s[ev[y]], m_b[ev[~y]]
        w = w * m
    cs = np.cumsum(np.where(y, w, 0.))[last]
    cb = np.cumsum(np.where(y, 0., w))[last]
    tpr, fpr = np.r_[0., cs / cs[-1]], np.r_[0., cb / cb[-1]]
    return np.array([trapezoid(tpr, fpr), np.interp(0.01, fpr, tpr), np.interp(0.5, tpr, fpr)])


def boson_key(f, rows=None):
    """One 16-byte key per event: the float32 bit patterns of bosonPt, bosonEta, bosonPhi, bosonM."""
    cols = [f[k][:] if rows is None else f[k][:][rows] for k in ('bosonPt', 'bosonEta', 'bosonPhi', 'bosonM')]
    return np.ascontiguousarray(np.stack(cols, axis=1).astype(np.float32)).view('V16').ravel()


def paired_summer26_nov26(rng, replicas):
    """NNjPF and ParTj (nov26) - NNj (summer26) on the same events: the summer26 NNj evaluation events
    (Herwig: test split; MG5+Pythia: all) that also pass the selection in nov26, matched by the
    boson variables (the same input events, written in input order). The two campaigns cluster
    slightly different particles, so the jets themselves can differ a little. Both jets of each
    event; the nov26 kWeight for both nets; shared bootstrap replicas."""
    print('\nPaired check: nov26 NNjPF and ParTj - summer26 NNj on the same events (VBF presel., both jets)')
    for gen in GENS:
        prep, n_ev = {}, {}
        for sample in SAMPLES[gen]:
            s26, _ = summer26_jets(sample)
            rows26 = s26['ev'][s26['slot'] == 0]                         # summer26 evaluation events
            with h5py.File(os.path.join(BASE, 'friends', 'summer26', f'{sample}.h5'), 'r') as f:
                key26 = boson_key(f, rows26)
                nnj = np.stack([f['NNj_jet0'][:][rows26], f['NNj_jet1'][:][rows26]], axis=1)
            with h5py.File(os.path.join(BASE, 'friends', 'nov26', f'{sample}.h5'), 'r') as f:
                key11 = boson_key(f)
                vbf = f['gate'][:] & (f['mjj'][:] > 400) & (np.abs(f['dYjj'][:]) > 2.5)
                nnjpf, partj, kw = f['NNjPF'][:], f['ParTj'][:], f['kWeight'][:].astype(np.float64)
            uniq, cnt = np.unique(key11, return_counts=True)
            lookup = dict(zip(key11.tolist(), range(len(key11))))
            row11 = np.array([lookup.get(k, -1) for k in key26.tolist()])
            dup = np.isin(key26, uniq[cnt > 1])
            ok = (row11 >= 0) & ~dup
            ok[ok] &= vbf[row11[ok]]
            print(f'  {sample:18s} {len(rows26)} summer26 events: {int((row11 >= 0).sum())} found in nov26 '
                  f'({int(dup.sum())} ambiguous), {int(ok.sum())} also pass the nov26 selection')
            r = row11[ok]
            n_ev[sample] = len(key11)
            prep[sample] = dict(r=r, nnj=nnj[ok], nnjpf=nnjpf[r], partj=partj[r], kw=kw[r])
        sig, bkg = (prep[s] for s in SAMPLES[gen])
        res = {}
        for net, col in (('NNj', 'nnj'), ('NNjPF', 'nnjpf'), ('ParTj', 'partj')):
            jets = [dict(score=d[col].T.ravel(), kw=np.r_[d['kw'], d['kw']], ev=np.r_[d['r'], d['r']]) for d in (sig, bkg)]
            res[net] = prepare(*jets, True)
        central = {net: metrics(p) for net, p in res.items()}
        reps = {net: [] for net in res}
        for _ in range(replicas):
            m_s = rng.poisson(1., n_ev[SAMPLES[gen][0]]).astype(np.float64)
            m_b = rng.poisson(1., n_ev[SAMPLES[gen][1]]).astype(np.float64)
            for net, p in res.items():
                reps[net].append(metrics(p, m_s, m_b))
        reps = {net: np.array(v) for net, v in reps.items()}
        for i, name in enumerate(METRICS):
            line = f'  {GEN_LABEL[gen]:10s} {name:10s} ' + '  '.join(f'{n} {central[n][i]:.4f}' for n in central)
            for net in ('NNjPF', 'ParTj'):
                d = central[net][i] - central['NNj'][i]
                s = (reps[net][:, i] - reps['NNj'][:, i]).std(ddof=1)
                line += f'   {net} - NNj {d:+.4f} ± {s:.4f}'
            print(line + '  (event-weighted)')


def make_figure():
    """ROC curves, transformer (ParTj) against DeepSets (NNjPF; summer26 NNj as reference), event-weighted
    -> figs/roc_summary.{pdf,png}. Rows: VBF-preselected jets, all gate jets (each with a ratio panel);
    columns: tested on Herwig, on MG5+Pythia."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    col_of = {'summer26': R.INK2, 'nov26': R.C_NET['NNjPF'], 'fall26': R.C_NET['ParTj']}    # grey, slots 1-2
    ls_of = {'NNj': '-', 'NNjPF': '-', 'ParTj': '--'}
    name_of = {'NNj': 'DeepSets (NNj)', 'NNjPF': 'DeepSets (NNjPF)', 'ParTj': 'transformer (ParTj)'}
    curves = [('summer26', 'NNj'), ('nov26', 'NNjPF'), ('nov26', 'ParTj'), ('fall26', 'NNjPF'), ('fall26', 'ParTj')]
    title_of = {'VBF presel.': 'VBF-preselected jets', 'gate jets': 'all gate jets'}
    fig, axes = plt.subplots(4, 2, sharex='col', figsize=(11.5, 13), facecolor=R.SURFACE,
                             gridspec_kw={'height_ratios': [3, 1.3, 3, 1.3]})
    for col, gen in enumerate(GENS):
        sig_name, bkg_name = SAMPLES[gen]
        for block, cat in enumerate(CATS):
            ax, axr = axes[2 * block, col], axes[2 * block + 1, col]
            roc = {}
            for campaign, net in curves:
                if campaign == 'summer26':
                    if cat != 'VBF presel.':
                        continue
                    sig, bkg = summer26_jets(sig_name)[0], summer26_jets(bkg_name)[0]
                else:
                    sig, bkg = (campaign_jets(campaign, net, s, cat)[0] for s in (sig_name, bkg_name))
                a, e, _ = metrics(prepare(sig, bkg, True))
                tpr, fpr = R.roc_points(sig, bkg, np.ones(len(sig['score']), bool),
                                        np.ones(len(bkg['score']), bool), True)
                roc[campaign, net] = (tpr, fpr)
                ax.plot(tpr, 1 / fpr, color=col_of[campaign], ls=ls_of[net], lw=1.7,
                        label=f'{campaign} {name_of[net]}: AUC {a:.3f}, ε_S(1%) {e:.3f}')
            for campaign in ('nov26', 'fall26'):
                (t1, f1), (t2, f2) = roc[campaign, 'NNjPF'], roc[campaign, 'ParTj']
                grid = np.linspace(max(t1[0], t2[0], 0.05), 1, 400)
                axr.plot(grid, R.rejection_at(t2, f2, grid) / R.rejection_at(t1, f1, grid),
                         color=col_of[campaign], lw=1.5, label=campaign)
            ax.axhline(1 / R.EPS_B, color=R.INK2, lw=0.9, ls=':')
            ax.text(0.98, 1 / R.EPS_B * 1.12, 'ε_B = 1 %', color=R.INK2, fontsize=8, ha='right', va='bottom',
                    transform=ax.get_yaxis_transform())
            ax.set_yscale('log')
            ax.set_title(f'{title_of[cat]}, tested on {GEN_LABEL[gen]}', fontsize=10, color=R.INK)
            ax.legend(frameon=False, fontsize=7, loc='upper right')
            axr.axhline(1, color=R.INK2, lw=0.8)
            axr.legend(frameon=False, fontsize=7.5, loc='best')
            if col == 0:
                ax.set_ylabel('background rejection 1 / ε_B', color=R.INK)
                axr.set_ylabel('rejection ratio\ntransformer / DeepSets', color=R.INK)
    for a in axes.flat:
        a.set_xlim(0, 1)
        a.grid(color=R.RULE, lw=0.4, alpha=0.7, which='both')
        a.tick_params(colors=R.INK2)
        for sp in a.spines.values():
            sp.set_color(R.RULE)
    for col in range(2):
        axes[-1, col].set_xlabel('signal efficiency ε_S (VBF Z jets)', color=R.INK)
    fig.suptitle('Transformer (ParTj, dashed) against DeepSets (NNjPF, solid; summer26 NNj in grey):\n'
                 'out-of-fold ROC curves, event-weighted', fontsize=11, color=R.INK)
    fig.tight_layout()
    out = os.path.join(BASE, 'figs', 'roc_summary')
    for ext in ('pdf', 'png'):
        fig.savefig(f'{out}.{ext}', dpi=150, facecolor=R.SURFACE)
    print(f'wrote {out}.pdf and .png')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--figure-only', action='store_true', help='make only the ROC figure figs/roc_summary.{pdf,png}')
    ap.add_argument('--paired-only', action='store_true', help='run only the paired summer26/nov26 check')
    ap.add_argument('--replicas', type=int, default=400)
    ap.add_argument('--seed', type=int, default=20261118)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    if args.paired_only:
        paired_summer26_nov26(rng, args.replicas)
        return
    if args.figure_only:
        make_figure()
        return

    # check: the summer26 NNj AUCs quoted in summer26.md (jet 0 only, unweighted)
    from sklearn.metrics import roc_auc_score
    for gen in GENS:
        (sig, _), (bkg, _) = (summer26_jets(s) for s in SAMPLES[gen])
        m_s, m_b = sig['slot'] == 0, bkg['slot'] == 0
        a = roc_auc_score(np.r_[np.ones(m_s.sum()), np.zeros(m_b.sum())], np.r_[sig['score'][m_s], bkg['score'][m_b]])
        print(f'check: summer26 NNj jet 0, {GEN_LABEL[gen]}, unweighted AUC {a:.4f}')

    res, stats, reps_all, permodel = {}, {}, {}, {}
    for campaign in ('summer26', 'nov26', 'fall26'):
        for gen in GENS:
            sig_name, bkg_name = SAMPLES[gen]
            prep, n_s, n_b = {}, None, None
            for camp, column, label in ROWS:
                if camp != campaign:
                    continue
                for cat in CATS:
                    if campaign == 'summer26':
                        if cat != 'VBF presel.':
                            continue
                        (sig, n_s), (bkg, n_b) = summer26_jets(sig_name), summer26_jets(bkg_name)
                    else:
                        (sig, n_s), (bkg, n_b) = (campaign_jets(campaign, column, s, cat) for s in (sig_name, bkg_name))
                    stats[campaign, label, gen, cat] = (len(sig['score']), len(bkg['score']),
                                                       C.n_eff(sig['kw']), C.n_eff(bkg['kw']))
                    for weighted in (True, False):
                        prep[label, cat, weighted] = prepare(sig, bkg, weighted)
                    if 'model' in sig:                  # each fold model on its own inference jets (weighted)
                        permodel[campaign, label, gen, cat] = np.array([
                            metrics(prepare(C.take(sig, sig['model'] == k), C.take(bkg, bkg['model'] == k), True))
                            for k in range(C.N_FOLDS)])
            reps = {key: [] for key in prep}
            for _ in range(args.replicas):
                m_s = rng.poisson(1., n_s).astype(np.float64)
                m_b = rng.poisson(1., n_b).astype(np.float64)
                for key, p in prep.items():
                    reps[key].append(metrics(p, m_s, m_b))
            for (label, cat, weighted), p in prep.items():
                reps_all[campaign, label, gen, cat, weighted] = np.array(reps[label, cat, weighted])
                res[campaign, label, gen, cat, weighted] = (metrics(p), reps_all[campaign, label, gen, cat, weighted].std(0, ddof=1))
            print(f'{campaign} {gen}: done ({time.time() - t0:.0f} s)')

    for cat in CATS:
        for weighted in (True, False):
            print(f'\n### {cat}, {"event-weighted" if weighted else "unweighted"} ({args.replicas} bootstrap replicas)\n')
            print('| Campaign | Net | ' + ' | '.join(f'{m} {GEN_LABEL[g]}' for g in GENS for m in METRICS) + ' |')
            print('|---|---|' + '---|' * 6)
            for campaign, column, label in ROWS:
                if (campaign, label, 'herwig', cat, weighted) not in res:
                    continue
                cells = []
                for gen in GENS:
                    v, s = res[campaign, label, gen, cat, weighted]
                    cells += [f'{v[0]:.4f} ± {s[0]:.4f}', f'{v[1]:.4f} ± {s[1]:.4f}', f'{v[2]:.4f} ± {s[2]:.4f}']
                print(f'| {campaign} | {label} | ' + ' | '.join(cells) + ' |')

    print('\n### ParTj - NNjPF on the same jets, event-weighted: pooled difference ± the larger of the bootstrap '
          'sigma and the standard error of the five per-model differences\n')
    print('| Jets | Tested on | Campaign | ' + ' | '.join(METRICS) + ' |')
    print('|---|---|---|---|---|---|')
    for cat in CATS:
        for gen in GENS:
            for campaign in ('nov26', 'fall26'):
                cells = []
                for i in range(len(METRICS)):
                    a, b = (campaign, 'ParTj', gen, cat, True), (campaign, 'NNjPF', gen, cat, True)
                    d = res[a][0][i] - res[b][0][i]
                    sb = (reps_all[a][:, i] - reps_all[b][:, i]).std(ddof=1)
                    pm = permodel[campaign, 'ParTj', gen, cat][:, i] - permodel[campaign, 'NNjPF', gen, cat][:, i]
                    se = pm.std(ddof=1) / np.sqrt(len(pm))
                    cells.append(f'{d:+.4f} ± {max(sb, se):.4f}')
                print(f'| {cat} | {GEN_LABEL[gen]} | {campaign} | ' + ' | '.join(cells) + ' |')

    print('\njets (VBF / QCD) and N_eff (VBF / QCD, event weights)')
    for (campaign, label, gen, cat), (ns, nb, es, eb) in stats.items():
        print(f'  {campaign:8s} {label:30s} {GEN_LABEL[gen]:10s} {cat:11s} {ns:8d} / {nb:7d}   N_eff {es:9.0f} / {eb:8.0f}')
    paired_summer26_nov26(rng, args.replicas)
    make_figure()
    print(f'({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
