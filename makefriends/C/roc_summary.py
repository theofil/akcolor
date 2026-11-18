#!/usr/bin/env python3
"""ROC summary: summer26 NNj against the nov26 NNj and NNjBC, and the beam-connection ablation
(lxplus; reads stored scores, runs no model).

    python3 roc_summary.py                  # tables for roc.summary.md + figures, 400 bootstrap replicas
    python3 roc_summary.py --replicas 50    # quicker, rougher uncertainties
    python3 roc_summary.py --paired-only    # only the summer26 / nov26 same-event comparison

Nets: summer26 NNj; the nov26 nets of nov26/common.NETS: NNj (same DeepSets, nov26 folds), NNjBC
(NNj + the four beam-connection variables bc_in, bc_out, bc_met, bc_dphi), the add-one nets
NNjBC<v> = NNj + v and the drop-one nets NNjBC_no<v> = NNjBC - v.
Figures of merit, VBF Z jets (signal) against QCD Z+jj jets (background):
  AUC
  eps_S(1%)    signal efficiency at 1 % background efficiency
  eps_B(50%)   background efficiency at 50 % signal efficiency
event-weighted (kWeight), with a bootstrap sigma: every event gets a Poisson(1) multiplicity,
shared by both of its jets, independently for the signal and the background sample; all nov26
nets see the same replicas.

Categories:
  VBF presel.  both stored jets of the events passing the summer26 selection: two jets with pT > 30 GeV,
               |eta| < 3, eta1 * eta2 < 0 (the nov26 gate), mjj > 400 GeV, |dYjj| > 2.5. The summer26
               friends hold only these events, so this is the one category both campaigns share.
  gate jets    both stored jets of gate events (nov26 only)
Scores:
  summer26 NNj   h5 columns NNj_jet0 / NNj_jet1 of one Herwig-trained model. Herwig: only the events of
                 the 15 % test split (summer26/NNj/split_indices_Z.npz, indices into [QCD, VBF]); the
                 model never saw either of their jets. MG5+Pythia: all events.
  nov26          out-of-fold h5 columns (nov26/save_scores.py), read with nov26/roc_scores.jets.
Ablation (nov26, same jets, paired): NNjBC - NNj, the add-one gains NNjBC<v> - NNj and the drop-one
losses NNjBC - NNjBC_no<v>, each +- the larger of the bootstrap sigma of the difference and the
standard error of the five per-fold-model differences (the second includes the training randomness).
Paired check (also alone: --paired-only): nov26 NNj and NNjBC - summer26 NNj on the same events.
summer26 and nov26 read the same generator files; events are matched by their boson variables.
Output: the tables of roc.summary.md (written by hand from this printout), figs/roc_summary.{pdf,png}
(ROC curves) and figs/nov26/bc_ablation.{pdf,png} (ablation).
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
SAMPLES = {'herwig': C.TRAIN_SAMPLES, 'pythia': C.TRANSFER_SAMPLES}
METRICS = ('AUC', 'eps_S(1%)', 'eps_B(50%)')
CATS = ('VBF presel.', 'gate jets')
MAIN = ('NNj', 'NNjBC')                                   # the nov26 nets of the main tables
ADD = {v: (f'NNjBC{v}', 'NNj') for v in C.BC_VARS}        # add-one gain: NNj + v - NNj
DROP = {v: ('NNjBC', f'NNjBC_no{v}') for v in C.BC_VARS}  # drop-one loss: NNjBC - (NNjBC - v)
ALL = ('NNjBC', 'NNj')


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


def nov26_jets(net, sample, cat):
    """Jets of category `cat` from the out-of-fold column `net` of a nov26 sample."""
    j = R.jets(sample, net)
    if j is None:
        sys.exit(f'nov26: no column {net} in {sample}.h5 (run nov26/save_scores.py --net all)')
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


def fmt(x, s):
    return f'{x:.4f} ± {s:.4f}'


def paired_summer26_nov26(rng, replicas):
    """nov26 NNj and NNjBC - summer26 NNj on the same events: the summer26 NNj evaluation events
    (Herwig: test split; MG5+Pythia: all) that also pass the selection in nov26, matched by the
    boson variables (the same input events, written in input order). The two campaigns cluster
    slightly different particles, so the jets themselves can differ a little. Both jets of each
    event; the nov26 kWeight for all nets; shared bootstrap replicas."""
    print('\n### Same events: summer26 NNj against nov26 NNj and NNjBC (VBF presel., both jets, event-weighted)\n')
    nets = ('NNj (summer26)', 'NNj', 'NNjBC')
    print('| Tested on | Metric | ' + ' | '.join(nets) + ' | NNj (nov26) − NNj (summer26) | NNjBC − NNj (summer26) '
          '| NNjBC − NNj (nov26) |')
    print('|---|---|---|---|---|---|---|---|')
    for gen in GENS:
        prep, n_ev = {}, {}
        for sample in SAMPLES[gen]:
            s26, _ = summer26_jets(sample)
            rows26 = s26['ev'][s26['slot'] == 0]                         # summer26 evaluation events
            with h5py.File(os.path.join(BASE, 'friends', 'summer26', f'{sample}.h5'), 'r') as f:
                key26 = boson_key(f, rows26)
                nnj26 = np.stack([f['NNj_jet0'][:][rows26], f['NNj_jet1'][:][rows26]], axis=1)
            with h5py.File(C.h5_path(sample), 'r') as f:
                key11 = boson_key(f)
                vbf = f['gate'][:] & (f['mjj'][:] > 400) & (np.abs(f['dYjj'][:]) > 2.5)
                cols = {net: f[net][:] for net in MAIN}
                kw = f['kWeight'][:].astype(np.float64)
            uniq, cnt = np.unique(key11, return_counts=True)
            lookup = dict(zip(key11.tolist(), range(len(key11))))
            row11 = np.array([lookup.get(k, -1) for k in key26.tolist()])
            dup = np.isin(key26, uniq[cnt > 1])
            ok = (row11 >= 0) & ~dup
            ok[ok] &= vbf[row11[ok]]
            print(f'<!-- {sample}: {len(rows26)} summer26 events, {int((row11 >= 0).sum())} found in nov26 '
                  f'({int(dup.sum())} ambiguous), {int(ok.sum())} also pass the nov26 selection -->')
            r = row11[ok]
            n_ev[sample] = len(key11)
            prep[sample] = dict(r=r, kw=kw[r], **{'NNj (summer26)': nnj26[ok]}, **{net: cols[net][r] for net in MAIN})
        sig, bkg = (prep[s] for s in SAMPLES[gen])
        res = {}
        for net in nets:
            jets = [dict(score=d[net].T.ravel(), kw=np.r_[d['kw'], d['kw']], ev=np.r_[d['r'], d['r']]) for d in (sig, bkg)]
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
            cells = [f'{central[n][i]:.4f}' for n in nets]
            for a, b in (('NNj', 'NNj (summer26)'), ('NNjBC', 'NNj (summer26)'), ('NNjBC', 'NNj')):
                d = central[a][i] - central[b][i]
                s = (reps[a][:, i] - reps[b][:, i]).std(ddof=1)
                cells.append(f'{d:+.4f} ± {s:.4f}')
            print(f'| {GEN_LABEL[gen]} | {name} | ' + ' | '.join(cells) + ' |')


def make_roc_figure():
    """ROC curves of nov26 NNj and NNjBC (summer26 NNj in grey), event-weighted -> figs/roc_summary.{pdf,png}.
    Rows: VBF-preselected jets, all gate jets (each with a ratio panel NNjBC / NNj);
    columns: tested on Herwig, on MG5+Pythia."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    style = {('summer26', 'NNj'): (R.INK2, ':', 'summer26 NNj'), ('nov26', 'NNj'): (R.C_NET['NNj'], '-', 'NNj'),
             ('nov26', 'NNjBC'): (R.C_NET['NNjBC'], '-', 'NNjBC (NNj + beam connection)')}
    title_of = {'VBF presel.': 'VBF-preselected jets', 'gate jets': 'all gate jets'}
    fig, axes = plt.subplots(4, 2, sharex='col', figsize=(11.5, 13), facecolor=R.SURFACE,
                             gridspec_kw={'height_ratios': [3, 1.3, 3, 1.3]})
    for col, gen in enumerate(GENS):
        sig_name, bkg_name = SAMPLES[gen]
        for block, cat in enumerate(CATS):
            ax, axr = axes[2 * block, col], axes[2 * block + 1, col]
            roc = {}
            for (campaign, net), (color, ls, label) in style.items():
                if campaign == 'summer26':
                    if cat != 'VBF presel.':
                        continue
                    sig, bkg = summer26_jets(sig_name)[0], summer26_jets(bkg_name)[0]
                else:
                    sig, bkg = (nov26_jets(net, s, cat)[0] for s in (sig_name, bkg_name))
                a, e, _ = metrics(prepare(sig, bkg, True))
                tpr, fpr = R.roc_points(sig, bkg, np.ones(len(sig['score']), bool),
                                        np.ones(len(bkg['score']), bool), True)
                roc[campaign, net] = (tpr, fpr)
                ax.plot(tpr, 1 / fpr, color=color, ls=ls, lw=1.8, label=f'{label}: AUC {a:.3f}, ε_S(1%) {e:.3f}')
            (t1, f1), (t2, f2) = roc['nov26', 'NNj'], roc['nov26', 'NNjBC']
            grid = np.linspace(max(t1[0], t2[0], 0.05), 1, 400)
            axr.plot(grid, R.rejection_at(t2, f2, grid) / R.rejection_at(t1, f1, grid), color=R.INK, lw=1.5)
            ax.axhline(1 / R.EPS_B, color=R.INK2, lw=0.9, ls=':')
            ax.text(0.98, 1 / R.EPS_B * 1.12, 'ε_B = 1 %', color=R.INK2, fontsize=8, ha='right', va='bottom',
                    transform=ax.get_yaxis_transform())
            ax.set_yscale('log')
            ax.set_title(f'{title_of[cat]}, tested on {GEN_LABEL[gen]}', fontsize=10, color=R.INK)
            ax.legend(frameon=False, fontsize=7.5, loc='upper right')
            axr.axhline(1, color=R.INK2, lw=0.8)
            if col == 0:
                ax.set_ylabel('background rejection 1 / ε_B', color=R.INK)
                axr.set_ylabel('rejection ratio\nNNjBC / NNj', color=R.INK)
    for a in axes.flat:
        a.set_xlim(0, 1)
        a.grid(color=R.RULE, lw=0.4, alpha=0.7, which='both')
        a.tick_params(colors=R.INK2)
        for sp in a.spines.values():
            sp.set_color(R.RULE)
    for col in range(2):
        axes[-1, col].set_xlabel('signal efficiency ε_S (VBF Z jets)', color=R.INK)
    fig.suptitle('NNjBC (NNj + beam connection) against NNj: out-of-fold ROC curves, event-weighted\n'
                 '(summer26 NNj dotted grey, VBF-preselected jets only)', fontsize=11, color=R.INK)
    fig.tight_layout()
    out = os.path.join(BASE, 'figs', 'roc_summary')
    for ext in ('pdf', 'png'):
        fig.savefig(f'{out}.{ext}', dpi=150, facecolor=R.SURFACE)
    print(f'wrote {out}.pdf and .png')


def make_ablation_figure(diff):
    """Add-one gains and drop-one losses per bc variable, with NNjBC - NNj, -> figs/nov26/bc_ablation.{pdf,png}.
    diff[(net, ref), gen, cat] = (central (3,), sigma (3,)). Rows: AUC, eps_S(1%); columns: category x generator."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cols = [(cat, gen) for cat in ('gate jets', 'VBF presel.') for gen in GENS]
    rows = ((0, 'Δ AUC'), (1, 'Δ ε_S(1 %)'))
    xlab = ['all four'] + [f'bc_{v}' for v in C.BC_VARS]
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.6), sharex=True, sharey='row', facecolor=R.SURFACE)
    for c, (cat, gen) in enumerate(cols):
        for r, (i, ylabel) in enumerate(rows):
            ax = axes[r, c]
            ax.axhline(0, color=R.INK2, lw=0.9)
            v, s = diff[ALL, gen, cat]
            ax.errorbar([0], [v[i]], yerr=[s[i]], fmt='D', ms=8, color=R.INK, mec=R.SURFACE, mew=1.5,
                        elinewidth=1.4, capsize=0, label='NNjBC − NNj')
            for key, off, color, marker, label in ((ADD, -0.14, R.C_NET['NNj'], 'o', 'add-one: (NNj + v) − NNj'),
                                                   (DROP, 0.14, R.C_NET['NNjBC'], 's', 'drop-one: NNjBC − (NNjBC − v)')):
                vals = [diff[key[v_], gen, cat] for v_ in C.BC_VARS]
                ax.errorbar(np.arange(1, 5) + off, [x[0][i] for x in vals], yerr=[x[1][i] for x in vals],
                            fmt=marker, ms=8, color=color, mec=R.SURFACE, mew=1.5, elinewidth=1.4, capsize=0,
                            label=label)
            ax.axvline(0.5, color=R.RULE, lw=0.8)
            if r == 0:
                ax.set_title(f'{"gate jets" if cat == "gate jets" else "VBF-preselected jets"}, '
                             f'{GEN_LABEL[gen]}', fontsize=10, color=R.INK)
            if c == 0:
                ax.set_ylabel(ylabel, color=R.INK)
            ax.set_xticks(range(5))
            ax.set_xticklabels(xlab, fontsize=8.5, color=R.INK2)
            ax.grid(axis='y', color=R.RULE, lw=0.4, alpha=0.7)
            ax.tick_params(colors=R.INK2)
            for sp in ax.spines.values():
                sp.set_color(R.RULE)
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc='upper center', ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.995))
    fig.suptitle('Beam-connection variables on top of NNj: paired out-of-fold differences, event-weighted '
                 '(± max of bootstrap σ and fold-model s.e.)', fontsize=10.5, color=R.INK, y=0.95)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    out = os.path.join(C.FIGS, 'bc_ablation')
    for ext in ('pdf', 'png'):
        fig.savefig(f'{out}.{ext}', dpi=150, facecolor=R.SURFACE)
    print(f'wrote {out}.pdf and .png')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--paired-only', action='store_true', help='run only the paired summer26/nov26 check')
    ap.add_argument('--replicas', type=int, default=400)
    ap.add_argument('--seed', type=int, default=20261118)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    if args.paired_only:
        paired_summer26_nov26(rng, args.replicas)
        return

    # check: the summer26 NNj AUCs quoted in summer26.md (jet 0 only, unweighted)
    from sklearn.metrics import roc_auc_score
    for gen in GENS:
        (sig, _), (bkg, _) = (summer26_jets(s) for s in SAMPLES[gen])
        m_s, m_b = sig['slot'] == 0, bkg['slot'] == 0
        a = roc_auc_score(np.r_[np.ones(m_s.sum()), np.zeros(m_b.sum())], np.r_[sig['score'][m_s], bkg['score'][m_b]])
        print(f'check: summer26 NNj jet 0, {GEN_LABEL[gen]}, unweighted AUC {a:.4f}')

    res, unweighted, stats, reps_all, permodel = {}, {}, {}, {}, {}
    for gen in GENS:
        sig_name, bkg_name = SAMPLES[gen]
        # summer26 NNj: VBF presel. only, its own replicas
        (sig, n_s), (bkg, n_b) = summer26_jets(sig_name), summer26_jets(bkg_name)
        p = prepare(sig, bkg, True)
        r = np.array([metrics(p, rng.poisson(1., n_s).astype(float), rng.poisson(1., n_b).astype(float))
                      for _ in range(args.replicas)])
        res['summer26', 'NNj', gen, 'VBF presel.'] = (metrics(p), r.std(0, ddof=1))
        unweighted['summer26', 'NNj', gen, 'VBF presel.'] = metrics(prepare(sig, bkg, False))
        stats['summer26', 'NNj', gen, 'VBF presel.'] = (len(sig['score']), len(bkg['score']))
        # nov26: every net, shared replicas
        prep = {}
        for net in C.NETS:
            for cat in CATS:
                (sig, n_s), (bkg, n_b) = (nov26_jets(net, s, cat) for s in (sig_name, bkg_name))
                prep[net, cat] = prepare(sig, bkg, True)
                if net in MAIN:
                    unweighted['nov26', net, gen, cat] = metrics(prepare(sig, bkg, False))
                    stats['nov26', net, gen, cat] = (len(sig['score']), len(bkg['score']))
                permodel[net, gen, cat] = np.array([
                    metrics(prepare(C.take(sig, sig['model'] == k), C.take(bkg, bkg['model'] == k), True))
                    for k in range(C.N_FOLDS)])
        reps = {key: [] for key in prep}
        for _ in range(args.replicas):
            m_s = rng.poisson(1., n_s).astype(np.float64)
            m_b = rng.poisson(1., n_b).astype(np.float64)
            for key, p in prep.items():
                reps[key].append(metrics(p, m_s, m_b))
        for (net, cat), p in prep.items():
            reps_all[net, gen, cat] = np.array(reps[net, cat])
            res['nov26', net, gen, cat] = (metrics(p), reps_all[net, gen, cat].std(0, ddof=1))
        print(f'{gen}: done ({time.time() - t0:.0f} s)')

    # ── main table ────────────────────────────────────────────────────────────
    for cat in CATS:
        print(f'\n### {cat}, event-weighted ({args.replicas} bootstrap replicas)\n')
        print('| Tested on | Campaign | Net | ' + ' | '.join(METRICS) + ' |')
        print('|---|---|---|---|---|---|')
        for gen in GENS:
            for campaign, net in (('summer26', 'NNj'), ('nov26', 'NNj'), ('nov26', 'NNjBC')):
                if (campaign, net, gen, cat) in res:
                    v, s = res[campaign, net, gen, cat]
                    print(f'| {GEN_LABEL[gen]} | {campaign} | {net} | ' + ' | '.join(fmt(v[i], s[i]) for i in range(3)) + ' |')
        print('\nunweighted: ' + '; '.join(f'{c} {n} {GEN_LABEL[g]} ' + ' '.join(f'{x:.4f}' for x in v)
                                           for (c, n, g, k), v in unweighted.items() if k == cat))

    # ── ablation ──────────────────────────────────────────────────────────────
    diff = {}
    for pair in [ALL] + list(ADD.values()) + list(DROP.values()):
        for gen in GENS:
            for cat in CATS:
                a, b = (pair[0], gen, cat), (pair[1], gen, cat)
                d = res[('nov26',) + a][0] - res[('nov26',) + b][0]
                sb = (reps_all[a] - reps_all[b]).std(0, ddof=1)
                pm = permodel[a] - permodel[b]
                se = pm.std(0, ddof=1) / np.sqrt(len(pm))
                diff[pair, gen, cat] = (d, np.maximum(sb, se))
    for cat in CATS:
        print(f'\n### Ablation, {cat}: paired differences ± max(bootstrap σ, fold-model s.e.), event-weighted\n')
        print('| Tested on | Variable | Add-one ΔAUC | Drop-one ΔAUC | Add-one Δε_S(1 %) | Drop-one Δε_S(1 %) '
              '| Add-one Δε_B(50 %) | Drop-one Δε_B(50 %) |')
        print('|---|---|---|---|---|---|---|---|')
        for gen in GENS:
            v, s = diff[ALL, gen, cat]
            print(f'| {GEN_LABEL[gen]} | all four (NNjBC − NNj) | ' + ' | '.join(
                f'{v[i]:+.4f} ± {s[i]:.4f} | —' for i in range(3)) + ' |')
            for var in C.BC_VARS:
                (va, sa), (vd, sd) = diff[ADD[var], gen, cat], diff[DROP[var], gen, cat]
                print(f'| {GEN_LABEL[gen]} | bc_{var} | ' + ' | '.join(
                    f'{va[i]:+.4f} ± {sa[i]:.4f} | {vd[i]:+.4f} ± {sd[i]:.4f}' for i in range(3)) + ' |')
    print('\nper-fold-model spread of NNj and NNjBC (event-weighted AUC, own inference jets)')
    for net in MAIN:
        for gen in GENS:
            for cat in CATS:
                a = permodel[net, gen, cat][:, 0]
                print(f'  {net:6s} {GEN_LABEL[gen]:10s} {cat:11s} ' + ' '.join(f'{x:.4f}' for x in a)
                      + f'  (rms {a.std():.4f})')
    print('\njets (VBF / QCD)')
    for (campaign, net, gen, cat), (ns, nb) in stats.items():
        print(f'  {campaign:8s} {net:6s} {GEN_LABEL[gen]:10s} {cat:11s} {ns:8d} / {nb:7d}')

    paired_summer26_nov26(rng, args.replicas)
    make_roc_figure()
    make_ablation_figure(diff)
    print(f'({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
