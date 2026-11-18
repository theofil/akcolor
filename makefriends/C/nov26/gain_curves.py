#!/usr/bin/env python3
"""What the constituents and the beam connection add beyond the jet kinematics: signal-efficiency
ratios at equal background efficiency (lxplus; reads the stored scores, runs no model).

    python3 nov26/gain_curves.py                         # tables + figs/nov26/gain_curves.{pdf,png}
    python3 nov26/gain_curves.py --replicas 50           # quicker, rougher uncertainties
    python3 nov26/gain_curves.py --nets NNj NNjBC --ref NNj --no-fig   # any subset of nets

Nets (out-of-fold h5 columns, save_scores.py), a 2 x 2 of constituents x beam connection (bc) on
the same jets, folds and training recipe:
  Kj     jet |eta|, m, pT only (NNj without the constituent branch)
  KjBC   Kj + the four bc values
  NNj    Kj + the constituents
  NNjBC  NNj + the four bc values
Figure of merit: R(eps_B) = eps_S(net) / eps_S(ref) at equal background efficiency eps_B. Both nets
then keep the same background, so at that working point R is the ratio of S/sqrt(B), of S/B and of
S/(delta B): it holds whether the analysis is limited by statistics or by a background systematic
(for S << B; S/B ~ 3 % after the VBF preselection). R^2 is the luminosity equivalent.
Printed: the ratios to Kj (NNj / Kj: constituents, KjBC / Kj: bc, NNjBC / Kj: both) and on top of
each other (NNjBC / NNj, NNjBC / KjBC), at eps_B = 10 %, 1 %, 0.3 %; the maximum of
eps_S / sqrt(eps_B) of each net (the "performance" of the summer26 table: the factor by which a cut
on the score multiplies S/sqrt(B) of the category); a population table.
Categories (VBF Z = signal, QCD Z+jj = background, event-weighted with kWeight):
  gate jets            both jets of gate events (roles A + B)
  VBF presel. jets     both jets of gate events with mjj > 400 GeV and |dYjj| > 2.5 (as roc_scores.jets)
  VBF presel. events   one entry per VBF-preselected event: logit(s0) + logit(s1) of its two jets,
                       the same combination for every net
Curves stop where fewer than 100 raw background entries remain above the cut.
Uncertainty: paired Poisson bootstrap. In every replica each event gets one Poisson(1) multiplicity
(independently for the signal and the background sample), shared by its two jets and by all nets.
The tables quote max(bootstrap sigma, standard error of the five per-fold-model values), as
roc_summary.py does; the second includes the training randomness. Both jets of an event belong to
the same sub-fold, so the event category has per-fold-model values too.
Population table: AUC and eps_S(1 %) of |eta|, pT and m alone and of the nets, in gate jets, gate jets
with |dYjj| > 2.5, gate jets with mjj > 400 GeV and VBF presel. jets, plus the AUC at fixed jet
(|eta|, pT): the signal-weighted mean of the AUCs in 15 x 20 (|eta|, log pT) bins with >= 50 raw
jets of each class (a cross-check that removes the jet kinematics).
"""
import argparse
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import roc_scores as R  # noqa: E402
from bootstrap_scores import event_index, trapezoid  # noqa: E402

NETS = ('Kj', 'KjBC', 'NNj', 'NNjBC')
REF = 'Kj'
ON_TOP = (('NNjBC', 'NNj'), ('NNjBC', 'KjBC'))
CATS = ('gate jets', 'VBF presel. jets', 'VBF presel. events')
POP_CATS = ('gate jets', 'gate, |dYjj| > 2.5', 'gate, mjj > 400', 'VBF presel. jets')
TABLE_EPS_B = (0.10, 0.01, 0.003)
MIN_RAW_B = 100
N_GRID = 80
KIN_VARS = (('|eta|', 'eta'), ('pT', 'pt'), ('m', 'm'))
ETA_BINS, PT_BINS = np.linspace(0., 3., 16), np.geomspace(30., 3000., 21)
MIN_BIN = 50
CAT_TITLE = {'gate jets': 'gate jets',
             'VBF presel. jets': 'VBF-preselected jets\n(mjj > 400 GeV, |ΔY_jj| > 2.5)',
             'VBF presel. events': 'VBF-preselected events, both jets\n(logit s₀ + logit s₁)'}
NET_LABEL = {'Kj': 'Kj: jet |η|, m, pT', 'KjBC': 'KjBC: Kj + beam connection',
             'NNj': 'NNj: Kj + constituents', 'NNjBC': 'NNjBC: NNj + beam connection'}


# ── data ──────────────────────────────────────────────────────────────────────

def load(sample, nets):
    """Per-jet arrays of `sample` in the jet order of roc_scores.jets: the out-of-fold score of every
    net in `nets` that has a column, kWeight, role, VBF flag, fold model, event row, slot, |eta|, pT,
    m, and the event's mjj and |dYjj|. None if no net has scores."""
    d, scores = None, {}
    for net in nets:
        j = R.jets(sample, net)
        if j is None:
            continue
        scores[net] = j['score'].astype(np.float64)
        if d is None:
            d = {k: j[k] for k in ('kw', 'role', 'vbf', 'model')}
    if d is None:
        return None
    rows, n_ev = event_index(sample)
    with h5py.File(C.h5_path(sample), 'r') as f:
        n_jets = f['nJets'][:]
        slots = np.r_[np.zeros((n_jets > 0).sum(), np.int64), np.ones((n_jets > 1).sum(), np.int64)]
        for name, key in (('eta', 'jetEta'), ('pt', 'jetPt'), ('m', 'jetM')):
            d[name] = f[key][:][rows, slots].astype(np.float64)
        mjj, dy = f['mjj'][:][rows].astype(np.float64), np.abs(f['dYjj'][:][rows]).astype(np.float64)
    d.update(eta=np.abs(d['eta']), mjj=mjj, dy=dy, score=scores, ev=rows, slot=slots, n_ev=n_ev,
             kw=d['kw'].astype(np.float64))
    return d


def jet_mask(d, cat):
    gate = d['role'] != C.ROLE_C
    return {'gate jets': gate, 'VBF presel. jets': gate & d['vbf'],
            'gate, |dYjj| > 2.5': gate & (d['dy'] > 2.5), 'gate, mjj > 400': gate & (d['mjj'] > 400)}[cat]


def logit(p):
    p = np.clip(p, 1e-7, 1 - 1e-7)
    return np.log(p / (1 - p))


def entries(d, cat, x):
    """(value, weight, event row, fold model) of the entries of category `cat`; x is a per-jet array
    (a net's score or a kinematic variable)."""
    if cat != 'VBF presel. events':
        m = jet_mask(d, cat)
        return x[m], d['kw'][m], d['ev'][m], d['model'][m]
    a = (d['role'] == C.ROLE_A) & d['vbf']                      # slot 0 of VBF-preselected events
    b = (d['role'] == C.ROLE_B) & d['vbf']                      # slot 1, the same events in the same order
    assert np.array_equal(d['ev'][a], d['ev'][b]) and np.array_equal(d['model'][a], d['model'][b])
    return logit(x[a]) + logit(x[b]), d['kw'][a], d['ev'][a], d['model'][a]


# ── ROC machinery: sort once, reweight per replica ────────────────────────────

def prepare(sig, bkg):
    """Signal and background entries sorted once by descending value, with the last index of each
    tie group (the roc_summary.prepare pattern, plus the fold model)."""
    x = np.r_[sig[0], bkg[0]]
    y = np.r_[np.ones(len(sig[0]), bool), np.zeros(len(bkg[0]), bool)]
    o = np.argsort(-x, kind='stable')
    last = np.r_[np.nonzero(np.diff(x[o]))[0], len(x) - 1]
    return dict(y=y[o], w=np.r_[sig[1], bkg[1]][o], ev=np.r_[sig[2], bkg[2]][o],
                model=np.r_[sig[3], bkg[3]][o], last=last)


def curve(p, w=None):
    """(eps_S, eps_B, raw background count) above each cut, starting at the tightest."""
    w = p['w'] if w is None else w
    cs = np.cumsum(np.where(p['y'], w, 0.))[p['last']]
    cb = np.cumsum(np.where(p['y'], 0., w))[p['last']]
    nb = np.cumsum(~p['y'])[p['last']]
    return np.r_[0., cs / cs[-1]], np.r_[0., cb / cb[-1]], np.r_[0, nb]


def auc_eff(p, w=None):
    es, eb, _ = curve(p, w)
    return trapezoid(es, eb), float(np.interp(0.01, eb, es))


def max_performance(p):
    """Maximum of eps_S / sqrt(eps_B) over the cuts that keep >= MIN_RAW_B raw background entries,
    and its eps_B."""
    es, eb, nb = curve(p)
    ok = nb >= MIN_RAW_B
    perf = np.where(ok, es / np.sqrt(np.where(ok, eb, 1.)), 0.)
    i = int(np.argmax(perf))
    return perf[i], eb[i]


def auc_w(xs, ws, xb, wb):
    """Weighted AUC of signal values xs against background values xb, ties handled as in curve()."""
    z = np.zeros
    return auc_eff(prepare((xs, ws, z(len(xs), int), z(len(xs), int)), (xb, wb, z(len(xb), int), z(len(xb), int))))[0]


def oriented_auc_eff(sig, bkg):
    """(AUC, eps_S(1 %), '>' or '<'): the variable is flipped if lower values are more VBF-like."""
    a, e = auc_eff(prepare(sig, bkg))
    if a >= 0.5:
        return a, e, '>'
    a, e = auc_eff(prepare((-sig[0],) + sig[1:], (-bkg[0],) + bkg[1:]))
    return a, e, '<'


def kin_bins(d, cat):
    """Jets of `cat` ordered by their (|eta|, pT) bin: (mask, order, bin edges into the ordered jets)."""
    m = jet_mask(d, cat)
    ie = np.clip(np.digitize(d['eta'][m], ETA_BINS) - 1, 0, len(ETA_BINS) - 2)
    ip = np.clip(np.digitize(d['pt'][m], PT_BINS) - 1, 0, len(PT_BINS) - 2)
    k = ie * (len(PT_BINS) - 1) + ip
    o = np.argsort(k, kind='stable')
    return m, o, np.searchsorted(k[o], np.arange((len(ETA_BINS) - 1) * (len(PT_BINS) - 1) + 1))


def conditional_auc(ds, db, bs, bb, xs, xb):
    """Signal-weighted mean AUC in the (|eta|, pT) bins with >= MIN_BIN raw jets of each class, and the
    fraction of the signal weight in those bins. bs, bb: kin_bins() of the signal and background."""
    (ms, os_, es_), (mb, ob, eb_) = bs, bb
    vs, ws = xs[ms][os_], ds['kw'][ms][os_]
    vb, wb = xb[mb][ob], db['kw'][mb][ob]
    num = cov = 0.
    for k in range(len(es_) - 1):
        s, b = slice(es_[k], es_[k + 1]), slice(eb_[k], eb_[k + 1])
        if es_[k + 1] - es_[k] < MIN_BIN or eb_[k + 1] - eb_[k] < MIN_BIN:
            continue
        num += ws[s].sum() * auc_w(vs[s], ws[s], vb[b], wb[b])
        cov += ws[s].sum()
    return num / cov, cov / ws.sum()


# ── analysis ──────────────────────────────────────────────────────────────────

def analyse(ds, db, nets, replicas, rng):
    """Per category and net: eps_S on an eps_B grid (figure) and at TABLE_EPS_B (tables), nominal and
    for every bootstrap replica; per-fold-model values at TABLE_EPS_B; the maximum performance.
    Plus the population table."""
    preps = {(cat, net): prepare(entries(ds, cat, ds['score'][net]), entries(db, cat, db['score'][net]))
             for cat in CATS for net in nets}
    out = {}
    for cat in CATS:
        r = dict(es={}, rep={}, fold={}, perf={})
        floor = 0.
        for net in nets:
            es, eb, nb = curve(preps[cat, net])
            floor = max(floor, eb[np.searchsorted(nb, MIN_RAW_B)])
        r['grid'] = np.geomspace(floor, 1., N_GRID)
        r['table'] = np.array([e if e >= floor else np.nan for e in TABLE_EPS_B])
        r['x'] = np.r_[r['grid'], np.nan_to_num(r['table'], nan=1.)]    # grid points, then table points
        for net in nets:
            p = preps[cat, net]
            es, eb, _ = curve(p)
            r['es'][net] = np.interp(r['x'], eb, es)
            r['fold'][net] = np.array([np.interp(r['table'], *curve(p, p['w'] * (p['model'] == k))[1::-1])
                                       for k in range(C.N_FOLDS)])
            r['rep'][net] = np.empty((replicas, len(r['x'])))
            r['perf'][net] = max_performance(p)
        out[cat] = r
    for i in range(replicas):
        m_s = rng.poisson(1., ds['n_ev']).astype(np.float64)
        m_b = rng.poisson(1., db['n_ev']).astype(np.float64)
        for (cat, net), p in preps.items():
            es, eb, _ = curve(p, p['w'] * np.where(p['y'], m_s[p['ev']], m_b[p['ev']]))
            out[cat]['rep'][net][i] = np.interp(out[cat]['x'], eb, es)
    for cat in CATS:                                                  # eps_S = 1 at eps_B = 1 (no cut)
        assert all(abs(out[cat]['es'][n][N_GRID - 1] - 1) < 1e-9 for n in nets)

    pop, signs = {}, {}
    for cat in POP_CATS:                                  # POP_CATS[0] = gate jets fixes each direction
        bs, bb = kin_bins(ds, cat), kin_bins(db, cat)
        inputs = [(name, ds[key], db[key]) for name, key in KIN_VARS]
        inputs += [(net, ds['score'][net], db['score'][net]) for net in nets]
        for name, xs, xb in inputs:
            if name not in signs:
                signs[name] = oriented_auc_eff(entries(ds, cat, xs), entries(db, cat, xb))[2]
            f = 1. if signs[name] == '>' else -1.
            a, e = auc_eff(prepare(entries(ds, cat, f * xs), entries(db, cat, f * xb)))
            ca, cov = conditional_auc(ds, db, bs, bb, f * xs, f * xb)
            pop[cat, name] = (a, e, signs[name] if name not in nets else '', ca, cov)
    return out, pop


def ratio_stats(r, net, ref):
    """R at TABLE_EPS_B (NaN below the curve floor), and max(bootstrap sigma, fold standard error)."""
    t = slice(N_GRID, None)
    nom = r['es'][net][t] / r['es'][ref][t]
    rep = r['rep'][net][:, t] / r['rep'][ref][:, t]
    fold = r['fold'][net] / r['fold'][ref]
    se = fold.std(axis=0, ddof=1) / np.sqrt(C.N_FOLDS)
    bad = np.isnan(r['table'])
    return np.where(bad, np.nan, nom), np.where(bad, np.nan, np.maximum(rep.std(axis=0), se))


def print_tables(res, pop, nets, ref):
    gens = [g for g in R.GENERATORS if any((g, c) in res for c in CATS)]
    pairs = [(n, ref) for n in nets if n != ref] + [p for p in ON_TOP if p[0] in nets and p[1] in nets and p[1] != ref]
    for cat in CATS:
        print(f'\n### eps_S ratio at equal eps_B: {cat}\n')
        head = ' | '.join(f'{g} {e:.1%}' for g in gens for e in TABLE_EPS_B)
        print(f'| Ratio | {head} |')
        print('|---|' + '---|' * (len(gens) * len(TABLE_EPS_B)))
        for net, base in pairs:
            cells = []
            for g in gens:
                nom, sig = ratio_stats(res[g, cat], net, base)
                cells += ['—' if np.isnan(v) else f'{v:.3f} ± {s:.3f}' for v, s in zip(nom, sig)]
            print(f'| {net} / {base} | ' + ' | '.join(cells) + ' |')
        floors = ', '.join(f'{g} {res[g, cat]["grid"][0]:.2%}' for g in gens)
        print(f'\n(curves start at eps_B = {floors}: {MIN_RAW_B} raw background entries; — = below that)')

    print('\n### Maximum of eps_S / sqrt(eps_B) (performance) and its eps_B\n')
    print('| Category | ' + ' | '.join(f'{n}, {g}' for g in gens for n in nets) + ' |')
    print('|---|' + '---|' * (len(gens) * len(nets)))
    for cat in CATS:
        cells = [f'{res[g, cat]["perf"][n][0]:.3f} ({res[g, cat]["perf"][n][1]:.0%})' for g in gens for n in nets]
        print(f'| {cat} | ' + ' | '.join(cells) + ' |')

    for g in gens:
        print(f'\n### Population dependence, {g}: AUC / eps_S(1 %) / AUC at fixed (|eta|, pT)\n')
        print('| Input | ' + ' | '.join(POP_CATS) + ' |')
        print('|---|' + '---|' * len(POP_CATS))
        for name in [v for v, _ in KIN_VARS] + list(nets):
            cells = []
            for cat in POP_CATS:
                a, e, sign, ca, cov = pop[g][cat, name]
                cells.append(f'{a:.4f} / {e:.4f} / {ca:.4f}')
            sign = pop[g][POP_CATS[0], name][2]
            label = f'{name} alone ({sign})' if sign else name
            print(f'| {label} | ' + ' | '.join(cells) + ' |')
        cov = min(pop[g][cat, n][4] for cat in POP_CATS for n in nets)
        print(f'\n(fixed (|eta|, pT): 15 x 20 bins with >= {MIN_BIN} raw jets per class, covering >= {cov:.1%} '
              f'of the signal weight; (>) / (<): larger / smaller values counted as VBF-like, fixed on gate jets '
              f'for every column, so a value below 0.5 means the opposite direction there)')


# ── figure ────────────────────────────────────────────────────────────────────

def make_figure(res, nets, ref):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    gens = [g for g in R.GENERATORS if any((g, c) in res for c in CATS)]
    fig, axes = plt.subplots(2, len(CATS), sharex='col', sharey='row', squeeze=False, facecolor=R.SURFACE,
                             figsize=(5.0 * len(CATS), 7.6))
    g_ = slice(0, N_GRID)
    for col, cat in enumerate(CATS):
        top, bot = axes[0, col], axes[1, col]
        for g in gens:
            r = res[g, cat]
            grid = r['grid']
            for net in nets:
                top.plot(grid, r['es'][net][g_] / np.sqrt(grid), color=R.C_NET[net], lw=1.6, ls=R.LS[g])
            for net in nets:
                if net == ref:
                    continue
                ratio = r['es'][net][g_] / r['es'][ref][g_]
                sd = (r['rep'][net][:, g_] / r['rep'][ref][:, g_]).std(axis=0)
                bot.fill_between(grid, ratio - sd, ratio + sd, color=R.C_NET[net], alpha=0.22, lw=0)
                bot.plot(grid, ratio, color=R.C_NET[net], lw=1.6, ls=R.LS[g])
        for ax in (top, bot):
            ax.axhline(1, color=R.INK2, lw=0.9, ls=':')
            ax.set_xscale('log')
            ax.set_xlim(min(res[g, cat]['grid'][0] for g in gens), 1)
        top.set_title(CAT_TITLE[cat], fontsize=10, color=R.INK)
        bot.set_xlabel('background efficiency ε_B (QCD Z+jj)', color=R.INK)
    # direct labels: above the ratio curves of the first generator, at a mid-range eps_B where the
    # curves are separated (the tight end is noisy, the loose end converges to 1)
    for col, cat in enumerate(CATS):
        bot = axes[1, col]
        r = res[gens[0], cat]
        x = np.sqrt(r['grid'][0] * 0.2)
        ys = sorted((np.interp(x, r['grid'], r['es'][n][g_] / r['es'][ref][g_]), n) for n in nets if n != ref)
        lo, hi = bot.get_ylim()
        sep, last = 0.075 * (hi - lo), -np.inf
        for y, n in ys:
            y = max(y + 0.02 * (hi - lo), last + sep)
            bot.text(x, y, f'{n} / {ref}', color=R.INK, fontsize=8, ha='center', va='bottom',
                     bbox=dict(boxstyle='round,pad=0.15', fc=R.SURFACE, ec='none', alpha=0.85))
            last = y
    axes[0, 0].set_ylabel('ε_S / √ε_B  (gain in S/√B from the cut)', color=R.INK)
    axes[1, 0].set_ylabel(f'ε_S(net) / ε_S({ref}) at equal ε_B', color=R.INK)
    handles = [Line2D([], [], color=R.C_NET[n], lw=1.8, label=NET_LABEL[n]) for n in nets]
    handles += [Line2D([], [], color=R.INK, lw=1.4, ls=R.LS[g], label=g) for g in gens]
    axes[0, 0].legend(handles=handles, frameon=False, fontsize=7.5, loc='upper left')
    axes[1, -1].text(0.97, 0.96, 'bands: ±1σ, paired bootstrap', transform=axes[1, -1].transAxes, fontsize=7.5,
                     color=R.INK2, ha='right', va='top')
    for a in axes.flat:
        a.grid(color=R.RULE, lw=0.4, alpha=0.7, which='both')
        a.tick_params(colors=R.INK2)
        for sp in a.spines.values():
            sp.set_color(R.RULE)
    fig.suptitle(f'nov26: signal efficiency at equal background efficiency, relative to {ref}', fontsize=11,
                 color=R.INK)
    fig.tight_layout()
    os.makedirs(C.FIGS, exist_ok=True)
    out = os.path.join(C.FIGS, 'gain_curves')
    for ext in ('pdf', 'png'):
        fig.savefig(f'{out}.{ext}', dpi=150, facecolor=R.SURFACE)
    print(f'\nwrote {out}.pdf and .png')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--nets', nargs='+', default=list(NETS), choices=list(C.NETS))
    ap.add_argument('--ref', default=REF, choices=list(C.NETS), help='denominator of the ratio curves')
    ap.add_argument('--replicas', type=int, default=200)
    ap.add_argument('--no-fig', action='store_true')
    args = ap.parse_args()
    nets = list(dict.fromkeys([args.ref] + args.nets))
    rng = np.random.default_rng(C.SEED)
    res, pop = {}, {}
    for g, (sig, bkg) in R.GENERATORS.items():
        ds, db = load(sig, nets), load(bkg, nets)
        if ds is None or db is None:
            print(f'{g}: no scores for {", ".join(nets)}; skipped')
            continue
        have = [n for n in nets if n in ds['score'] and n in db['score']]
        if args.ref not in have:
            sys.exit(f'{g}: no scores for the reference {args.ref} (run save_scores.py --net {args.ref})')
        missing = sorted(set(nets) - set(have))
        if missing:
            print(f'{g}: no scores for {", ".join(missing)}; left out')
        nets = have
        out, pop[g] = analyse(ds, db, nets, args.replicas, rng)
        for cat, r in out.items():
            res[g, cat] = r
        print(f'{g}: done ({args.replicas} replicas)', flush=True)
    print_tables(res, pop, nets, args.ref)
    if not args.no_fig:
        make_figure(res, nets, args.ref)


if __name__ == '__main__':
    main()
