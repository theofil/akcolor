#!/usr/bin/env python3
"""Beam-connection (bc) distributions in dijet-gate events: VBF Z vs QCD Z+jj, Herwig vs MG5+Pythia.

    python3 nov26/bc_plots.py                       # friends/nov26/*.h5 -> figs/nov26/bc_distributions.{pdf,png}
    python3 nov26/bc_plots.py --folder fast_tmp_nov26 --out <path stem>

Both stored jets of every gate event enter, weighted by kWeight, and every curve is normalised
to unit area. Values beyond the plotted range are added to the last bin. bc_dphi is shown as
|bc_dphi|; jets with an empty region (bc_dphi = -99) are left out of that panel only.
"""
import argparse
import os

import h5py
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # makefriends/C

C_VBF, C_QCD = '#2a78d6', '#eb6834'                                  # categorical slots 1-2
INK, INK2, RULE, SURFACE = '#0b0b0b', '#52514e', '#c3c2b7', '#fcfcfb'
SAMPLES = [  # (file stem, label, colour, line style): colour = process, line style = generator
    ('VBFZ_herwig', 'VBF Z, Herwig', C_VBF, '-'),
    ('VBFZ_mg5_pythia', 'VBF Z, MG5+Pythia', C_VBF, '--'),
    ('QCDZjj_herwig', 'QCD Z+jj, Herwig', C_QCD, '-'),
    ('QCDZjj_mg5_pythia', 'QCD Z+jj, MG5+Pythia', C_QCD, '--'),
]
VARS = [  # (dataset, axis label, lo, hi, bins)
    ('jetBCin', 'bc_in [GeV]', 0.0, 40.0, 40),
    ('jetBCout', 'bc_out [GeV]', 0.0, 150.0, 50),
    ('jetBCmet', 'bc_met [GeV]', 0.0, 80.0, 40),
    ('jetBCdphi', '|bc_dphi| [rad]', 0.0, np.pi, 32),
]


def load(path):
    """Per-jet bc values and kWeight for the stored jets of gate events."""
    with h5py.File(path, 'r') as h:
        gate, n_jets, kw = h['gate'][:], h['nJets'][:], h['kWeight'][:]
        valid = gate[:, None] & (np.arange(2)[None, :] < n_jets[:, None])
        w = np.broadcast_to(kw[:, None], valid.shape)[valid]
        return {v: h[v][:][valid] for v, *_ in VARS}, w, int(gate.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--folder', default=os.path.join(BASE, 'friends', 'nov26'))
    ap.add_argument('--out', default=os.path.join(BASE, 'figs', 'nov26', 'bc_distributions'))
    args = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.8), facecolor=SURFACE)
    for stem, label, color, ls in SAMPLES:
        vals, w, n_gate = load(os.path.join(args.folder, f'{stem}.h5'))
        print(f'{stem:18s} gate events {n_gate}, jets {len(w)}')
        for ax, (v, xlabel, lo, hi, nb) in zip(axes.flat, VARS):
            x, ww = vals[v], w
            if v == 'jetBCdphi':
                keep = x != -99
                x, ww = np.abs(x[keep]), ww[keep]
            edges = np.linspace(lo, hi, nb + 1)
            h, _ = np.histogram(np.clip(x, lo, np.nextafter(hi, lo)), bins=edges, weights=ww)
            ax.stairs(h / (h.sum() * np.diff(edges)), edges, color=color, ls=ls, lw=1.8,
                      label=label if v == 'jetBCin' else None)
    for ax, (v, xlabel, lo, hi, nb) in zip(axes.flat, VARS):
        ax.set_xlim(lo, hi)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(xlabel + ('' if v == 'jetBCdphi' else '  (last bin includes overflow)'), color=INK)
        ax.set_ylabel('normalised density', color=INK)
        ax.tick_params(colors=INK2, labelsize=9)
        ax.grid(color=RULE, lw=0.4, alpha=0.6)
        for sp in ax.spines.values():
            sp.set_color(RULE)
    axes[1, 1].set_xticks(np.arange(5) * np.pi / 4, ['0', 'π/4', 'π/2', '3π/4', 'π'])
    fig.legend(loc='upper center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle('nov26 beam connection: both leading jets of dijet-gate events, kWeight-weighted, '
                 'unit area', y=0.995, fontsize=10.5, color=INK)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.08, hspace=0.32, wspace=0.22)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{args.out}.{ext}', dpi=150, facecolor=SURFACE)
    print(f'wrote {args.out}.pdf and .png')


if __name__ == '__main__':
    main()
