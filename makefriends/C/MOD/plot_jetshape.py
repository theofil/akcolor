#!/usr/bin/env python3
"""2D jet-shape maps analogous to the `jetShape`/`jetShapeRaw` TH2Poly histograms
written by makefriends.cpp -- the average per-constituent energy density around
the jet axis in (Delta-y, Delta-phi), summed over both jets of every dijet pair.

makefriends.cpp fills these directly from FastJet constituents at clustering time
and writes them as ROOT TH2Poly (annular-sector bins) into the .friend.root file;
root_to_hdf5.py never carries them into the .h5 conversion, and MOD writes .h5
only (no ROOT step at all -- see MOD/README.md), so there is no ROOT file to pull
them from here. This reproduces the same annular-sector geometry, pull weight and
color-map convention as summer26/plots.py's raw TH2Poly->Draw('COLZ') (128 phi
wedges x 60 rho rings out to r=0.5, linear z-scale, ROOT-style rainbow palette),
built directly from what's already stored in the friend h5 files (jcsDEta/jcsDPhi/
jcsW/jcsPt) instead of reprocessing raw Zenodo data. Two disclosed approximations
relative to the original:
  1. Position/weight use pseudorapidity Delta-eta (jcsDEta, already sign-flipped
     per jet exactly like the original's Delta-y flip) instead of rapidity
     Delta-y -- negligible difference for relativistic PF candidates, and the
     sign convention is identical either way (a jet's eta and y always share sign).
  2. jcsPt/jcsDEta/jcsDPhi/jcsW are already truncated to the top NC_MAX=80
     constituents by pt (see make_mod_friends.py); makefriends.cpp's jetShape
     uses the untruncated constituent list. The dropped tail is the softest
     constituents (pt-sorted), which carry negligible weight in either
     jetShape's or jetShapeRaw's pt-weighted sum.

jetShape   weight = jcsW                     (pt fraction times radial factor,
                                                same formula as makefriends.cpp's
                                                w = (pt_i/pt_jet)^PV_A * r^PV_B,
                                                PV_A=PV_B=1, PV_C=0)
jetShapeRaw weight = jcsPt / jetPt[jet slot]  (pt fraction only)

Note: under the linear z-scale, both maps show a faint plus/cross-shaped
enhancement along Delta-eta=0 and Delta-phi=0, in both DATA and MC. Confirmed
(via a plain rectangular sanity histogram straight from jcsDEta/jcsDPhi/jcsW,
independent of this script's polar-binning code) to be present in the raw
data itself, not a plotting artifact. Most likely a genuine CMS
detector-reconstruction effect (calorimeter/tracker projective geometry) --
MOD's jets are real 2011 data + Pythia6 GEANT-simulated SIM, unlike
summer26's pure generator-truth jets (no detector simulation), which is why
summer26's jetShape looks smoothly radially symmetric and MOD's doesn't.

Processed one sample (one h5 file) at a time to keep peak memory bounded -- the
full combined-MC constituent arrays are ~40GB if loaded all at once, which is
exactly what OOM-killed make_mod_friends.py before it was fixed to stream.

Usage:
    python3 plot_jetshape.py
"""

import os
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FRIENDS_DIR = os.path.join(HERE, "friends")
FIGS_DIR = os.path.join(HERE, "figs")

DATA_SAMPLE = "CMS2011A_data"
QCD_SAMPLES = ["QCD_170-300", "QCD_300-470", "QCD_470-600", "QCD_600-800",
               "QCD_800-1000", "QCD_1000-1400", "QCD_1400-1800", "QCD_1800-inf"]

# Same annular-sector binning as makefriends.cpp's TH2Poly (128 phi wedges x
# 60 rho rings, r in [0, 0.5]); histogram2d silently drops entries with
# r>0.5, matching TH2Poly::Fill()'s behavior of not counting points that fall
# outside every defined polygon.
NPHI = 128
NRHO = 60
RHO_MAX = 0.5
PHI_BINS = np.linspace(0., 2. * np.pi, NPHI + 1)
RHO_BINS = np.linspace(0., RHO_MAX, NRHO + 1)


def accumulate_sample(sample):
    """Return (counts_shape, counts_raw, n_jets) for one sample, both jet slots
    combined, using only what's already in the friend h5 file."""
    path = os.path.join(FRIENDS_DIR, f"{sample}.h5")
    with h5py.File(path, "r") as f:
        jcsDEta = f["jcsDEta"][:]   # (N, 2, 80)
        jcsDPhi = f["jcsDPhi"][:]
        jcsPt = f["jcsPt"][:]
        jcsW = f["jcsW"][:]
        jetPt = f["jetPt"][:]       # (N, 2)

    n = jcsPt.shape[0]
    valid = jcsPt > 0  # zero-padded slots (beyond each jet's real constituent count)

    deta = jcsDEta[valid]
    dphi = jcsDPhi[valid]
    w_shape = jcsW[valid]

    jetPt_b = np.broadcast_to(jetPt[:, :, None], jcsPt.shape)
    w_raw = (jcsPt / jetPt_b)[valid]

    # (deta, dphi) -> polar (rho, theta) around the jet axis, same convention
    # as makefriends.cpp filling TH2Poly's Cartesian (dYs, dPhi) plane with
    # annular-sector bins parametrized by angle/radius.
    rho = np.sqrt(deta * deta + dphi * dphi)
    theta = np.arctan2(dphi, deta) % (2. * np.pi)

    counts_shape, _, _ = np.histogram2d(rho, theta, bins=[RHO_BINS, PHI_BINS], weights=w_shape)
    counts_raw, _, _ = np.histogram2d(rho, theta, bins=[RHO_BINS, PHI_BINS], weights=w_raw)

    del jcsDEta, jcsDPhi, jcsPt, jcsW, jetPt, deta, dphi, w_shape, w_raw, valid, jetPt_b, rho, theta
    return counts_shape, counts_raw, 2 * n  # 2 jets per event


def accumulate_group(samples):
    tot_shape = np.zeros((NRHO, NPHI))
    tot_raw = np.zeros_like(tot_shape)
    tot_n = 0
    for s in samples:
        print(f"  accumulating {s} ...")
        cs, cr, n = accumulate_sample(s)
        tot_shape += cs
        tot_raw += cr
        tot_n += n
    return tot_shape, tot_raw, tot_n


def draw_side_by_side(map_data, map_mc, title, outpath):
    # Bin-edge grid in Cartesian (x=Delta-eta, y=Delta-phi) space so pcolormesh
    # renders true annular sectors on ordinary linear axes -- same visual as
    # ROOT's TH2Poly->Draw('COLZ') on a square (-0.5,0.5)x(-0.5,0.5) canvas.
    R, THETA = np.meshgrid(RHO_BINS, PHI_BINS, indexing="ij")
    X = R * np.cos(THETA)
    Y = R * np.sin(THETA)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    vmax = max(map_data.max(), map_mc.max(), 1e-12)
    norm = mcolors.Normalize(vmin=0., vmax=vmax)  # linear z-scale, matching ROOT (no SetLogz)

    for ax, m, label in zip(axes, [map_data, map_mc], ["CMS2011A data", "QCD-MC (combined 8 bins)"]):
        im = ax.pcolormesh(X, Y, m, norm=norm, cmap="jet")
        ax.set_xlabel(r"$\Delta\eta_i$")
        ax.set_ylabel(r"$\Delta\phi_i$")
        ax.set_title(label, fontsize=12)
        ax.set_aspect("equal")
        ax.set_xlim(-RHO_MAX, RHO_MAX)
        ax.set_ylim(-RHO_MAX, RHO_MAX)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(title)
    plt.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print("  Saved", outpath)


def main():
    os.makedirs(FIGS_DIR, exist_ok=True)

    print(f"Accumulating {DATA_SAMPLE} ...")
    data_shape, data_raw, data_n = accumulate_group([DATA_SAMPLE])
    print(f"Accumulating combined QCD-MC ({len(QCD_SAMPLES)} samples) ...")
    mc_shape, mc_raw, mc_n = accumulate_group(QCD_SAMPLES)

    # normalize to average per-jet energy density (so DATA/MC are on a common scale
    # despite very different total statistics)
    data_shape_avg = data_shape / data_n
    mc_shape_avg = mc_shape / mc_n
    data_raw_avg = data_raw / data_n
    mc_raw_avg = mc_raw / mc_n

    print(f"  data: {data_n} jets, MC: {mc_n} jets")

    draw_side_by_side(data_shape_avg, mc_shape_avg,
                       "jetShape: avg. pull-weighted energy density per jet",
                       os.path.join(FIGS_DIR, "jetShape_data_vs_mc.pdf"))
    draw_side_by_side(data_raw_avg, mc_raw_avg,
                       "jetShapeRaw: avg. pT-fraction energy density per jet",
                       os.path.join(FIGS_DIR, "jetShapeRaw_data_vs_mc.pdf"))


if __name__ == "__main__":
    main()
