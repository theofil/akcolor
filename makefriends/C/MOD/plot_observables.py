#!/usr/bin/env python3
"""DATA-vs-combined-QCD-MC overlay plots for mjj, jet pT, jet eta, the pull angle
theta_21, and the leading-jet signed pull-vector angle (SPVA), from
MOD/friends/*.h5 files.

Style matches summer26/plots.py: MC drawn as a step histogram outline (same
apply_style/FIG_SIZE helpers), DATA drawn as black circle markers with
statistical error bars (ROOT marker style 20 equivalent) -- no standalone,
data-only plots are produced.

theta_21 is the classical jet-pull observable (Gallicchio & Schwartz, arXiv:1001.5027):
the angle between the SUBLEADING jet's own pull vector and the direction pointing from
the subleading jet toward the LEADING jet, in the (rapidity, phi) plane. It tests whether
the subleading jet's radiation is "pulled" toward its color partner.

All the inputs needed are already in the friend file:
  - jetPVA (raw/unsigned pull-vector angle, per jet slot, atan2(pv_phi, pv_y) around the
    jet's own axis)
  - dYjj, dPhijj  ( = y_lead - y_sub ,  deltaPhi(phi_lead, phi_sub) -- computed from the
    pt-ordered pair before any slot randomization, so this is exactly the direction
    FROM subleading TO leading jet)
  - leadJetIndex  (which slot, 0 or 1, holds the leading jet)

theta_21 = jetPVA[slot_sub] - atan2(dPhijj, dYjj), wrapped to (-pi, pi].

Usage:
    python3 plot_observables.py
"""

import os
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/afs/cern.ch/user/t/theofil/work/akcolor/analysis/Winter25/python/NN")
from style import apply_style, FIG_SIZE  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FRIENDS_DIR = os.path.join(HERE, "friends")
FIGS_DIR = os.path.join(HERE, "figs")

DATA_SAMPLE = "CMS2011A_data"
QCD_SAMPLES = ["QCD_170-300", "QCD_300-470", "QCD_470-600", "QCD_600-800",
               "QCD_800-1000", "QCD_1000-1400", "QCD_1400-1800", "QCD_1800-inf"]

NEEDED = ["mjj", "jetPt", "jetEta", "jetPVA", "jetSPVA", "dYjj", "dPhijj",
          "leadJetIndex", "kWeight"]

MC_COLOR = "tomato"
DATA_COLOR = "black"


def delta_phi(p1, p2):
    return np.arctan2(np.sin(p1 - p2), np.cos(p1 - p2))


def wrap_pi(a):
    return np.arctan2(np.sin(a), np.cos(a))


def load_observables(sample):
    """Read just the small per-event branches needed here (not the big [N,2,80]
    constituent arrays) and compute the derived theta_21 / leading-jet SPVA."""
    path = os.path.join(FRIENDS_DIR, f"{sample}.h5")
    with h5py.File(path, "r") as f:
        d = {k: f[k][:] for k in NEEDED}

    n = len(d["mjj"])
    idx = np.arange(n)
    lead = d["leadJetIndex"]
    sub_slot = 1 - lead

    dir_sub_to_lead = np.arctan2(d["dPhijj"], d["dYjj"])
    d["theta21"] = wrap_pi(d["jetPVA"][idx, sub_slot] - dir_sub_to_lead)
    d["spva_lead"] = d["jetSPVA"][idx, lead]
    d["n"] = n
    return d


def load_combined_qcd(samples):
    """Concatenate several samples' observables, keeping each sample's own kWeight
    (already cross-section-normalized per pT-hat bin, see MOD/README.md) so the
    combination is a properly weighted inclusive QCD-MC prediction."""
    parts = [load_observables(s) for s in samples]
    out = {}
    for key in ["mjj", "jetPt", "jetEta", "theta21", "spva_lead", "kWeight"]:
        out[key] = np.concatenate([p[key] for p in parts], axis=0)
    out["n"] = sum(p["n"] for p in parts)
    return out


def weighted_density_with_err(values, weights, bins):
    """Normalized histogram + Sumw2-style stat error (sqrt(sum w^2)/norm), the
    same convention summer26/plots.py uses for ROOT histograms with Sumw2."""
    counts, _ = np.histogram(values, bins=bins, weights=weights)
    sumw2, _ = np.histogram(values, bins=bins, weights=weights ** 2)
    bin_width = bins[1] - bins[0]
    norm = counts.sum() * bin_width
    if norm <= 0:
        z = np.zeros_like(counts, dtype=float)
        return z, z
    return counts / norm, np.sqrt(sumw2) / norm


def overlay_plot(data_vals, data_w, mc_vals, mc_w, bins, xlabel, ylabel, outpath,
                  log_y=False, xlim=None, ylim=None):
    centers = (bins[:-1] + bins[1:]) / 2

    mc_hist, mc_err = weighted_density_with_err(np.clip(mc_vals, bins[0], bins[-1]), mc_w, bins)
    data_hist, data_err = weighted_density_with_err(np.clip(data_vals, bins[0], bins[-1]), data_w, bins)

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.hist(centers, bins=bins, weights=mc_hist, histtype="step", linewidth=3,
            color=MC_COLOR, label=f"QCD-MC (combined 8 bins, N={len(mc_vals)})")
    ax.errorbar(centers, mc_hist, yerr=mc_err, fmt="none",
                color=MC_COLOR, capsize=3, capthick=1, elinewidth=1)
    ax.errorbar(centers, data_hist, yerr=data_err, fmt="o", color=DATA_COLOR,
                markersize=4, markerfacecolor=DATA_COLOR, capsize=3, capthick=1,
                elinewidth=1, label=f"CMS2011A data (N={len(data_vals)})")
    apply_style(ax, xlabel=xlabel, ylabel=ylabel, title="",
                xlim=xlim or (bins[0], bins[-1]), ylim=ylim, log_y=log_y, legend_loc="upper right")
    plt.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print("  Saved", outpath)


def make_comparison_plots():
    os.makedirs(FIGS_DIR, exist_ok=True)
    print(f"Loading {DATA_SAMPLE} ...")
    data = load_observables(DATA_SAMPLE)
    print(f"Loading and combining {len(QCD_SAMPLES)} QCD-MC samples ...")
    mc = load_combined_qcd(QCD_SAMPLES)
    print(f"  data: {data['n']} events, combined MC: {mc['n']} events")

    overlay_plot(data["mjj"], data["kWeight"], mc["mjj"], mc["kWeight"],
                 np.arange(0, 4000, 20),
                 r"$m_{jj}$ [GeV]", "Normalized density [1/GeV]",
                 os.path.join(FIGS_DIR, "mjj_data_vs_mc.pdf"), log_y=True, xlim=(0, 4000))

    data_pt = np.concatenate([data["jetPt"][:, 0], data["jetPt"][:, 1]])
    data_pt_w = np.concatenate([data["kWeight"], data["kWeight"]])
    mc_pt = np.concatenate([mc["jetPt"][:, 0], mc["jetPt"][:, 1]])
    mc_pt_w = np.concatenate([mc["kWeight"], mc["kWeight"]])
    overlay_plot(data_pt, data_pt_w, mc_pt, mc_pt_w,
                 np.arange(300, 2000, 20),
                 r"$p_T^{\,j}$ [GeV]", "Normalized density [1/GeV]",
                 os.path.join(FIGS_DIR, "jetPt_data_vs_mc.pdf"), log_y=True, xlim=(300, 2000))

    data_eta = np.concatenate([data["jetEta"][:, 0], data["jetEta"][:, 1]])
    mc_eta = np.concatenate([mc["jetEta"][:, 0], mc["jetEta"][:, 1]])
    overlay_plot(data_eta, data_pt_w, mc_eta, mc_pt_w,
                 np.arange(-2, 2.05, 0.05),
                 r"$\eta^{\,j}$", "Normalized density", xlim=(-2, 2),
                 outpath=os.path.join(FIGS_DIR, "jetEta_data_vs_mc.pdf"))

    overlay_plot(data["theta21"], data["kWeight"], mc["theta21"], mc["kWeight"],
                 np.linspace(-np.pi, np.pi, 41),
                 r"$\theta_{21}$ [rad]", "Normalized density [1/rad]",
                 os.path.join(FIGS_DIR, "theta21_data_vs_mc.pdf"), xlim=(-np.pi, np.pi),
                 ylim=(0.14, 0.18))

    overlay_plot(data["spva_lead"], data["kWeight"], mc["spva_lead"], mc["kWeight"],
                 np.linspace(-np.pi, np.pi, 41),
                 r"$\theta_s^{\,\mathrm{lead}}$ [rad]", "Normalized density [1/rad]",
                 os.path.join(FIGS_DIR, "jetSPVA_lead_data_vs_mc.pdf"), xlim=(-np.pi, np.pi))

    print(f"Wrote 5 comparison figures to {FIGS_DIR}")


if __name__ == "__main__":
    make_comparison_plots()
