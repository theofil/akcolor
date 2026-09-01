# MOD — session synopsis

## What we did

**Goal:** replicate the `makefriends.cpp` friends-building logic against public
Zenodo CMS Open Data instead of the group's own generator ntuples, since the two
input formats are structurally incompatible (see `README.md` "Showstoppers").

1. **Source data** (all downloaded to `MOD/raw/`, ~85 GB total):
   - Real data: [Zenodo 3340205](https://zenodo.org/records/3340205) — CMS 2011A,
     Jet300 HLT, pT>375 GeV, 18 files.
   - 8 companion Pythia6 QCD-MC records, binned in hard-parton p̂T (170-300 through
     1800-inf GeV). Each record has both `GEN*` (truth-only, unused) and `SIM*`
     (detector-level, used) files — `download_zenodo.py` fetches only `SIM*`.

2. **Converter** (`make_mod_friends.py`): reads `jets_i/jets_f/pfcs/pfcs_index`
   directly (no FastJet reclustering — MOD only stores the two hardest jets'
   own constituents), reproduces `fillPV`/pull-vector/`jcs*` arrays and
   `mjj/dYjj/dPhijj/ptjj` exactly as `makefriends.cpp` does, writes
   `MOD/friends/<sample>.h5` in the same schema as `root_to_hdf5.py`'s output
   (minus the jet-slot-2 axis: only 2 jets ever exist here, per instruction).
   `bosonM/Y/Eta/Pt/Phi` are always `-99` — no boson truth exists in this dataset
   family at all.

   **Bug hit and fixed:** first full run OOM-killed partway through the 470-600 GeV
   bin (this container has a ~35GB memory cgroup limit, separate from the 57GB
   host RAM `free` reports) because the converter buffered a whole sample's
   events in Python lists before writing. Rewrote to stream each input file
   straight into resizable HDF5 datasets (`append_h5`); verified byte-identical
   output against the pre-fix version on a known-good sample before rerunning.

3. **Final friend files, all 9 samples, no NaNs, no further crashes:**

   | Sample | Events |
   |---|---:|
   | `CMS2011A_data` | 455,359 |
   | `QCD_170-300` | 124 |
   | `QCD_300-470` | 544,619 |
   | `QCD_470-600` | 3,290,625 |
   | `QCD_600-800` | 3,791,301 |
   | `QCD_800-1000` | 3,916,666 |
   | `QCD_1000-1400` | 1,952,884 |
   | `QCD_1400-1800` | 1,989,694 |
   | `QCD_1800-inf` | 995,393 |

4. **Plots** (`MOD/figs/`, all DATA-vs-combined-QCD-MC overlays, styled like
   `summer26/plots.py`: MC = step histogram, DATA = black circle markers +
   Sumw2 stat error bars; no standalone data-only plots are produced):
   - `mjj_data_vs_mc.pdf`, `jetPt_data_vs_mc.pdf`, `jetEta_data_vs_mc.pdf` —
     data/MC agree closely across ~4-5 orders of magnitude.
   - `theta21_data_vs_mc.pdf` (y-axis zoomed to 0.14-0.18) — both flat at
     ~0.16/rad (≈1/2π, i.e. no strong pull/color correlation), but the zoom
     reveals a real, not-hugely-significant shape difference: MC bumps up
     around θ21≈0 where data dips, comparable in size to data's error bars.
   - `jetSPVA_lead_data_vs_mc.pdf` — mild U-shape, both data and MC agree,
     consistent with the weak discriminating power (AUC 0.547) this variable
     showed in `studies/O_jet.md`.
   - `jetShape_data_vs_mc.pdf` / `jetShapeRaw_data_vs_mc.pdf` — 2D analogs of
     `makefriends.cpp`'s `TH2Poly` jet-shape maps, built directly from the
     already-stored `jcsDEta/jcsDPhi/jcsW` (no ROOT file exists here to hold a
     `TH2Poly`, and root_to_hdf5.py never carries those into `.h5` either way).
     **2026-08-01: restyled to match `summer26/plots.py`'s raw `TH2Poly->Draw('COLZ')`**
     — same annular-sector geometry (128 phi wedges × 60 rho rings, r<0.5,
     rendered as true circular wedges via a curvilinear `pcolormesh` grid
     instead of the old rectangular Δη-Δφ grid), linear z-scale, `jet`
     rainbow colormap. Pull weight (`jcsW`) already matched `makefriends.cpp`
     exactly (`w = (pt_i/pt_jet)·r`, verified against the C++ source) — no
     change needed there. Remaining disclosed approximation: pseudorapidity Δη
     used in place of rapidity Δy (negligible difference, same sign convention).
     Data/MC agree well; sharp circular jet-radius profile as expected.
     **Real feature, not a bug**: under the linear scale a faint plus/cross
     enhancement is visible along Δη=0 and Δφ=0 in both DATA and MC — confirmed
     (via an independent plain-rectangular sanity histogram straight from
     `jcsDEta/jcsDPhi/jcsW`) to be present in the raw data itself, not an
     artifact of the polar-binning code. Most likely a genuine CMS
     detector-reconstruction effect (calorimeter/tracker projective geometry):
     MOD's jets are real 2011 data + Pythia6 GEANT-simulated SIM, unlike
     summer26's pure generator-truth jets (no detector simulation at all),
     which is why summer26's jetShape looks smoothly radially symmetric and
     MOD's doesn't. Data/MC agreeing on it is a good cross-check, not a concern.

   `theta_21` definition (not previously used anywhere in this codebase, so this
   is the convention adopted): pull angle of the **subleading** jet relative to
   the direction toward the **leading** jet — `jetPVA[sub_slot] -
   atan2(dPhijj, dYjj)`, wrapped to (-π,π]. Flip if a different convention was
   intended.

5. **Found but not yet used:** the SIM `SIM*.h5` files also carry particle-level
   truth — `gens`/`gens_index` (truth particles, same columns as `pfcs`),
   `gen_jet_pt/y/phi/m/eta` (truth-matched jet kinematics), and `hard_pid`/
   `hard_pt/y/phi` (hard-parton flavor/kinematics — a quark/gluon truth label).
   None of this is in the friend files yet; deliberately out of scope until now.

## What's next

**Unfold `jetSPVA` (leading-jet signed pull angle) and `theta_21` to particle level.**
Method chosen: **Iterative Bayesian Unfolding (IBU / D'Agostini)**, not OmniFold
(ruled out for now — much heavier build, would need the same GPU/condor training
infra as the `summer26/NN*` classifiers).

Plan (nothing implemented yet):
1. **Truth extraction (prerequisite, not started):** compute truth-level
   `jetPVA_gen/jetSPVA_gen/jetPVM_gen` and `theta21_gen` per SIM sample from
   `gens`/`gen_jet_*`, using the same `fillPV`-equivalent formula already in
   `make_mod_friends.py` — fed truth particles instead of `pfcs`. Must be paired
   jet-by-jet with the existing reco values (not just an independent truth
   histogram), since IBU needs the (truth bin, reco bin) pair per simulated jet.
   `theta21_gen` additionally needs truth-level `dYjj_gen`/`dPhijj_gen`.
2. **Response matrix:** `R[i,j] = P(reco bin i | truth bin j)` from the paired
   SIM truth/reco values, weighted by `kWeight`, same binning as the current
   plots (`np.linspace(-π, π, 41)`).
3. **D'Agostini iteration:** unfold the DATA reco histogram through `R`, few
   iterations (stop by χ² convergence or a small fixed count to avoid noise
   amplification — the usual IBU bias/variance tradeoff).
4. **Uncertainty:** Poisson toy/bootstrap resampling of the DATA histogram,
   rerun IBU per toy, take the spread (not analytic propagation).
5. **Validate before trusting it:** closure test (unfold SIM reco with its own
   response matrix, check it recovers SIM truth), and check response-matrix
   diagonal dominance for both observables — `theta_21`'s smearing could be
   broad given how weak the raw discriminating power already looked
   (`studies/O_jet.md`: `|θ_s|` AUC only 0.547), which may force coarser bins.
