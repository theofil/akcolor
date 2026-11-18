# summer26 analysis — VBF H→inv vs QCD dijet

Analysis code for the summer26 campaign (Herwig: campaign-00023/00030; MadGraph5+Pythia: pythia-dipole-recoil, since 2026-09-27).
Prerequisite: complete Steps 1–3 from [README.md](README.md) to build `makefriends`
and produce `friends/summer26/*.{friend.root,h5}`.

---

## Campaign details

**Sources**:
- campaign-00023: `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/`
- campaign-00030: `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00030/merged/` (QCDHjj Herwig only)
- pythia-dipole-recoil: `/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged/<PROC>/<PROC>_mg5_pythia_dipolerecoil.root` — **all six MG5+Pythia samples** since 2026-09-27 (see changelog below); replaces the buggy campaign-00023/00031 Pythia files

| Process | σ Herwig (pb) | σ MG5+Pythia (pb) | Herwig file | MG5+Pythia file |
|---------|--------------|-------------------|-------------|----------------|
| VBF H→inv | 2.97189 ± 0.00325 | 2.88990 ± 0.00083 | `VBFH_herwig.root` (camp-00023) | `VBFH_mg5_pythia.root` (dipole-recoil) |
| VBF W→inv | 7.494 ± 0.00632 | 7.20333 ± 0.00198 | `VBFW_herwig.root` (camp-00023) | `VBFW_mg5_pythia.root` (dipole-recoil) |
| VBF Z→inv | 1.099 ± 0.00126 | 1.08581 ± 0.00030 | `VBFZ_herwig.root` (camp-00023) | `VBFZ_mg5_pythia.root` (dipole-recoil) |
| QCD H+jj | 5.06051 ± 0.06030 | 4.982 ± 0.0008692 | `QCDHjj_herwig.root` (camp-00030) | `QCDHjj_mg5_pythia.root` (dipole-recoil) |
| QCD W+jj | 1733.3 ± 1.84 | 1646.805 ± 0.510 | `QCDWjj_herwig.root` (camp-00023) | `QCDWjj_mg5_pythia.root` (dipole-recoil) |
| QCD Z+jj | 340.2 ± 0.316 | 316.363 ± 0.098 | `QCDZjj_herwig.root` (camp-00023) | `QCDZjj_mg5_pythia.root` (dipole-recoil) |

All 12 files available. The dipole-recoil Pythia files contain the same events with the same `evweight` normalization as the campaign-00023/00031 files they replace (`validation.json`: `normalization_matches_original: true`), so all cross sections are unchanged. `QCDHjj_mg5_pythia` has 409 621 events (xs = 4.982 ± 0.0008692 pb).

---

## jet2 + boson 4-vector branches (added 2026-07-06)

The friend trees and h5 files were extended with the 3rd pt-ordered jet ("jet2") and the
generator-level boson 4-vector. **Friends/h5 produced before this change lack these branches
and must be regenerated** (`./run.summer26.sh` + h5 conversion).

### jet2 branches

`NJETMAX` is now 3. New branches, same conventions as jet0/jet1:
`jetPt2, jetEta2, jetPhi2, jetM2, jetSPVA2, jetPVA2, jetPVM2, jetNC2` and the constituent
arrays `jcsPt2, jcsDEta2, jcsDEtaRaw2, jcsDPhi2, jcsM2, jcsW2` (each `[80]`).

- jet2 is the **pt-ordered 3rd jet** passing the same jet cuts (pt > 30 GeV, |η| < 3);
  it is **never randomized** — the leading/subleading swap (`leadJetIndex`) stays confined
  to jets 0/1 — so jet2 is always softer than both leading jets.
- When only 2 jets pass the cuts, all jet2 slots are zero (validity mask: `jcsPt2 > 0`,
  `jetPt2 > 0`), and `nJets` = 2; otherwise `nJets` = 3.
- The dijet event selection is unchanged (two leading jets only).
- h5 per-jet/per-constituent datasets now have shapes **(N, 3)** and **(N, 3, 80)**
  (was (N, 2) / (N, 2, 80)). Backward compatible: all `summer26/NN*/dataset.py` index
  `[:, 0]` / `[:, 1]` only.

### Boson 4-vector branches

New scalar event branches `bosonM, bosonY, bosonEta, bosonPt, bosonPhi` (M/pT 2dp, Y/Eta/Phi 3dp;
−99 sentinel when absent). The boson type is passed explicitly per sample via
`makefriends --boson H|W|Z` (no auto-detect); without `--boson` the branches are −99.

The `Data` tree stores **final-state particles only** (no status codes, no intermediate
resonances), so the boson is reconstructed as:

| Channel | What's in the record | Recipe |
|---------|---------------------|--------|
| **H** (VBFH, QCDHjj) | Higgs kept stable (models H→inv): exactly one pid 25/event, m = 125.00 exactly | Take the pid 25 particle directly |
| **W** (VBFW, QCDWjj) | No pid ±24; W→τν_τ with the τ kept stable | p4(τ) + p4(ν_τ), flavor-matched: pid 15 ↔ −16, pid −15 ↔ 16 |
| **Z** (VBFZ, QCDZjj) | No pid 23; Z→ν_e ν̄_e, always electron flavor (pid ±12) | p4(ν_e) + p4(ν̄_e), i.e. a (12, −12) pair |

Ambiguity resolution (extra ν's from hadron decays in ~13% of Z events; a rare 2nd shower τ
in W events): among flavor-consistent candidates, pick the pair with invariant mass **closest
to the pole mass** (m_W = 80.379, m_Z = 91.1876 GeV). Expected peaks: H = 125 exactly;
W median ≈ 80.4 GeV; Z median ≈ 91.2 GeV. Do **not** use the `theETmiss[4]` branch — its
components do not match the ν system.

---

## mjj>400 GeV & |ΔY_jj|>2.5 preselection (added 2026-09-01)

The dijet selection in `makefriends.cpp` was tightened to a standard VBF-style
preselection: **`mjj > 400 GeV`** and **`|ΔY_jj| > 2.5`**, applied in addition
to (not instead of) the existing opposite-sign-η requirement — see "Dijet
selection applied in `makefriends.cpp`" below for the full 4-cut list.
Previously the selection was `mjj > 0` (no real mass cut) with no ΔY cut.

**Everything downstream was regenerated from scratch under the new
selection**: all 12 friend ROOT/h5 files, all 18 trained networks, all
inference ROC curves, all per-channel event-level score columns, all
diagnostic plots, and the SR significance table in "Signal/background
estimation from ROC curves" below. Any such artifact with a timestamp before
2026-09-01 used the old, looser selection and is stale.

**Effect on efficiency** (from `figs/summer26/checksum.txt`): the cut is far
more aggressive on background than signal, as expected for a VBF topology
cut — QCD weighted efficiency drops to **~2.3–5.3%** (was ~15–26%), VBF
signal efficiency to **~14.7–22.1%** (was ~26–42%).

**Effect on raw signal statistics**: the preselection's ~30–40% loss of raw
signal (removed as an unavoidable side effect of removing far more
background — see `raw S₀` in "Results" below for the post-preselection
counts) is a real cost, not a bug. This is a real effect of trading some raw
signal statistics for a cleaner, more physically realistic VBF-topology
sample — not a regression.

---

## pythia-dipole-recoil fix: all MG5+Pythia samples replaced (2026-09-27)

A bug was found in all six MG5+Pythia samples (campaign-00023 and, for
QCDHjj, campaign-00031). They are replaced by the **pythia-dipole-recoil**
re-shower of the same events (`/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged/`, also reachable as
the `friends/summer26/pythia-dipole-recoil` symlink). Every entry count and
Σ`evweight` matches the old files, so all `--xs`/`--boson` values are unchanged.

**What was redone:**
- The 6 `*_mg5_pythia.{friend.root,h5}` files, rebuilt with
  `./run.summer26.sh --pythiaOnly`. The old buggy files are kept in
  `friends/summer26/pythia_prefix_bug/`.
- All Pythia inference, **with no retraining** (all nets are Herwig-trained):
  every `infer_gpu.sub` (incl. `Kj_W`, now added to `infer.txt`), plus the
  `NNj_jet*`/`NNjB`/`NNjjBj` score columns in the Pythia h5 files.
- Every plot, `checksum.txt`, `pythia_transfer.txt` and
  `Z_performance_resilience.txt`.

The Herwig friends/h5 files and their score columns were not touched.

**Effect:**
- Fewer Pythia events pass the preselection: 5–17% fewer than with the old
  files (e.g. VBFH 171 085 → 142 565, QCDZjj 32 309 → 30 452; see
  `checksum.txt`).
- The Herwig-trained nets transfer noticeably better when they use the 3rd
  jet. Pythia AUC for `NNjjBj`: Z 0.837→0.867, W 0.827→0.864, H 0.869→0.896
  (Herwig 0.870/0.875/0.927); `KjjBj` Z 0.828→0.864.
- Nets that see only the two leading jets move down slightly (e.g. `NNj_H`
  0.773→0.757, `NNkin_H` 0.804→0.786).
- The Pythia-transfer and resilience tables below were regenerated.

**Herwig results are unchanged:** every Herwig AUC is identical. The one
visible Herwig change is a GPU-nondeterminism artifact. Re-running
inference reproduces each Herwig AUC to ~1e-8, but the decimated
(≤5000-point) ROC grid in `roc_data_*.npz` moves by a few points. For
`KjB`-Z only, the S/√B scan therefore lands on a neighbouring grid point:
58.50 → 58.48, raw S/B 20789/911 → 20771/910. All other rows of
`SB_estimate.txt` are identical.

To reproduce the Pythia side only:
`./run.summer26.sh --pythiaOnly`, which stages into `friends/summer26/pythia_stage/`.
Move the staged files into `friends/summer26/` by hand. The staging step is
needed because `root_to_hdf5.py` rewrites every h5 in its folder and would
drop the Herwig score columns. Then rerun the `infer.txt` jobs.

---

## Step 2 — Run event reconstruction

### Quick test (100 events per sample, foreground)

```bash
./run.summer26.sh --goFast
# or a custom number of events:
./run.summer26.sh --goFast 500
```

### Full production (20 parallel jobs × 50k events, then hadd)

```bash
./run.summer26.sh
# MG5+Pythia samples only, staged in friends/summer26/pythia_stage/:
./run.summer26.sh --pythiaOnly
```

**Output** (all 12 files):

| File | Description |
|------|-------------|
| `friends/summer26/VBFH_herwig.friend.root` | VBF H Herwig friend tree |
| `friends/summer26/VBFH_herwig.h5` | VBF H Herwig HDF5 |
| `friends/summer26/VBFH_mg5_pythia.friend.root` | VBF H MG5+Pythia friend tree |
| `friends/summer26/VBFH_mg5_pythia.h5` | VBF H MG5+Pythia HDF5 |
| `friends/summer26/VBFW_herwig.friend.root` | VBF W Herwig friend tree |
| `friends/summer26/VBFW_herwig.h5` | VBF W Herwig HDF5 |
| `friends/summer26/VBFW_mg5_pythia.friend.root` | VBF W MG5+Pythia friend tree |
| `friends/summer26/VBFW_mg5_pythia.h5` | VBF W MG5+Pythia HDF5 |
| `friends/summer26/VBFZ_herwig.friend.root` | VBF Z Herwig friend tree |
| `friends/summer26/VBFZ_herwig.h5` | VBF Z Herwig HDF5 |
| `friends/summer26/VBFZ_mg5_pythia.friend.root` | VBF Z MG5+Pythia friend tree |
| `friends/summer26/VBFZ_mg5_pythia.h5` | VBF Z MG5+Pythia HDF5 |
| `friends/summer26/QCDHjj_herwig.friend.root` | QCD H+jj Herwig friend tree |
| `friends/summer26/QCDHjj_herwig.h5` | QCD H+jj Herwig HDF5 |
| `friends/summer26/QCDHjj_mg5_pythia.friend.root` | QCD H+jj MG5+Pythia friend tree (dipole-recoil, xs=4.982 ± 0.0008692 pb) |
| `friends/summer26/QCDHjj_mg5_pythia.h5` | QCD H+jj MG5+Pythia HDF5 (dipole-recoil) |
| `friends/summer26/QCDWjj_herwig.friend.root` | QCD W+jj Herwig friend tree |
| `friends/summer26/QCDWjj_herwig.h5` | QCD W+jj Herwig HDF5 |
| `friends/summer26/QCDWjj_mg5_pythia.friend.root` | QCD W+jj MG5+Pythia friend tree |
| `friends/summer26/QCDWjj_mg5_pythia.h5` | QCD W+jj MG5+Pythia HDF5 |
| `friends/summer26/QCDZjj_herwig.friend.root` | QCD Z+jj Herwig friend tree |
| `friends/summer26/QCDZjj_herwig.h5` | QCD Z+jj Herwig HDF5 |
| `friends/summer26/QCDZjj_mg5_pythia.friend.root` | QCD Z+jj MG5+Pythia friend tree |
| `friends/summer26/QCDZjj_mg5_pythia.h5` | QCD Z+jj MG5+Pythia HDF5 |

---

## Step 5 — Train the neural networks

All NNs use the LCG_106_cuda environment. Training is submitted to a GPU node
via HTCondor; the full submission list (all networks × Z/H/W) is in `train.txt`.
**Plain `condor_submit` fails from `/eos` paths** ("Standard batch schedds
cannot use /eos paths directly") — submit via the EosSubmit `bigbird24`
schedd instead: `condor_submit -name bigbird24.cern.ch <sub file>`. To submit
the base (Z) trainings:

```bash
condor_submit -name bigbird24.cern.ch summer26/NNkin/train_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNj/train_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjB/train_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjj/train_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjB/train_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjBj/train_gpu.sub
```

Or run interactively (falls back to CPU):

```bash
cd summer26/NNkin   && python train.py
cd summer26/NNj     && python train.py
cd summer26/NNjB    && python train.py
cd summer26/NNjj    && python train.py
cd summer26/NNjjB   && python train.py
cd summer26/NNjjBj  && python train.py
```

The `_H` / `_W` variant directories (generated by `summer26/make_variants.py`)
train the same architectures on the H and W samples. Exception: `NNjB_H`,
`NNjjB_H` and `NNjjBj_H` are excluded from generation and maintained manually —
they keep `bosonM`, which the base dirs dropped; the script enforces this
invariant on every run and fails loudly if violated.

### NNkin (`summer26/NNkin/`)

Kinematic-only MLP. 12 input features: `dPhijj, dYjj, mjj, ptjj` plus
`eta, m, phi, pt` for each of the two leading jets.

Architecture: **12 → 128 → 64 → 32 → 1** (BatchNorm + ReLU + Dropout(0.3) at each hidden layer).

Trained on **QCDZjj_herwig + VBFZ_herwig**.

### NNj (`summer26/NNj/`)

DeepSets over the leading-jet constituents using only low-level features — no
pull-vector jet scalars (NC, |t⃗|, θ_s removed) and no constituent weight (w
removed).

**Jet scalars** (3, leading jet): |η|, m, pT  
**Constituent features** (80 max × 3): Δη, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP**: 3 → 64 → 64 (ReLU, shared across constituents)
- masked-sum pooling → (B, 64)
- **rho MLP**: 67 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3))

Trained on **QCDZjj_herwig + VBFZ_herwig**.

### NNjB (`summer26/NNjB/`)

Identical inputs to NNj plus the generator-boson pT, η, φ as additional
jet-level scalars. `bosonM` is excluded (on-shell generation artifact in the
Herwig QCD Z/W samples — see tuesday.md); only the `_H` variant keeps it
(harmless constant 125).

**Jet-level scalars** (7 = 4 jet + 3 boson): |η|, m, pT, φ, boson pT, boson η, boson φ  
**Constituent features** (80 max × 3): Δη, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP**: 3 → 64 → 64 (ReLU, shared across constituents)
- masked-sum pooling → (B, 64)
- **rho MLP**: 71 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3)); 72 in `_H` (bosonM kept)

Trained on **QCDZjj_herwig + VBFZ_herwig**.

### NNjj (`summer26/NNjj/`)

DeepSets over the constituents of **both** the leading (index 0) and sub-leading
(index 1) jets, using only raw low-level features — no pull-vector scalars, no
constituent weight, and signed η (not |η|). The phi MLP is shared across all
constituents from both jets.

**Jet scalars** (4 per jet × 2 jets = 8 total): η (signed), m, pT, φ  
**Constituent features** (80 max × 3 per jet): Δη_raw, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP** (shared): 3 → 64 → 64 (ReLU, applied independently per jet)
- masked-sum pooling per jet → pool₀ (B, 64), pool₁ (B, 64)
- concat [pool₀, pool₁, jet₀ scalars, jet₁ scalars] → (B, 136)
- **rho MLP**: 136 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3))

Trained on **QCDZjj_herwig + VBFZ_herwig**.

### NNjjB (`summer26/NNjjB/`)

Extension of NNjj that adds the generator-boson pT, η, φ as event-level
scalars. `bosonM` is excluded (on-shell generation artifact in the Herwig QCD
Z/W samples — see tuesday.md); only the `_H` variant keeps it (harmless
constant 125).

**Jet scalars** (4 per jet × 2 jets = 8 total): η (signed), m, pT, φ  
**Boson scalars** (3, event-level): boson pT, η, φ  
**Constituent features** (80 max × 3 per jet): Δη_raw, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP** (shared): 3 → 64 → 64 (ReLU, applied independently per jet)
- masked-sum pooling per jet → pool₀ (B, 64), pool₁ (B, 64)
- concat [pool₀, pool₁, jet₀ scalars, jet₁ scalars, boson scalars] → (B, 139)
- **rho MLP**: 139 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3)); 140 in `_H` (bosonM kept)

Trained on **QCDZjj_herwig + VBFZ_herwig**.

### NNjjBj (`summer26/NNjjBj/`)

Extension of NNjjB that adds the **3rd pt-ordered jet** — its 4-momentum
(pT, η, φ, m) as event-level scalars **and its constituents** through the shared
phi MLP — when a 3rd jet exists. Jets have pT > 30 GeV, so the all-zero rows
written when no 3rd jet passes the cuts (`jetPt[:,2] == 0`, ~75% of QCDZjj
events) unambiguously encode absence; the StandardScaler is fit on all training
rows so absence maps to a constant code, the jet-2 constituent mask is all-False
(zero pooled vector), and the jcsPt/jetPt division is guarded against 0/0.
`bosonM` is excluded as in NNjjB; only the `_H` variant keeps it.

**Jet scalars** (4 per jet × 2 jets = 8 total): η (signed), m, pT, φ
**3rd-jet scalars** (4, event-level): η (signed), m, φ, pT — zeros when absent
**Boson scalars** (3, event-level): boson pT, η, φ
**Constituent features** (80 max × 3 per jet, all 3 jets): Δη_raw, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP** (shared): 3 → 64 → 64 (ReLU, applied independently per jet, all 3 jets)
- masked-sum pooling per jet → pool₀, pool₁, pool₂ (B, 64) each; pool₂ = 0 when no 3rd jet
- concat [pool₀, pool₁, pool₂, jet₀ scalars, jet₁ scalars, jet₂ scalars, boson scalars] → (B, 207)
- **rho MLP**: 207 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3)); 208 in `_H` (bosonM kept)

Trained on **QCDZjj_herwig + VBFZ_herwig**.

**Outputs saved by training** (all NNs use `_{W/Z/H}` suffix matching `PROCESS`):

| File | Contents |
|------|----------|
| `best_model_{Z/H}.pt` | Model weights at best validation loss |
| `scaler_{Z/H}.pkl` | Fitted `StandardScaler` dict |
| `split_indices_{Z/H}.npz` | `train_idx`, `val_idx`, `test_idx` |
| `figs/summer26/{NN}_{Z/H}_loss_curve.pdf` | Train/val loss vs epoch |

---

## Step 6 — Run inference

Submit all inference jobs to HTCondor (commands are in `infer.txt`), via the
EosSubmit `bigbird24` schedd (plain `condor_submit` fails from `/eos` paths —
see Step 5):

```bash
condor_submit -name bigbird24.cern.ch summer26/NNkin/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNkin_H/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNkin_W/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNj/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNj/save_scores.sub
condor_submit -name bigbird24.cern.ch summer26/NNj_H/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNj_H/save_scores.sub
condor_submit -name bigbird24.cern.ch summer26/NNj_W/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNj_W/save_scores.sub
condor_submit -name bigbird24.cern.ch summer26/NNjB/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjB/save_scores.sub
condor_submit -name bigbird24.cern.ch summer26/NNjB_H/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjB_H/save_scores.sub
condor_submit -name bigbird24.cern.ch summer26/NNjB_W/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjB_W/save_scores.sub
condor_submit -name bigbird24.cern.ch summer26/NNjj/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjj_H/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjj_W/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjB/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjB_H/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjB_W/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjBj/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjBj/save_scores.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjBj_H/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjBj_H/save_scores.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjBj_W/infer_gpu.sub
condor_submit -name bigbird24.cern.ch summer26/NNjjBj_W/save_scores.sub
```

**Caution**: within a channel (Z/H/W), the `save_scores.sub` jobs of `NNj`,
`NNjB` and `NNjjBj` all write to the *same* 4 h5 files for that channel — do
not submit same-channel `save_scores` jobs concurrently (risks
concurrent-write corruption on EOS FUSE); submit in channel-disjoint batches
or one at a time instead. `infer.txt` carries the same note.

Each job requests 2 CPUs, 8 GB RAM, 1 GPU (`microcentury` flavour). Each `inference.py`
evaluates on the Herwig test split (training domain) and the full MG5+Pythia dataset
(generator transfer test), producing `figs/summer26/roc_<nn>_{W,Z,H}.pdf`. It also
saves the model ROC curves (decimated to ≤5000 points) to `roc_data_{W,Z,H}.npz` in
its own directory; `summer26/plots.py` reads these to draw the all-NN summary
comparison `figs/summer26/roc_summary_{W,Z,H}.pdf` (model curves only, no
single-observable baselines). Each net's Herwig curve is marked with a ★ at its
`S/√B`-optimal operating point (see "Signal/background estimation from ROC
curves" below).

---

## kWeight factors and Run 3 statistical uncertainty

### Where `evweight` lives

The `evweight` branch is stored in the **original campaign ROOT files on EOS** — not in the friend trees. The TTree in these files is named **`Data`**. The exact file paths per sample are:

| Sample | EOS file |
|--------|----------|
| `VBFH_herwig` | `campaign-00023/merged/VBFH_herwig.root` |
| `VBFH_mg5_pythia` | `pythia-dipole-recoil/merged/VBFH/VBFH_mg5_pythia_dipolerecoil.root` |
| `VBFW_herwig` | `campaign-00023/merged/VBFW_herwig.root` |
| `VBFW_mg5_pythia` | `pythia-dipole-recoil/merged/VBFW/VBFW_mg5_pythia_dipolerecoil.root` |
| `VBFZ_herwig` | `campaign-00023/merged/VBFZ_herwig.root` |
| `VBFZ_mg5_pythia` | `pythia-dipole-recoil/merged/VBFZ/VBFZ_mg5_pythia_dipolerecoil.root` |
| `QCDWjj_herwig` | `campaign-00023/merged/QCDWjj_herwig.root` |
| `QCDWjj_mg5_pythia` | `pythia-dipole-recoil/merged/QCDWjj/QCDWjj_mg5_pythia_dipolerecoil.root` |
| `QCDZjj_herwig` | `campaign-00023/merged/QCDZjj_herwig.root` |
| `QCDZjj_mg5_pythia` | `pythia-dipole-recoil/merged/QCDZjj/QCDZjj_mg5_pythia_dipolerecoil.root` |
| `QCDHjj_herwig` | `campaign-00030/merged/QCDHjj_herwig.root` |
| `QCDHjj_mg5_pythia` | `pythia-dipole-recoil/merged/QCDHjj/QCDHjj_mg5_pythia_dipolerecoil.root` |

where `campaign-00023`, `campaign-00030` are under
`/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/`
and `pythia-dipole-recoil` is `/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/`
(replaces the buggy campaign-00023/00031 Pythia files since 2026-09-27).

### kWeight formula

`makefriends` computes `kWeight` for each event (line 308 of `makefriends.cpp`):

```
kWeight_i = evweight_i × σ / sumWtot
```

where:

- **`evweight_i`** — per-event generator weight read from branch `evweight` in TTree `Data`
  (passed as `--genWeight evweight` for all 12 samples in `run.summer26.sh`)
- **`σ`** — process cross-section in pb, passed via `--xs` (exact values in the table below)
- **`sumWtot`** — sum of `evweight` over **all** events in the input file **before any selection**,
  computed at startup in `makefriends.cpp` line 129 via
  `ROOT::RDataFrame::Sum<double>("evweight").GetValue()`

### Per-sample values

Exact `σ` values (from `run.summer26.sh`) and measured `evweight` statistics (from reading all
original EOS files; histograms in `figs/summer26/evweight_{W,Z,H}.pdf`).
`kWeight` statistics are from the **post-selection** `events` tree in the friend ROOT files
(`friends/summer26/*.friend.root`); histograms in `figs/summer26/kWeight_{W,Z,H}.pdf`.
Run `python summer26/weights.py` to regenerate both sets of PDFs and reprint this table.

Post-selection (`N sel.`, `kWeight mean`, `kWeight std`) columns reflect the
current `mjj>400 GeV` + `|ΔY_jj|>2.5` preselection (added 2026-09-01); `evweight`
columns are selection-independent (computed over all events before any cut) and
unchanged from before. MG5+Pythia rows updated 2026-09-27 for the
dipole-recoil samples (same `evweight`, fewer events pass the preselection).

| Sample | N events | σ `--xs` (pb) | evweight min | evweight max | evweight mean | evweight std | N sel. | kWeight mean | kWeight std |
|--------|----------|--------------|-------------|-------------|--------------|-------------|--------|-------------|------------|
| `VBFH_herwig` | 1 000 000 | 2.97189 | 1.000 | 424.94 | 1.0056 | 0.622 | 146 901 | 2.9668e-06 | 1.3108e-06 |
| `VBFH_mg5_pythia` | 1 000 000 | 2.88990 | 2.8870 | 2.8920 | 2.8900 | 0.0015 | 142 565 | 2.8899e-06 | 1.466e-09 |
| `VBFW_herwig` | 1 000 000 | 7.494 | 1.000 | 91.604 | 1.0060 | 0.192 | 176 406 | 7.4755e-06 | 9.3813e-07 |
| `VBFW_mg5_pythia` | 1 000 000 | 7.20333 | 7.1962 | 7.2120 | 7.2033 | 0.0051 | 176 575 | 7.2033e-06 | 5.0422e-09 |
| `VBFZ_herwig` | 1 000 000 | 1.099 | 1.000 | 46.794 | 1.0183 | 0.194 | 205 250 | 1.1017e-06 | 2.182e-07 |
| `VBFZ_mg5_pythia` | 1 000 000 | 1.08581 | 1.0846 | 1.0868 | 1.0856 | 0.00058 | 207 827 | 1.0858e-06 | 5.7941e-10 |
| `QCDHjj_herwig` | 1 000 000 | 5.06051 | 1.000 | 12.842 | 1.0002 | 0.026 | 23 340 | 5.0596e-06 | ~0 |
| `QCDHjj_mg5_pythia` | 409 621 | 4.982 | 4.9818 | 4.9818 | 4.9818 | ~0 | 19 945 | 1.2162e-05 | 0 |
| `QCDWjj_herwig` | 1 000 000 | 1733.3 | 0.00146 | 6.3698 | 0.11363 | 0.011 | 22 603 | 0.0017315 | 4.4768e-05 |
| `QCDWjj_mg5_pythia` | 1 000 000 | 1646.805 | 1644.5 | 1649.8 | 1646.8 | 1.73 | 28 712 | 0.0016468 | 1.7307e-06 |
| `QCDZjj_herwig` | 1 000 000 | 340.2 | 0.06666 | 19.326 | 0.06679 | 0.025 | 22 821 | 0.00033969 | 6.7361e-06 |
| `QCDZjj_mg5_pythia` | 1 000 000 | 316.363 | 316.00 | 316.72 | 316.37 | 0.232 | 30 452 | 0.00031636 | 2.3191e-07 |

### Character of the weights by generator

**MG5+Pythia**: `evweight` is essentially constant across events and equals the process
cross-section σ (std ≈ 0). This means `sumWtot ≈ N × σ` and therefore
`kWeight_i ≈ σ/N` — all events contribute equally to the absolute normalization.

**Herwig VBF** (`VBFH`, `VBFW`, `VBFZ`, `QCDHjj`): `evweight` has a minimum of 1 but a wide
upward tail (up to 425 for `VBFH_herwig`). These are per-event parton-shower importance
weights. The mean is close to 1 so `sumWtot ≈ N`, and `kWeight_i = evweight_i × σ / N`.

**Herwig QCD** (`QCDWjj`, `QCDZjj`): `evweight` mean is well below 1 (0.114 for
`QCDWjj_herwig`, 0.067 for `QCDZjj_herwig`), reflecting Herwig's inclusive QCD sampling
strategy where most generated events carry low weight.

### Dijet selection applied in `makefriends.cpp`

Only events passing all four cuts enter the output `events` TTree (and thus contribute to histograms):

1. ≥ 2 anti-kT (R=0.4) jets with pT > 30 GeV and |η| < 3.0
2. Two leading jets with opposite-sign pseudorapidity: η₀ × η₁ < 0
3. Dijet invariant mass: mjj > 400 GeV
4. Dijet rapidity separation: |ΔY_jj| > 2.5

Cuts 3–4 were added 2026-09-01 (previously `mjj > 0`, no ΔY cut — see the
changelog section above); cut 2 (opposite-sign η) is unchanged and applies
in addition to, not instead of, the new mjj/ΔY cuts.
Summing `kWeight` over all selected events gives `σ_eff` (effective cross-section after selection, in pb).

---

## Signal/background estimation from ROC curves

Script: `summer26/sb_estimate.py` (numpy + h5py only, no ROOT/torch/matplotlib needed).
**Herwig only** — every net is trained on Herwig, and this table does not evaluate
MadGraph5+Pythia at all (no generator-transfer numbers are computed here). The only
figure of merit used is the background-dominated significance `S/√B`. For every (net, channel)
combination — the same 18 `roc_data_{H,W,Z}.npz` files used by the
`roc_summary_{H,W,Z}.pdf` figures — the script scans the net's own Herwig ROC
curve (signal efficiency `tpr` vs background efficiency `fpr`, no additional
mjj/ΔYjj cuts) and reports two operating points per row, at **L = 300 fb⁻¹**: the
*initial*, unoptimized point at the trivial no-NN-cut baseline (`_0`), and the
*final* point — the ROC operating point that **maximizes `S/√B`**, subject to a
`raw B ≥ 10` statistical-validity floor (same guard convention as the removed
`SR_optimization.txt` study). `S` and `B` come from each sample's post-selection
`Σ kWeight` (pb, see kWeight section above) × `L` × the ROC efficiency at that point.
`analyze_one(net, channel)` is imported directly by `plots.py`, which marks the
*final* point with a ★ on the corresponding Herwig curve in `roc_summary_{H,W,Z}.pdf`
(see "Step 6" above). Raw output: `figs/summer26/SB_estimate.txt`.

### Branching-ratio rescale (decay-chain-aware)

The samples were generated with the boson decayed in an "invisible-like" mode: H is kept
fully **stable** (BSM-invisible-Higgs convention — its tabulated cross section is
therefore the true **inclusive** VBF/QCD production rate, unaffected by any decay BR); Z
is reconstructed specifically as **Z→ν_eν̄_e** and W as **W→τν_τ** (real SM decay chains —
see "Boson 4-vector branches" above — so their tabulated cross sections already have that
one channel's real branching fraction baked in by the matrix element). To reinterpret
these yields for a realistic H→γγ / Z→ee,μμ / W→eν,μν search, each channel's `S` and `B`
are rescaled by:

| Channel | Tabulated σ represents | Rescale formula | BR values used (PDG) | Factor |
|---------|------------------------|------------------|-----------------------|-------:|
| **H** | inclusive VBF/QCD H production | `BR(H→γγ)` | BR(H→γγ) = 2.27×10⁻³ (mH=125.09 GeV) | **0.00227** |
| **Z** | exclusive Z→ν_eν̄_e production | `[BR(Z→ee)+BR(Z→μμ)] / BR(Z→ν_eν̄_e)` | BR(Z→ee)=3.363%, BR(Z→μμ)=3.366%, BR(Z→νν̄)ₜₒₜₐₗ=20.00% ÷ 3 flavors = 6.667% | **1.0093** |
| **W** | exclusive W→τν_τ production | `[BR(W→eν)+BR(W→μν)] / BR(W→τν)` | BR(W→eν)=10.71%, BR(W→μν)=10.63%, BR(W→τν)=11.38% | **1.8752** |

The **same** per-channel factor is applied to both the VBF signal and the QCD background,
since QCD H/W/Z+jj is only a real background to a γγ/ee,μμ/eν,μν search once its boson
also decays to that same final state.

**Rescale-invariance caveat:** because the identical factor multiplies `S` and `B` within
a channel, it cancels exactly in the `S/B` ratio and does not shift which ROC point
maximizes `S/√B` (`significance_final = √factor × significance_unrescaled` at every
point, so the argmax is unchanged — verified numerically by `sb_estimate.py`'s
`verify_rescale_invariance()`). The scan therefore runs once on the un-rescaled yields;
the rescale is applied only to the winning point's final `S`, `B` before display. In other
words, this whole exercise changes the absolute event counts to a physically realistic
decay channel, but says nothing about which NN or working point is best — that ordering
is exactly what it would be without any BR at all.

**Statistical floor caveat:** `S/√B` has no `+S` term discouraging a near-zero-background
cut, so the scan does push there — the `raw B ≥ 10` guard binds for 6 of the 18 rows below
(`NNkin`-H, `NNjjB`-H, `NNjjBj`-H, `NNjjB`-W, `NNjjB`-Z, `NNjjBj`-Z), each marked `†` on
`raw B` in the table below.

**Decimated-ROC-grid caveat:** `roc_data_{H,W,Z}.npz` stores only the decimated
(≤5000-point) `fpr`/`tpr` arrays saved by each net's `inference.py` — the underlying
sklearn `roc_curve` **threshold values are discarded and never saved**. The table below
therefore reports the operating point as a **(signal efficiency, background efficiency)**
pair, not an NN-score cut value, and the true continuous optimum may differ slightly from
the grid optimum.

**No detector effects, no other SM backgrounds:** `S` and `B` here are MC-truth,
particle-level yields from exactly two samples per channel (VBF-produced vs
QCD-produced boson+dijet) — there is no detector simulation and no competing SM
background process (e.g. real diphoton QCD for H→γγ, Drell-Yan for Z→ll, W+jets/ttbar
for W→lν). The resulting significance values (up to a few hundred for the W/Z channels,
driven by their much larger QCD cross sections) are **not** discovery-significance
estimates for a real search — they measure the NNs' relative discrimination power only.

### Results (6 nets × 3 channels, Herwig only, S/√B initial vs. final at L = 300 fb⁻¹)

Regenerated 2026-09-02. Two points are
reported per row: `_0` (initial, the unoptimized no-NN-cut baseline, identical
across the 6 nets within a channel) and the unsuffixed *final* columns (the
`S/√B`-maximizing ROC operating point). The underlying `mjj>400 GeV` +
`|ΔY_jj|>2.5` preselection and yields are unchanged from the 2026-09-01 version —
only the table's shape changed. The *final* point of each
row is also marked with a ★ on the corresponding Herwig curve in
`roc_summary_{H,W,Z}.pdf` (see "Step 6" above).

Subscript legend: `_0` = initial, unoptimized baseline (full preselected sample, no
NN cut, same for every net within a channel); unsuffixed columns = final, the ROC
operating point maximizing `S/√B` — the only optimization objective. `†` on
`raw B` marks a row where the `raw B ≥ 10` floor changed the unconstrained
optimum.

| Net | Ch | S/√B₀ | S₀ | B₀ | raw S₀ | raw B₀ | S/√B | S | B | raw S | raw B |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `NNkin` | **H** | 33.10 | 296.8 | 80.4 | 22036 | 3501 | **97.90** | 53.5 | 0.3 | 3972 | 13† |
| `NNj` | **H** | 33.10 | 296.8 | 80.4 | 22036 | 3501 | **56.53** | 133.8 | 5.6 | 9937 | 244 |
| `NNjB` | **H** | 33.10 | 296.8 | 80.4 | 22036 | 3501 | **61.84** | 40.9 | 0.4 | 3033 | 19 |
| `NNjj` | **H** | 33.10 | 296.8 | 80.4 | 22036 | 3501 | **119.73** | 72.6 | 0.4 | 5389 | 16 |
| `NNjjB` | **H** | 33.10 | 296.8 | 80.4 | 22036 | 3501 | **131.48** | 82.2 | 0.4 | 6100 | 17† |
| `NNjjBj` | **H** | 33.10 | 296.8 | 80.4 | 22036 | 3501 | **155.10** | 74.3 | 0.2 | 5519 | 10† |
| `NNkin` | **W** | 158.10 | 741 869.6 | 22 017 470.1 | 26461 | 3391 | **180.77** | 526 410.3 | 8 479 745.2 | 18776 | 1306 |
| `NNj` | **W** | 158.10 | 741 869.6 | 22 017 470.1 | 26461 | 3391 | **169.62** | 626 023.4 | 13 622 132.8 | 22329 | 2098 |
| `NNjB` | **W** | 158.10 | 741 869.6 | 22 017 470.1 | 26461 | 3391 | **218.54** | 455 478.4 | 4 343 759.2 | 16246 | 669 |
| `NNjj` | **W** | 158.10 | 741 869.6 | 22 017 470.1 | 26461 | 3391 | **198.73** | 407 956.8 | 4 213 900.9 | 14551 | 649 |
| `NNjjB` | **W** | 158.10 | 741 869.6 | 22 017 470.1 | 26461 | 3391 | **281.25** | 379 892.4 | 1 824 508.7 | 13550 | 281† |
| `NNjjBj` | **W** | 158.10 | 741 869.6 | 22 017 470.1 | 26461 | 3391 | **377.13** | 142 536.8 | 142 844.1 | 5084 | 22† |
| `NNkin` | **Z** | 44.69 | 68 473.3 | 2 347 365.1 | 30788 | 3423 | **54.46** | 39 807.8 | 534 209.0 | 17899 | 779 |
| `NNj` | **Z** | 44.69 | 68 473.3 | 2 347 365.1 | 30788 | 3423 | **48.44** | 51 321.6 | 1 122 593.3 | 23076 | 1637 |
| `NNjB` | **Z** | 44.69 | 68 473.3 | 2 347 365.1 | 30788 | 3423 | **59.20** | 32 223.9 | 296 249.4 | 14489 | 432 |
| `NNjj` | **Z** | 44.69 | 68 473.3 | 2 347 365.1 | 30788 | 3423 | **57.75** | 40 550.7 | 493 063.3 | 18233 | 719 |
| `NNjjB` | **Z** | 44.69 | 68 473.3 | 2 347 365.1 | 30788 | 3423 | **76.84** | 31 941.5 | 172 812.2 | 14362 | 252 |
| `NNjjBj` | **Z** | 44.69 | 68 473.3 | 2 347 365.1 | 30788 | 3423 | **108.93** | 16 633.5 | 23 315.9 | 7479 | 34 |

`NNjjBj` has the best final `S/√B` in every channel (H: 155.10, W: 377.13,
Z: 108.93), and the ranking broadly tracks AUC: the boson-aware, both-jets nets
(`NNjjB`, `NNjjBj`) consistently beat the single-jet nets (`NNkin`, `NNj`, `NNjB`).
The gain from `_0` (no cut) to the best net's final point is largest, in relative
terms, for **H** (33.10 → 155.10, ×4.7) even though H's baseline background is
already the smallest of the three in absolute terms (`B₀`=80.4 pb-equivalent
yield) — `S/√B` has no `+S` term capping the payoff of driving `B` further toward
zero, and H's NN can push `B` from 80.4 down to 0.2 while keeping a meaningful
fraction of signal. **W** (158.10 → 377.13, ×2.4) and **Z** (44.69 → 108.93, ×2.4)
gain less in relative terms because their much larger baseline background
(`B₀` in the millions) can't be suppressed to quite the same relative degree
without eating unacceptably far into signal too.

### Z_only_table: value of the jet constituents (Z channel)

Added 2026-09-02: `Kj` is a Z-channel-only, feedforward net on just the
leading jet's own four-momentum `(|η|, m, φ, pT)` — no constituents, no
second jet, no event-level kinematics. It is trained and evaluated exactly
like the other nets above and inserted just above `NNj` (which consumes the
same jet plus its ≤80 constituents), isolating what the constituent
information is actually worth. This Z-channel-only comparison is referred to
as the `Z_only_table` from here on.

Extended 2026-09-03 with `KjB` and `KjjBj` — the same feedforward,
no-constituents idea applied one level up the family: `KjB` uses exactly
`NNjB`'s jet-level + boson scalars (`|η|, m, pT, φ` of jet 0 plus
`bosonPt, bosonEta, bosonPhi`, 7 inputs) with no constituent branch, and
`KjjBj` uses exactly `NNjjBj`'s scalars (both leading jets' four-momenta,
the 3rd jet's four-momentum zero-padded when absent, and the boson scalars,
15 inputs total) with no constituent branches at all. Each is inserted
directly above its DeepSets counterpart (`KjB` above `NNjB`, `KjjBj` above
`NNjjBj`), trained for both Z and W (not H).

Narrowed 2026-09-07 to just the constituent-ablation family (`Kj`/`NNj`,
`KjB`/`NNjB`, `KjjBj`/`NNjjBj`) — `NNkin`, `NNjj`, `NNjjB` dropped from this
table since they duplicate the main Results table above and add nothing to
the constituents comparison; `S/(S+B)` (purity) and `eff_sig`/`eff_bkg`
added as context columns.

The `_0` unoptimized baseline (full preselected sample, no NN cut,
`eff_sig=eff_bkg=1`) is identical for every net in the Z channel, so it is
given once here rather than repeated per row: S/√B₀ = 44.69,
S₀ = 68 473.3, B₀ = 2 347 365.1, raw S₀ = 30788, raw B₀ = 3423.

| Net | S/√B | S/(S+B) | eff_sig | eff_bkg | S | B | raw S | raw B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `Kj` | **46.93** | 0.036 | 0.8694 | 0.6857 | 59 532.7 | 1 609 484.7 | 26768 | 2347 |
| `NNj` | **48.44** | 0.044 | 0.7495 | 0.4782 | 51 321.6 | 1 122 593.3 | 23076 | 1637 |
| `KjB` | **58.48** | 0.069 | 0.6746 | 0.2658 | 46 195.2 | 624 043.9 | 20771 | 910 |
| `NNjB` | **59.20** (was 57.48) | 0.098 | 0.4706 | 0.1262 | 32 223.9 | 296 249.4 | 14489 | 432 |
| `KjjBj` | **100.35** | 0.377 | 0.2427 | 0.0117 | 16 620.1 | 27 430.5 | 7473 | 40 |
| `NNjjBj` | **108.93** (was 102.09) | 0.416 | 0.2429 | 0.0099 | 16 633.5 | 23 315.9 | 7479 | 34 |

`S/√B`, `S`, `B`, `raw S`, `raw B` are all at the ROC operating point
maximizing `S/√B` (the only optimization objective); `S/(S+B)` (purity) and
`eff_sig`/`eff_bkg` are shown for context only. `raw S`/`raw B` are the
raw MC event counts actually passing the cut in the signal region (Herwig
test split). `†` on `raw B`
marks a row where the `raw B ≥ 10` statistical-validity floor changed the
unconstrained optimum — the floor is a check on the *total* raw background
count passing the chosen operating point, not a per-bin requirement on any
downstream histogram (e.g. the pull-angle plots below).

**Added jet-φ (2026-09-03):** `NNjB`, `NNjj`, `NNjjB`, `NNjjBj` were extended
with each leading jet's own azimuthal angle `jetPhi` as a jet-level scalar
(`NNj` and `Kj` deliberately excluded — `Kj` already had it, `NNj` stays a
minimal single-jet, no-φ baseline) and retrained from scratch. The
per-constituent Δη convention was investigated at the same time
(`jcsDEta`/`jcsDEtaRaw`, sign-flipped for single-jet nets vs. raw for
multi-jet nets, matching the sign convention used when filling the
`jetShape`/`jetShapeRaw` diagnostic histograms) and confirmed to already be
exactly as intended — no change made there. The `(was ...)` values above are
each net's pre-φ `S/√B` from the previous table, so the retraining
effect is entirely attributable to φ. All four improve: `NNjB` sees the
largest relative gain (57.48→59.20), and for
`NNjjB`/`NNjjBj` the extra discriminating power also relaxes the raw-B floor
constraint — both previously carried a `†` (315→252 raw B for `NNjjB`,
41→34 for `NNjjBj`) and no longer need it, meaning the unconstrained
`S/√B`-optimal cut now naturally lands above the statistical-validity floor.

Constituents clearly help: `Kj` (no constituents) trails `NNj` (same leading
jet, plus its ≤80 constituents) on final `S/√B` (46.93 vs 48.44) despite
starting from the same jet four-momentum — the constituent information is
worth a real, if modest, gain.

The same pattern holds one level up: `KjB` (58.48) trails `NNjB` (59.20), and
`KjjBj` (100.35) trails `NNjjBj` (108.93) — small, consistent gaps, since
`KjB`/`KjjBj` already start from richer inputs (boson scalars, and for
`KjjBj` a second and third jet too) than `Kj`/`NNj` do, so there is less
headroom left for the constituents to add on top.

### Z_only_table (Pythia transfer, no retraining)

Added 2026-09-07: `summer26/pythia_transfer.py` supplements the
`Z_only_table` above with a true generator-transfer test. For each net, the
exact NN-score threshold that produced its tabulated Herwig `S/√B`-optimal
`(eff_sig, eff_bkg)` operating point is recovered (`roc_data_{channel}.npz`
only stores a decimated ROC curve, not thresholds, so the cut is recovered
by re-scoring the Herwig test split at full resolution and matching the
tabulated operating point — same technique as the pull-vector-angle section
below), then that *same, unchanged* score cut is applied to the full
MG5+Pythia sample (never seen in training, and never used to pick the cut).
No retraining and no re-optimization on Pythia: this tests whether the
working point chosen on Herwig alone still holds up when the same events
are generated by a different generator.

| Net | S/√B | S/(S+B) | eff_sig | eff_bkg | S | B | raw S | raw B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `Kj` | **42.54** | 0.031 | 0.8373 | 0.6203 | 57 215.1 | 1 809 397.0 | 174018 | 18888 |
| `NNj` | **43.30** | 0.039 | 0.6772 | 0.3914 | 46 270.8 | 1 141 698.1 | 140731 | 11918 |
| `KjB` | **54.53** | 0.063 | 0.6461 | 0.2247 | 44 151.7 | 655 628.6 | 134286 | 6844 |
| `NNjB` | **54.31** | 0.093 | 0.4213 | 0.0963 | 28 788.7 | 280 970.0 | 87560 | 2933 |
| `KjjBj` | **83.64** | 0.340 | 0.1983 | 0.0090 | 13 551.4 | 26 248.1 | 41216 | 274 |
| `NNjjBj` | **88.56** | 0.381 | 0.1869 | 0.0071 | 12 769.2 | 20 787.8 | 38837 | 217 |

Regenerated 2026-09-27 for the pythia-dipole-recoil samples. The old buggy
Pythia gave 44.33 / 45.24 / 56.54 / 56.53 / 71.39 / 73.21 for the same rows.

Every net loses significance under generator transfer, and the tight-cut
nets still lose more than the loose-cut ones. `S/√B` drops ~7–11% for
`Kj`/`NNj`/`KjB`/`NNjB` and ~17–19% for `KjjBj`/`NNjjBj`.

Compared with the buggy Pythia (~3–7% and ~29–33%), the gap has narrowed
from both ends:
- **Loose-cut nets lose a few % more.** Part of this is an overall
  normalization effect: fewer dipole-recoil Pythia events pass the
  preselection (e.g. VBFZ σ_eff 0.240 → 0.226 pb), so both `S` and `B` shrink.
- **Tight-cut 3rd-jet nets recover a lot.** `NNjjBj` goes 73.21 → 88.56 and
  `eff_sig` 0.159 → 0.187, i.e. the corrected Pythia's 3rd-jet radiation is
  much closer to Herwig's.

The ranking is essentially preserved (`NNjjBj` best, `Kj` worst). `KjB`/`NNjB`
still swap, by a small margin (54.53 vs 54.31). The two feedforward/DeepSets
pairs remain close under transfer, just as they were on Herwig.

### Z_only_table (performance & resilience, Herwig vs Pythia)

Added 2026-09-14: recasts the two `Z_only_table` sections above with a
pure efficiency-ratio figure of merit, **performance** = `eps_sig/sqrt(eps_bkg)`
— unlike `S/√B`, independent of cross-section normalization, so directly
comparable across generators — and introduces **resilience**, the symmetric
sensitivity of a fixed tagger's `(eps_sig, eps_bkg)` to swapping Herwig for
Pythia at the identical Herwig-chosen cut (same transferred-cut convention as
the Pythia-transfer table above — no re-optimization on Pythia):

`zeta(Herwig,Pythia) = { 0.5 * [ ((effS_hw-effS_py)/(effS_hw+effS_py))^2 + ((effB_hw-effB_py)/(effB_hw+effB_py))^2 ] }^-1/2`

Larger `zeta` = more resilient (less sensitive) to the generator swap. The
"no tagging" baseline performance is exactly `1.00` for either generator,
since `eps_sig=eps_bkg=1` at the no-NN-cut point on the same
`mjj>400 GeV` + `|ΔY_jj|>2.5`-preselected sample used throughout — stated
once here rather than as a repeated table column.

**Operating point = maximum Herwig performance (changed 2026-11-15).** Each
net's Herwig operating point is the point of its Herwig test-split ROC grid
(`roc_data_Z.npz`, the same grid used everywhere above) with the largest
`performance`, subject to a raised `raw B >= 20` floor (up from the
`raw B >= 10` used elsewhere in this document). Before 2026-11-15 this table
reused the `S/√B`-optimal point. Both choices give the same point: for a
fixed sample `S/√B = (σ_S·√L/√σ_B) × eps_sig/sqrt(eps_bkg)`, a constant
times `performance`, so both peak at the same ROC point. Checked for all 6
nets (`same_sqrtB = 1` in the output file): identical grid points and cuts.
The floor binds for none of them (Herwig raw B 34–2347, no `†`), and with no
floor at all the same points are selected. Rerunning on 2026-11-15
reproduced every number below exactly; only the selection criterion and the
column order (Herwig next to Pythia) changed.

**Same cut for Herwig and Pythia.** For each net, the NN-score cut `c`
defining the Herwig point is recovered by re-scoring the Herwig test split at
full resolution (pass = `score >= c`, the `roc_curve` convention). That one
`c` is applied to the re-scored Herwig test split (`eff_sig`/`eff_bkg (Hw)`)
and to the net's inference on the full MG5+Pythia sample
(`eff_sig`/`eff_bkg (Py)`), and ζ is computed from these four efficiencies.
Applying `c` to the Herwig test split reproduces the chosen grid point
exactly (distance 0 for all 6 nets). Cuts (`cut` column of the output file):
`Kj` 0.3638, `NNj` 0.4420, `KjB` 0.5422, `NNjB` 0.6633, `KjjBj` 0.9111,
`NNjjBj` 0.9101.

**Statistical errors:** `eff_sig`/`eff_bkg` are each a binomial proportion
`k/N` (MC events passing the fixed cut, out of the total signal/background
population in that generator's sample), so their statistical error is the
standard normal-approximation binomial error `sigma_eff = sqrt(eff*(1-eff)/N)`
— the same formula used for efficiency plots throughout HEP. This is only
approximate for the smallest-count rows (e.g. `NNjjBj`'s Herwig `raw B ≈ 34`,
where the true binomial interval is mildly asymmetric), but adequate here.
`performance` and `resilience` are then deterministic functions of these four
(independent) efficiencies, so their errors are obtained by standard linear
(delta-method) error propagation, `sigma_f^2 = sum_i (df/dxi)^2 * sigma_xi^2`;
`resilience`'s partial derivatives are unwieldy to write out by hand (nested
ratios inside a `-1/2` power over 4 variables), so they're evaluated
numerically via central finite differences rather than derived analytically.

| Net | eff_sig (Hw) | eff_sig (Py) | eff_bkg (Hw) | eff_bkg (Py) | Performance (Hw) | Performance (Py) | Resilience ζ |
|---|---:|---:|---:|---:|---:|---:|---:|
| `Kj` | 0.8694 ± 0.0019 | 0.8373 ± 0.0008 | 0.6857 ± 0.0079 | 0.6203 ± 0.0028 | **1.050 ± 0.007** | **1.063 ± 0.003** | 26.44 ± 2.87 |
| `NNj` | 0.7495 ± 0.0025 | 0.6772 ± 0.0010 | 0.4782 ± 0.0085 | 0.3914 ± 0.0028 | **1.084 ± 0.010** | **1.082 ± 0.004** | 12.62 ± 0.96 |
| `KjB` | 0.6746 ± 0.0027 | 0.6461 ± 0.0010 | 0.2658 ± 0.0076 | 0.2247 ± 0.0024 | **1.308 ± 0.019** | **1.363 ± 0.008** | 16.35 ± 2.76 |
| `NNjB` | 0.4706 ± 0.0028 | 0.4213 ± 0.0011 | 0.1262 ± 0.0057 | 0.0963 ± 0.0017 | **1.325 ± 0.031** | **1.358 ± 0.012** | 9.74 ± 1.47 |
| `KjjBj` | 0.2427 ± 0.0024 | 0.1983 ± 0.0009 | 0.0117 ± 0.0018 | 0.0090 ± 0.0005 | **2.245 ± 0.178** | **2.091 ± 0.064** | 8.60 ± 3.43 |
| `NNjjBj` | 0.2429 ± 0.0024 | 0.1869 ± 0.0009 | 0.0099 ± 0.0017 | 0.0071 ± 0.0005 | **2.437 ± 0.209** | **2.214 ± 0.076** | 6.74 ± 2.25 |

Pythia columns and ζ regenerated 2026-09-27 for the pythia-dipole-recoil
samples. The `KjB` Herwig point moved by one ROC grid step (see the
changelog near the top). Old buggy-Pythia ζ: 22.50 / 10.27 / 15.36 / 8.18 /
7.72 / 6.07.

The errors show which `performance` shifts are statistically meaningful
and which are noise:
- `KjB`'s Herwig→Pythia increase (1.308±0.019 → 1.363±0.008) is significant.
- `KjjBj`'s decrease (2.245±0.178 → 2.091±0.064) is now within ~1σ. It was
  2.245 → 1.731 with the buggy Pythia. The tighter
cuts carry visibly larger relative errors (`NNjjBj` performance is ~9%
uncertain vs `Kj`'s <1%) — a direct consequence of their much smaller raw
background counts.

Unlike `S/√B`, which drops under generator transfer for every net (see
above), this pure efficiency-ratio `performance` metric *increases* or stays
flat from Herwig to Pythia for the four loosest-cut nets
(`Kj`/`NNj`/`KjB`/`NNjB`, e.g. `KjB` 1.308→1.363, `NNj` 1.084→1.082). It
only *decreases* for the two tightest-cut ones (`KjjBj`/`NNjjBj`, e.g.
`NNjjBj` 2.437→2.214, much milder than the 2.437→1.775 seen with the buggy
Pythia) — comparing `eps_sig` and
`eps_bkg` on equal footing (rather than weighting them by each generator's
own total cross section) reveals a more nuanced, cut-dependent picture than
the uniform decline seen in absolute `S/√B`.

Resilience is highest for the simplest, loosest-cut net (`Kj`, ζ=26.44) and
lowest for the richest, tightest-cut one (`NNjjBj`, ζ=6.74), consistent with
the larger relative `S/√B` drops already seen for `KjjBj`/`NNjjBj` in the
Pythia-transfer table above. Every ζ is higher than with the buggy Pythia,
so the corrected shower sits closer to Herwig. Within each
constituent-ablation pair, the feedforward net (no constituents) is
consistently *more* resilient than its DeepSets counterpart — `Kj` 26.44 vs
`NNj` 12.62, `KjB` 16.35 vs `NNjB` 9.74, `KjjBj` 8.60 vs `NNjjBj` 6.74
(the last pair within errors) — suggesting the per-constituent substructure
information that helps performance on Herwig is also the part most sensitive
to the Herwig/Pythia parton-shower difference.

Output: `figs/summer26/Z_performance_resilience.txt`, written by
`summer26/pythia_transfer.py` (which also writes `pythia_transfer.txt`). It
re-scores with every net's saved model, so run it on a GPU node:
`condor_submit -name bigbird24.cern.ch summer26/pythia_transfer_gpu.sub`
(log in `summer26/logs/`; the job also prints the operating-point check
table).

### W_only_table: value of the jet constituents (W channel)

Added 2026-09-03: the W-channel analog of the `Z_only_table` above, using
`KjB`/`KjjBj`'s `_W` variants.

`Kj_W` added 2026-09-07 (trained via `make_variants.py`'s `VARIANTS` list,
same mechanism as `KjB_W`/`KjjBj_W`) so this table can mirror the Z-channel
row set exactly; narrowed at the same time to just the constituent-ablation
family (`Kj`/`NNj`, `KjB`/`NNjB`, `KjjBj`/`NNjjBj`) — `NNkin`, `NNjj`,
`NNjjB` dropped, matching `Z_only_table`; `S/(S+B)` (purity) and
`eff_sig`/`eff_bkg` added as context columns.

Unoptimized baseline (`_0`, identical for every net in the W channel):
S/√B₀ = 158.10, S₀ = 741 869.6, B₀ = 22 017 470.1,
raw S₀ = 26461, raw B₀ = 3391.

| Net | S/√B | S/(S+B) | eff_sig | eff_bkg | S | B | raw S | raw B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `Kj` | **163.74** | 0.041 | 0.8531 | 0.6786 | 632 892.3 | 14 940 194.2 | 22574 | 2301 |
| `NNj` | **169.62** | 0.044 | 0.8438 | 0.6187 | 626 023.4 | 13 622 132.8 | 22329 | 2098 |
| `KjB` | **210.33** | 0.085 | 0.6397 | 0.2312 | 474 543.1 | 5 090 444.3 | 16926 | 784 |
| `NNjB` | **218.54** | 0.095 | 0.6140 | 0.1973 | 455 478.4 | 4 343 759.2 | 16246 | 669 |
| `KjjBj` | **346.59** | 0.497 | 0.1641 | 0.0056 | 121 733.8 | 123 365.4 | 4342 | 19† |
| `NNjjBj` | **377.13** | 0.499 | 0.1921 | 0.0065 | 142 536.8 | 142 844.1 | 5084 | 22† |

Same columns/objective/`†` convention as the `Z_only_table`. Same
qualitative pattern as Z: `KjB` (210.33) trails `NNjB` (218.54) and `KjjBj`
(346.59) trails `NNjjBj` (377.13) — constituents help in both channels, by a
similar relative margin. `Kj` (163.74) trailing `NNj` (169.62) and `KjjBj`
(346.59) trailing `NNjjBj` (377.13) hold the same relative pattern as Z too.

### W_only_table (Pythia transfer, no retraining)

Added 2026-09-07: same generator-transfer test as the `Z_only_table (Pythia
transfer...)` table above — each net's Herwig-chosen score cut, applied
unchanged to the full MG5+Pythia sample.

| Net | S/√B | S/(S+B) | eff_sig | eff_bkg | S | B | raw S | raw B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `Kj` | **144.58** | 0.035 | 0.8093 | 0.6031 | 579 071.8 | 16 042 120.0 | 142898 | 17316 |
| `NNj` | **147.39** | 0.037 | 0.7813 | 0.5409 | 559 045.1 | 14 387 510.0 | 137956 | 15530 |
| `KjB` | **183.98** | 0.073 | 0.6032 | 0.2069 | 431 594.8 | 5 503 014.1 | 106505 | 5940 |
| `NNjB` | **185.24** | 0.079 | 0.5614 | 0.1768 | 401 708.8 | 4 702 575.7 | 99130 | 5076 |
| `KjjBj` | **257.88** | 0.424 | 0.1265 | 0.0046 | 90 521.2 | 123 215.6 | 22338 | 133 |
| `NNjjBj` | **291.71** | 0.467 | 0.1359 | 0.0042 | 97 264.3 | 111 172.0 | 24002 | 120 |

Regenerated 2026-09-27 for the pythia-dipole-recoil samples. The old buggy
Pythia gave 154.07 / 156.67 / 194.17 / 196.61 / 208.70 / 221.41.

Same pattern as Z: every net loses significance under transfer, more so for
the tight-cut nets (`KjjBj`/`NNjjBj`, ~23–26% drop, down from ~40% with the
buggy Pythia) than the loose-cut ones (`Kj`/`NNj`, ~12–13%, up from ~6–8%,
partly from the lower Pythia σ_eff after preselection). Unlike Z's
`KjB`/`NNjB` flip, W's ranking is fully preserved under transfer (no swaps).

Output: `figs/summer26/pythia_transfer.txt`.

### Pull-vector angle inside the S/√B-optimal signal region

Added 2026-09-02 (`summer26/spva_signal_region.py`): for `NNj`, `NNjj`, and
`NNjjBj`, the leading (true pT-highest, via `leadJetIndex` — not the
randomized jet slot 0 used for NN training) jet's pull-vector angle
`|θ_s|` (`jetSPVA`), shown separately for VBF and QCD Herwig events that pass
each net's own S/√B-optimal score cut from the table above. The cut is
recovered by re-scoring the Herwig test split with each net's saved model and
matching the resulting ROC curve to the already-tabulated `(eff_sig, eff_bkg)`
operating point; the resulting passing raw `S`/`B` counts reproduce the `raw S`/
`raw B` columns above exactly (23076/1637 for `NNj`, 18233/719 for `NNjj`,
7479/34 for `NNjjBj`), confirming the recovered threshold is the same
operating point.

Absolute normalization (updated 2026-09-02): each histogram's raw-MC passing-
event shape is rescaled so its integral equals `analyze_one()`'s BR-rescaled
`S` (VBF) or `B` (QCD) at `L = 300 fb⁻¹` — the same totals in the table above,
not a separate re-derivation — giving an `Events / rad` y-axis instead of a
unit-normalized shape comparison. Plotted log-`y`: `S` and `B` differ by up to
~22× in absolute scale across the three nets' cuts, which would otherwise
flatten the signal curve to near-zero on a linear axis.

Caught in the process: `jetSPVA` in the friend h5 files is the *signed* pull
angle in `[-π, π]` (~50/50 split), not `|θ_s|` directly — unlike
`spva_plots.py`'s `hLeadJetSPVA` ROOT histograms, which already take `abs()`
at C++ fill time. `spva_signal_region.py`'s first pass skipped that `abs()`,
silently binning only the ~50% of events with a positive raw angle into
`[0, π]` (the rest fell outside the binning range and were dropped from the
histogram, though the displayed absolute total was unaffected — see the
script's `load_leadjet_spva()` docstring). Fixed; raw histogrammed counts now
match the passing `S`/`B` counts almost exactly (23068/1637, 19200/806,
7694/41 vs. the tabulated 23076/1637, 19208/806, 7697/41 — the handful of
missing events are the genuine `-99` no-pull-vector sentinel, NC=1 jets).

- **NNj** (ε_sig=0.750, ε_bkg=0.478, S=51 321.6, B=1 122 593.3): good raw-MC
  stats on both sides. VBF stays peaked toward small `|θ_s|`; QCD is flatter
  with a mild excess at large angle — the shape difference the pull-angle
  variable is designed to capture survives a loose cut, even though QCD
  outnumbers VBF by ×22 in absolute yield.
- **NNjj** (S=40 550.7, B=493 063.3) and **NNjjBj** (S=16 633.5, B=23 315.9):
  same qualitative VBF/QCD shape difference, though `NNjjBj`'s very tight cut
  (ε_bkg=0.010, only 34 raw background events) makes the QCD histogram
  statistics-limited — visible as large error bars rather than a smooth curve.

Every plot's title reports only the expected event yield for VBF and QCD at
`L = 300 fb⁻¹`, with the raw MC count passing the cut in parenthesis (no
efficiencies or other text); a black `VBF + QCD` step line adds the two
absolute-yield histograms bin-by-bin, showing the total shape an inclusive
counting experiment in the signal region would actually see.

Extended 2026-09-03 to `NNjB` across all three of its channel variants:
**Z** (ε_sig=0.471, ε_bkg=0.126, S=32 223.9, B=296 249.4, raw S/B=14489/432)
and **W** (ε_sig=0.614, ε_bkg=0.197, S=455 478.4, B=4 343 759.2,
raw S/B=16246/669) each get their own plot from that channel's own model and
own S/√B-optimal cut, and a third, **V = W+Z**, adds the two channels'
already-absolute (`L = 300 fb⁻¹`) histograms bin-by-bin (VBF=487 702.3,
QCD=4 640 008.6) — the combined yield if each channel is run inclusively at
its own best `NNjB` working point. Errors on the `V` histogram are added in
quadrature; raw counts and yields are added linearly.

Output: `figs/summer26/jetSPVA_SR_{NNj,NNjj,NNjjBj}_Z.pdf`,
`figs/summer26/jetSPVA_SR_NNjB_{Z,W,V}.pdf`.

---

## Event-level NN score columns in the h5 files (per-channel, 2026-07-08)

Each production h5 file carries per-event sigmoid scores written back by the
`save_scores.py` script of the corresponding net directory (HTCondor GPU jobs,
`save_scores.sub`):

| Column | Net | Inputs |
|--------|-----|--------|
| `NNj_jet0` | NNj | leading jet scalars + constituents (score of jet 0) |
| `NNj_jet1` | NNj | same net applied to the subleading jet |
| `NNjB` | NNjB | leading jet + constituents + boson scalars (one score/event) |
| `NNjjBj` | NNjjBj | jets 0+1 with constituents, 3rd-jet 4-mom + constituents, boson scalars (one score/event) |

**Per-channel convention:** every net is trained on its own channel's Herwig
samples and applied only to that channel's 4 files — Herwig **and** MG5+Pythia
(generator-transfer inference is intended; cross-process inference is not).
Concretely: `NNjB_H/best_model_H.pt` scores `QCDHjj_*` / `VBFH_*`, the base
(Z) dirs score the Z files, `_W` dirs the W files. The `_H` nets keep the
`bosonM` input; the base/W nets exclude it (Herwig on-shell generation
artifact, commit 9ce67b6).

This replaces the original NNj setup where the **Z-trained** net scored all
12 files; the `NNj_jet0`/`NNj_jet1` columns were rewritten per-channel on
2026-07-08 (H- and W-file scores changed, Z unchanged). Validation: for every
net/channel, the AUC recomputed from the stored scores on the full MG5+Pythia
files reproduces `auc_mg5` in the dir's `roc_data_{P}.npz` to 4 decimals.

`plots.py` draws the score distributions for all three event-level columns
in the same two variants as the mjj figures:
`figs/summer26/{NNj_jet0,NNjB,NNjjBj}_{H,W,Z}[_abs].pdf`
(normalized / absolute log-y).

It also draws a 2D `NNj_jet0` vs `NNj_jet1` score-correlation histogram
(kWeight-weighted, 50×50 bins over [0,1]²) for the Z channel's Herwig
samples, signal and background separately:
`figs/summer26/NNj_jet0_vs_jet1_{VBFZ_herwig,QCDZjj_herwig}.pdf` (added
2026-09-01).

Writes go through a small retry loop (`write_dataset(s)` in `save_scores.py`):
EOS FUSE occasionally fails HDF5 metadata operations on the GPU nodes with
`bad object header version number`; if a job dies there anyway, delete the
stale column interactively from lxplus and resubmit (creates don't hit it).

---

## Input files

| Full path | Size (MB) |
|-----------|-----------|
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/VBFH_herwig.root` | 14869 |
| `/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged/VBFH/VBFH_mg5_pythia_dipolerecoil.root` | 14655 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/VBFW_herwig.root` | 15304 |
| `/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged/VBFW/VBFW_mg5_pythia_dipolerecoil.root` | 15098 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/VBFZ_herwig.root` | 15371 |
| `/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged/VBFZ/VBFZ_mg5_pythia_dipolerecoil.root` | 15197 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/QCDWjj_herwig.root` | 16220 |
| `/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged/QCDWjj/QCDWjj_mg5_pythia_dipolerecoil.root` | 16368 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/QCDZjj_herwig.root` | 16198 |
| `/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged/QCDZjj/QCDZjj_mg5_pythia_dipolerecoil.root` | 16368 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00030/merged/QCDHjj_herwig.root` | 17538 |
| `/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged/QCDHjj/QCDHjj_mg5_pythia_dipolerecoil.root` | 7168 |

Total: ~180 GB across 12 files (sizes in MiB; Pythia rows are the dipole-recoil files since 2026-09-27).

---

## File map

```
makefriends/C/
├── run.summer26.sh         # full production script for summer26
├── friends/summer26/       # 12x {_herwig,_mg5_pythia}.{friend.root,h5}; pythia_prefix_bug/ = old buggy Pythia
├── figs/summer26/          # diagnostic PDFs
└── summer26/               # analysis code
    ├── plots.py            # diagnostic plots (adapt from summer25/plots.py)
    ├── NNkin/
    │   ├── dataset.py      # loads VBFZ_herwig + QCDZjj_herwig h5; 12-feature kinematic matrix
    │   ├── model.py        # KinNN: 4-layer MLP (12→128→64→32→1)
    │   ├── train.py        # training loop, early stopping, loss curve
    │   ├── inference.py    # ROC on Herwig test + MG5 transfer test
    │   ├── style.py        # matplotlib style helpers
    │   ├── train_gpu.sh    # LCG_106_cuda wrapper
    │   └── train_gpu.sub   # HTCondor submission
    ├── NNj/
    │   ├── dataset.py      # jet scalars (3: |η|,m,pT) + constituents (3: Δη,ΔΦ,pT/pT_jet)
    │   ├── model.py        # JetNN: phi(3→64→64) + masked pool + rho(67→128→64→1)
    │   ├── train.py        # training loop
    │   ├── inference.py    # ROC on Herwig test + MG5 transfer test
    │   ├── save_scores.py  # writes NNj_jet0/NNj_jet1 sigmoid scores into this channel's 4 h5 files
    │   ├── style.py        # matplotlib style helpers
    │   ├── train_gpu.sh    # LCG_106_cuda wrapper
    │   └── train_gpu.sub   # HTCondor submission
    ├── NNjB/
    │   ├── dataset.py      # jet-level scalars (7: |η|,m,pT + boson pT,η,φ,M) + constituents (3: Δη,ΔΦ,pT/pT_jet)
    │   ├── model.py        # JetNN: phi(3→64→64) + masked pool + rho(71→128→64→1)
    │   ├── train.py        # training loop
    │   ├── inference.py    # ROC on Herwig test + MG5 transfer test
    │   ├── save_scores.py  # writes NNjB event-level score into this channel's 4 h5 files
    │   ├── style.py        # matplotlib style helpers
    │   ├── train_gpu.sh    # LCG_106_cuda wrapper
    │   └── train_gpu.sub   # HTCondor submission
    ├── NNjj/
    │   ├── dataset.py      # jet scalars (3 per jet × 2 jets) + constituents (3 per constituent × 2 jets)
    │   ├── model.py        # JetNN: shared phi(3→64→64) + masked pool per jet + rho(134→128→64→1)
    │   ├── train.py        # training loop; outputs best_model_Z.pt, scaler_Z.pkl, split_indices_Z.npz
    │   ├── inference.py    # ROC on Herwig test + MG5 transfer test
    │   ├── style.py        # matplotlib style helpers
    │   ├── train_gpu.sh    # LCG_106_cuda wrapper
    │   └── train_gpu.sub   # HTCondor submission
    ├── NNjjB/
    │   ├── dataset.py      # NNjj inputs + boson 4-vector (4: pT,η,φ,M) as event-level scalars
    │   ├── model.py        # JetNN: shared phi(3→64→64) + masked pool per jet + rho(138→128→64→1)
    │   ├── train.py        # training loop
    │   ├── inference.py    # ROC on Herwig test + MG5 transfer test
    │   ├── style.py        # matplotlib style helpers
    │   ├── train_gpu.sh    # LCG_106_cuda wrapper
    │   └── train_gpu.sub   # HTCondor submission
    ├── NNjjBj/
    │   ├── dataset.py      # NNjjB inputs + 3rd-jet 4-momentum (4: η,m,φ,pT) + constituents (zeros when absent)
    │   ├── model.py        # JetNN: shared phi(3→64→64) + masked pool per jet (×3) + rho(205→128→64→1)
    │   ├── train.py        # training loop
    │   ├── inference.py    # ROC on Herwig test + MG5 transfer test
    │   ├── save_scores.py  # writes NNjjBj event-level score into this channel's 4 h5 files
    │   ├── style.py        # matplotlib style helpers
    │   ├── train_gpu.sh    # LCG_106_cuda wrapper
    │   └── train_gpu.sub   # HTCondor submission
    └── <NN>_H/, <NN>_W/    # per-process variants generated by make_variants.py
                            # (PROCESS="H"/"W", QCD{H,W}jj + VBF{H,W} samples);
                            # save_scores.py copies are hand-maintained per channel
                            # (make_variants.py skips them by design)
```
