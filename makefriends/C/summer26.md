# summer26 analysis — VBF H→inv vs QCD dijet

Analysis code for the summer26 campaign (MadGraph5+Pythia samples, campaign-00022).
Prerequisite: complete Steps 1–3 from [README.md](README.md) to build `makefriends`
and produce `friends/summer26/*.{friend.root,h5}`.

---

## Campaign details

**Sources**:
- campaign-00023: `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/`
- campaign-00030: `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00030/merged/` (QCDHjj Herwig only)
- campaign-00031: `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00031/merged/` (QCDHjj MG5+Pythia only)

| Process | σ Herwig (pb) | σ MG5+Pythia (pb) | Herwig file | MG5+Pythia file |
|---------|--------------|-------------------|-------------|----------------|
| VBF H→inv | 2.97189 ± 0.00325 | 2.88990 ± 0.00083 | `VBFH_herwig.root` (camp-00023) | `VBFH_mg5_pythia.root` (camp-00023) |
| VBF W→inv | 7.494 ± 0.00632 | 7.20333 ± 0.00198 | `VBFW_herwig.root` (camp-00023) | `VBFW_mg5_pythia.root` (camp-00023) |
| VBF Z→inv | 1.099 ± 0.00126 | 1.08581 ± 0.00030 | `VBFZ_herwig.root` (camp-00023) | `VBFZ_mg5_pythia.root` (camp-00023) |
| QCD H+jj | 5.06051 ± 0.06030 | 4.982 ± 0.0008692 | `QCDHjj_herwig.root` (camp-00030) | `QCDHjj_mg5_pythia.root` (camp-00031) |
| QCD W+jj | 1733.3 ± 1.84 | 1646.805 ± 0.510 | `QCDWjj_herwig.root` (camp-00023) | `QCDWjj_mg5_pythia.root` (camp-00023) |
| QCD Z+jj | 340.2 ± 0.316 | 316.363 ± 0.098 | `QCDZjj_herwig.root` (camp-00023) | `QCDZjj_mg5_pythia.root` (camp-00023) |

All 12 files available. `QCDHjj_mg5_pythia` comes from campaign-00031 (xs = 4.982 ± 0.0008692 pb).

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

**Effect on NN significance** (counter-intuitive): the per-net/channel
significance estimates in "Results" below are *lower* than under the old
inclusive selection, even though signal purity `S/(S+B)` improved
substantially. At the significance-maximizing operating point signal already
dominates background (S≫B), so significance is close to
signal-count-limited (∝√S) rather than background-limited; the
preselection's ~30–40% loss of raw signal (removed as an unavoidable side
effect of removing far more background) costs more significance than the
much larger background rejection recovers. This is a real effect of trading
some peak toy significance for a cleaner, more physically realistic
VBF-topology sample — not a regression.

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
| `friends/summer26/QCDHjj_mg5_pythia.friend.root` | QCD H+jj MG5+Pythia friend tree (from camp-00031, xs=4.982 ± 0.0008692 pb) |
| `friends/summer26/QCDHjj_mg5_pythia.h5` | QCD H+jj MG5+Pythia HDF5 (from camp-00031) |
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

**Jet-level scalars** (6 = 3 jet + 3 boson): |η|, m, pT, boson pT, boson η, boson φ  
**Constituent features** (80 max × 3): Δη, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP**: 3 → 64 → 64 (ReLU, shared across constituents)
- masked-sum pooling → (B, 64)
- **rho MLP**: 70 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3)); 71 in `_H` (bosonM kept)

Trained on **QCDZjj_herwig + VBFZ_herwig**.

### NNjj (`summer26/NNjj/`)

DeepSets over the constituents of **both** the leading (index 0) and sub-leading
(index 1) jets, using only raw low-level features — no pull-vector scalars, no
constituent weight, and signed η (not |η|). The phi MLP is shared across all
constituents from both jets.

**Jet scalars** (3 per jet × 2 jets = 6 total): η (signed), m, pT  
**Constituent features** (80 max × 3 per jet): Δη_raw, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP** (shared): 3 → 64 → 64 (ReLU, applied independently per jet)
- masked-sum pooling per jet → pool₀ (B, 64), pool₁ (B, 64)
- concat [pool₀, pool₁, jet₀ scalars, jet₁ scalars] → (B, 134)
- **rho MLP**: 134 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3))

Trained on **QCDZjj_herwig + VBFZ_herwig**.

### NNjjB (`summer26/NNjjB/`)

Extension of NNjj that adds the generator-boson pT, η, φ as event-level
scalars. `bosonM` is excluded (on-shell generation artifact in the Herwig QCD
Z/W samples — see tuesday.md); only the `_H` variant keeps it (harmless
constant 125).

**Jet scalars** (3 per jet × 2 jets = 6 total): η (signed), m, pT  
**Boson scalars** (3, event-level): boson pT, η, φ  
**Constituent features** (80 max × 3 per jet): Δη_raw, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP** (shared): 3 → 64 → 64 (ReLU, applied independently per jet)
- masked-sum pooling per jet → pool₀ (B, 64), pool₁ (B, 64)
- concat [pool₀, pool₁, jet₀ scalars, jet₁ scalars, boson scalars] → (B, 137)
- **rho MLP**: 137 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3)); 138 in `_H` (bosonM kept)

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

**Jet scalars** (3 per jet × 2 jets = 6 total): η (signed), m, pT
**3rd-jet scalars** (4, event-level): η (signed), m, φ, pT — zeros when absent
**Boson scalars** (3, event-level): boson pT, η, φ
**Constituent features** (80 max × 3 per jet, all 3 jets): Δη_raw, ΔΦ, pT/pT_jet

Architecture:
- **phi MLP** (shared): 3 → 64 → 64 (ReLU, applied independently per jet, all 3 jets)
- masked-sum pooling per jet → pool₀, pool₁, pool₂ (B, 64) each; pool₂ = 0 when no 3rd jet
- concat [pool₀, pool₁, pool₂, jet₀ scalars, jet₁ scalars, jet₂ scalars, boson scalars] → (B, 205)
- **rho MLP**: 205 → 128 → 64 → 1 (BatchNorm + ReLU + Dropout(0.3)); 206 in `_H` (bosonM kept)

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
single-observable baselines).

---

## kWeight factors and Run 3 statistical uncertainty

### Where `evweight` lives

The `evweight` branch is stored in the **original campaign ROOT files on EOS** — not in the friend trees. The TTree in these files is named **`Data`**. The exact file paths per sample are:

| Sample | EOS file |
|--------|----------|
| `VBFH_herwig` | `campaign-00023/merged/VBFH_herwig.root` |
| `VBFH_mg5_pythia` | `campaign-00023/merged/VBFH_mg5_pythia.root` |
| `VBFW_herwig` | `campaign-00023/merged/VBFW_herwig.root` |
| `VBFW_mg5_pythia` | `campaign-00023/merged/VBFW_mg5_pythia.root` |
| `VBFZ_herwig` | `campaign-00023/merged/VBFZ_herwig.root` |
| `VBFZ_mg5_pythia` | `campaign-00023/merged/VBFZ_mg5_pythia.root` |
| `QCDWjj_herwig` | `campaign-00023/merged/QCDWjj_herwig.root` |
| `QCDWjj_mg5_pythia` | `campaign-00023/merged/QCDWjj_mg5_pythia.root` |
| `QCDZjj_herwig` | `campaign-00023/merged/QCDZjj_herwig.root` |
| `QCDZjj_mg5_pythia` | `campaign-00023/merged/QCDZjj_mg5_pythia.root` |
| `QCDHjj_herwig` | `campaign-00030/merged/QCDHjj_herwig.root` |
| `QCDHjj_mg5_pythia` | `campaign-00031/merged/QCDHjj_mg5_pythia.root` |

where `campaign-00023`, `campaign-00030`, `campaign-00031` are all under
`/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/`.

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
unchanged from before.

| Sample | N events | σ `--xs` (pb) | evweight min | evweight max | evweight mean | evweight std | N sel. | kWeight mean | kWeight std |
|--------|----------|--------------|-------------|-------------|--------------|-------------|--------|-------------|------------|
| `VBFH_herwig` | 1 000 000 | 2.97189 | 1.000 | 424.94 | 1.0056 | 0.622 | 146 901 | 2.9668e-06 | 1.3108e-06 |
| `VBFH_mg5_pythia` | 1 000 000 | 2.88990 | 2.8870 | 2.8920 | 2.8900 | 0.0015 | 171 085 | 2.8899e-06 | 1.4685e-09 |
| `VBFW_herwig` | 1 000 000 | 7.494 | 1.000 | 91.604 | 1.0060 | 0.192 | 176 406 | 7.4755e-06 | 9.3813e-07 |
| `VBFW_mg5_pythia` | 1 000 000 | 7.20333 | 7.1962 | 7.2120 | 7.2033 | 0.0051 | 191 379 | 7.2033e-06 | 5.0483e-09 |
| `VBFZ_herwig` | 1 000 000 | 1.099 | 1.000 | 46.794 | 1.0183 | 0.194 | 205 250 | 1.1017e-06 | 2.182e-07 |
| `VBFZ_mg5_pythia` | 1 000 000 | 1.08581 | 1.0846 | 1.0868 | 1.0856 | 0.00058 | 220 628 | 1.0858e-06 | 5.7927e-10 |
| `QCDHjj_herwig` | 1 000 000 | 5.06051 | 1.000 | 12.842 | 1.0002 | 0.026 | 23 340 | 5.0596e-06 | ~0 |
| `QCDHjj_mg5_pythia` | 409 621 | 4.982 | 4.9818 | 4.9818 | 4.9818 | ~0 | 21 818 | 1.2162e-05 | 0 |
| `QCDWjj_herwig` | 1 000 000 | 1733.3 | 0.00146 | 6.3698 | 0.11363 | 0.011 | 22 603 | 0.0017315 | 4.4768e-05 |
| `QCDWjj_mg5_pythia` | 1 000 000 | 1646.805 | 1644.5 | 1649.8 | 1646.8 | 1.73 | 30 340 | 0.0016468 | 1.7338e-06 |
| `QCDZjj_herwig` | 1 000 000 | 340.2 | 0.06666 | 19.326 | 0.06679 | 0.025 | 22 821 | 0.00033969 | 6.7361e-06 |
| `QCDZjj_mg5_pythia` | 1 000 000 | 316.363 | 316.00 | 316.72 | 316.37 | 0.232 | 32 309 | 0.00031636 | 2.3195e-07 |

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
For every (net, channel) combination — the same 18 `roc_data_{H,W,Z}.npz` files used by
the `roc_summary_{H,W,Z}.pdf` figures — the script scans the net's own ROC curve
(signal efficiency `tpr` vs background efficiency `fpr`, no additional mjj/ΔYjj cuts) and
picks the operating point that **maximizes significance `S/√(S+B)`** at **L = 300 fb⁻¹**,
subject to a `raw B ≥ 10` statistical-validity floor (same guard convention as the removed
`SR_optimization.txt` study). `S` and `B` come from each sample's post-selection
`Σ kWeight` (pb, see kWeight section above) × `L` × the ROC efficiency at that point. Both
the Herwig held-out test split (training domain) and the full MG5+Pythia dataset
(independent generator-transfer test) are scanned and reported side by side. Raw output:
`figs/summer26/SB_estimate.txt`.

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
maximizes significance (`significance_final = √factor × significance_unrescaled` at every
point, so the argmax is unchanged — verified numerically by `sb_estimate.py`'s
`verify_rescale_invariance()`). The scan therefore runs once on the un-rescaled yields;
the rescale is applied only to the winning point's final `S`, `B` before display. In other
words, this whole exercise changes the absolute event counts to a physically realistic
decay channel, but says nothing about which NN or working point is best — that ordering
is exactly what it would be without any BR at all.

**Statistical floor caveat:** the `raw B ≥ 10` guard never actually binds for any of the
36 rows below (`floor` column in `SB_estimate.txt` is 0 throughout) — unlike the old
S/(S+B)-purity-optimized table (which sat at this floor almost everywhere), a significance
objective doesn't push toward the tightest possible cut, so the guard is inactive here by
construction rather than by coincidence.

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

### Results (18 nets × 2 generators, best ROC operating point at L = 300 fb⁻¹)

Regenerated 2026-09-01 under the `mjj>400 GeV` + `|ΔY_jj|>2.5` preselection
(see the changelog section near the top of this document); numbers below are
**not** comparable to the pre-2026-09-01 version of this table (S, B, and
significance all shift because the preselected sample is smaller — see the
"Effect on NN significance" note in the changelog).

| Net | Ch | AUC (HW) | ε_sig* (HW) | ε_bkg* (HW) | S (HW) | B (HW) | S/(S+B) (HW) | S/√(S+B) (HW) | raw S (HW) | raw B (HW) | AUC (MG5) | ε_sig* (MG5) | ε_bkg* (MG5) | S (MG5) | B (MG5) | S/(S+B) (MG5) | S/√(S+B) (MG5) | raw S (MG5) | raw B (MG5) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `NNkin` | **H** | 0.8558 | 0.961 | 0.632 | 285.1 | 50.8 | 0.849 | **15.56** | 21171 | 2212 | 0.8037 | 0.934 | 0.627 | 314.6 | 113.3 | 0.735 | **15.21** | 159839 | 13679 |
| `NNj` | **H** | 0.8153 | 0.978 | 0.800 | 290.2 | 64.4 | 0.818 | **15.41** | 21547 | 2802 | 0.7733 | 0.953 | 0.739 | 320.7 | 133.5 | 0.706 | **15.05** | 162978 | 16118 |
| `NNjB` | **H** | 0.8278 | 0.974 | 0.749 | 289.1 | 60.3 | 0.828 | **15.47** | 21464 | 2623 | 0.7825 | 0.951 | 0.716 | 320.0 | 129.3 | 0.712 | **15.10** | 162625 | 15617 |
| `NNjj` | **H** | 0.9019 | 0.959 | 0.473 | 284.6 | 38.0 | 0.882 | **15.84** | 21129 | 1655 | 0.8566 | 0.924 | 0.459 | 311.0 | 82.9 | 0.790 | **15.67** | 158029 | 10011 |
| `NNjjB` | **H** | 0.9104 | 0.949 | 0.396 | 281.6 | 31.9 | 0.898 | **15.90** | 20906 | 1388 | 0.8652 | 0.922 | 0.434 | 310.4 | 78.4 | 0.798 | **15.74** | 157739 | 9470 |
| `NNjjBj` | **H** | 0.9284 | 0.967 | 0.421 | 287.0 | 33.8 | 0.895 | **16.02** | 21307 | 1473 | 0.8692 | 0.913 | 0.428 | 307.3 | 77.3 | 0.799 | **15.67** | 156153 | 9331 |
| `NNkin` | **W** | 0.7241 | 0.748 | 0.429 | 555 175.6 | 9 453 682.2 | 0.055 | **175.48** | 19802 | 1456 | 0.7400 | 0.611 | 0.257 | 473 707.7 | 7 211 371.9 | 0.062 | **170.88** | 116897 | 7784 |
| `NNj` | **W** | 0.6842 | 0.844 | 0.619 | 626 023.4 | 13 622 132.8 | 0.044 | **165.85** | 22329 | 2098 | 0.6886 | 0.799 | 0.554 | 619 296.5 | 15 567 817.8 | 0.038 | **153.93** | 152824 | 16804 |
| `NNjB` | **W** | 0.7821 | 0.641 | 0.231 | 475 187.9 | 5 090 444.3 | 0.085 | **201.42** | 16949 | 784 | 0.7800 | 0.634 | 0.231 | 491 817.6 | 6 486 899.6 | 0.070 | **186.17** | 121366 | 7002 |
| `NNjj` | **W** | 0.7515 | 0.561 | 0.201 | 416 255.5 | 4 434 660.0 | 0.086 | **188.99** | 14847 | 683 | 0.7572 | 0.519 | 0.164 | 402 597.0 | 4 616 426.8 | 0.080 | **179.71** | 99349 | 4983 |
| `NNjjB` | **W** | 0.8412 | 0.477 | 0.070 | 354 127.0 | 1 545 313.4 | 0.186 | **256.95** | 12631 | 238 | 0.8384 | 0.523 | 0.092 | 405 287.8 | 2 581 975.0 | 0.136 | **234.49** | 100013 | 2787 |
| `NNjjBj` | **W** | 0.8755 | 0.465 | 0.047 | 345 015.2 | 1 038 866.2 | 0.249 | **293.28** | 12306 | 160 | 0.8260 | 0.414 | 0.057 | 321 448.7 | 1 602 733.0 | 0.167 | **231.73** | 79324 | 1730 |
| `NNkin` | **Z** | 0.7458 | 0.581 | 0.228 | 39 807.8 | 534 209.0 | 0.069 | **52.54** | 17899 | 779 | 0.7630 | 0.560 | 0.187 | 40 599.2 | 577 654.5 | 0.066 | **51.63** | 123481 | 6030 |
| `NNj` | **Z** | 0.6847 | 0.780 | 0.518 | 53 385.5 | 1 216 542.7 | 0.042 | **47.37** | 24004 | 1774 | 0.7045 | 0.739 | 0.450 | 53 589.3 | 1 391 639.7 | 0.037 | **44.58** | 162990 | 14527 |
| `NNjB` | **Z** | 0.7679 | 0.585 | 0.207 | 40 081.4 | 486 205.6 | 0.076 | **55.25** | 18022 | 709 | 0.7822 | 0.671 | 0.255 | 48 652.8 | 787 928.4 | 0.058 | **53.19** | 147976 | 8225 |
| `NNjj` | **Z** | 0.7609 | 0.624 | 0.235 | 42 719.1 | 552 724.6 | 0.072 | **55.36** | 19208 | 806 | 0.7741 | 0.536 | 0.157 | 38 884.9 | 486 072.8 | 0.074 | **53.67** | 118267 | 5074 |
| `NNjjB` | **Z** | 0.8324 | 0.564 | 0.111 | 38 618.0 | 259 904.0 | 0.129 | **70.68** | 17364 | 379 | 0.8463 | 0.469 | 0.065 | 34 039.8 | 200 406.8 | 0.145 | **70.30** | 103531 | 2092 |
| `NNjjBj` | **Z** | 0.8707 | 0.446 | 0.040 | 30 515.9 | 93 263.7 | 0.247 | **86.74** | 13721 | 136 | 0.8381 | 0.375 | 0.038 | 27 218.8 | 117 734.2 | 0.188 | **71.49** | 82785 | 1229 |

`NNjjBj` has the best significance in 4 of the 6 channel/generator blocks
(H-Herwig 16.02, W-Herwig 293.28, Z-Herwig 86.74, Z-MG5 71.49); `NNjjB` edges
it out in the other 2 (H-MG5: 15.74 vs 15.67, a 0.4% difference; W-MG5:
234.49 vs 231.73, a 1.2% difference). In every channel the ranking broadly
tracks AUC: the boson-aware, both-jets nets (`NNjjB`, `NNjjBj`) consistently
beat the single-jet nets (`NNkin`, `NNj`, `NNjB`). Compared to the old
inclusive-selection table, absolute significance is lower across the board
(e.g. H-Herwig `NNjjBj`: 19.82 → 16.02) even though purity `S/(S+B)`
improved (same cell: 0.852 → 0.895) — see the "Effect on NN significance"
note in the preselection changelog above for why.

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
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/VBFH_mg5_pythia.root` | 15912 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/VBFW_herwig.root` | 15304 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/VBFW_mg5_pythia.root` | 15981 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/VBFZ_herwig.root` | 15371 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/VBFZ_mg5_pythia.root` | 16027 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/QCDWjj_herwig.root` | 16220 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/QCDWjj_mg5_pythia.root` | 16774 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/QCDZjj_herwig.root` | 16198 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged/QCDZjj_mg5_pythia.root` | 16776 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00030/merged/QCDHjj_herwig.root` | 17538 |
| `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00031/merged/QCDHjj_mg5_pythia.root` | 7380 |

Total: ~178 GB across 12 files.

---

## File map

```
makefriends/C/
├── run.summer26.sh         # full production script for summer26
├── friends/summer26/       # 12x {_herwig,_mg5_pythia}.{friend.root,h5}
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
