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
To submit the base (Z) trainings:

```bash
condor_submit summer26/NNkin/train_gpu.sub
condor_submit summer26/NNj/train_gpu.sub
condor_submit summer26/NNjB/train_gpu.sub
condor_submit summer26/NNjj/train_gpu.sub
condor_submit summer26/NNjjB/train_gpu.sub
```

Or run interactively (falls back to CPU):

```bash
cd summer26/NNkin  && python train.py
cd summer26/NNj    && python train.py
cd summer26/NNjB   && python train.py
cd summer26/NNjj   && python train.py
cd summer26/NNjjB  && python train.py
```

The `_H` / `_W` variant directories (generated by `summer26/make_variants.py`)
train the same architectures on the H and W samples. Exception: `NNjB_H` and
`NNjjB_H` are excluded from generation and maintained manually — they keep
`bosonM`, which the base dirs dropped; the script enforces this invariant on
every run and fails loudly if violated.

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

**Outputs saved by training** (all NNs use `_{W/Z/H}` suffix matching `PROCESS`):

| File | Contents |
|------|----------|
| `best_model_{Z/H}.pt` | Model weights at best validation loss |
| `scaler_{Z/H}.pkl` | Fitted `StandardScaler` dict |
| `split_indices_{Z/H}.npz` | `train_idx`, `val_idx`, `test_idx` |
| `figs/summer26/{NN}_{Z/H}_loss_curve.pdf` | Train/val loss vs epoch |

---

## Step 6 — Run inference

Submit all inference jobs to HTCondor (commands are in `infer.txt`):

```bash
condor_submit summer26/NNkin/infer_gpu.sub
condor_submit summer26/NNkin_H/infer_gpu.sub
condor_submit summer26/NNkin_W/infer_gpu.sub
condor_submit summer26/NNj/infer_gpu.sub
condor_submit summer26/NNj/save_scores.sub
condor_submit summer26/NNj_H/infer_gpu.sub
condor_submit summer26/NNj_W/infer_gpu.sub
condor_submit summer26/NNjB/infer_gpu.sub
condor_submit summer26/NNjB_H/infer_gpu.sub
condor_submit summer26/NNjB_W/infer_gpu.sub
condor_submit summer26/NNjj/infer_gpu.sub
condor_submit summer26/NNjj_H/infer_gpu.sub
condor_submit summer26/NNjj_W/infer_gpu.sub
condor_submit summer26/NNjjB/infer_gpu.sub
condor_submit summer26/NNjjB_H/infer_gpu.sub
condor_submit summer26/NNjjB_W/infer_gpu.sub
```

Each job requests 2 CPUs, 8 GB RAM, 1 GPU (`microcentury` flavour). Each `inference.py`
evaluates on the Herwig test split (training domain) and the full MG5+Pythia dataset
(generator transfer test), producing `figs/summer26/roc_<nn>_{W,Z,H}.pdf`.

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

| Sample | N events | σ `--xs` (pb) | evweight min | evweight max | evweight mean | evweight std | N sel. | kWeight mean | kWeight std |
|--------|----------|--------------|-------------|-------------|--------------|-------------|--------|-------------|------------|
| `VBFH_herwig` | 1 000 000 | 2.97189 | 1.000 | 424.94 | 1.0056 | 0.622 | 255 574 | 2.9735e-06 | 1.4187e-06 |
| `VBFH_mg5_pythia` | 1 000 000 | 2.88990 | 2.8870 | 2.8920 | 2.8900 | 0.0015 | 335 031 | 2.8899e-06 | 1.4679e-09 |
| `VBFW_herwig` | 1 000 000 | 7.494 | 1.000 | 91.604 | 1.0060 | 0.192 | 330 597 | 7.475e-06 | 9.7132e-07 |
| `VBFW_mg5_pythia` | 1 000 000 | 7.20333 | 7.1962 | 7.2120 | 7.2033 | 0.0051 | 369 846 | 7.2033e-06 | 5.0482e-09 |
| `VBFZ_herwig` | 1 000 000 | 1.099 | 1.000 | 46.794 | 1.0183 | 0.194 | 378 566 | 1.1004e-06 | 2.0376e-07 |
| `VBFZ_mg5_pythia` | 1 000 000 | 1.08581 | 1.0846 | 1.0868 | 1.0856 | 0.00058 | 416 673 | 1.0858e-06 | 5.7907e-10 |
| `QCDHjj_herwig` | 1 000 000 | 5.06051 | 1.000 | 12.842 | 1.0002 | 0.026 | 146 109 | 5.0596e-06 | 2.083e-09 |
| `QCDHjj_mg5_pythia` | 409 621 | 4.982 | 4.9818 | 4.9818 | 4.9818 | ~0 | 106 059 | 1.2162e-05 | 0 |
| `QCDWjj_herwig` | 1 000 000 | 1733.3 | 0.00146 | 6.3698 | 0.11363 | 0.011 | 148 885 | 0.0017311 | 0.00010734 |
| `QCDWjj_mg5_pythia` | 1 000 000 | 1646.805 | 1644.5 | 1649.8 | 1646.8 | 1.73 | 200 332 | 0.0016468 | 1.7313e-06 |
| `QCDZjj_herwig` | 1 000 000 | 340.2 | 0.06666 | 19.326 | 0.06679 | 0.025 | 150 848 | 0.00033961 | 6.1153e-06 |
| `QCDZjj_mg5_pythia` | 1 000 000 | 316.363 | 316.00 | 316.72 | 316.37 | 0.232 | 210 073 | 0.00031636 | 2.3246e-07 |

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

Only events passing all three cuts enter the output `events` TTree (and thus contribute to histograms):

1. ≥ 2 anti-kT (R=0.4) jets with pT > 30 GeV and |η| < 3.0
2. Two leading jets with opposite-sign pseudorapidity: η₀ × η₁ < 0
3. Dijet invariant mass: mjj > 0

Summing `kWeight` over all selected events gives `σ_eff` (effective cross-section after selection, in pb).

### Expected events at Run 3 (L = 300 fb⁻¹) and Poisson uncertainty

For a histogram bin of width Δm (e.g. 10 GeV in mjj):

```
N_bin = L × Σ_{i ∈ bin} kWeight_i
```

with L = 300 000 pb⁻¹ = 300 fb⁻¹. The Poisson statistical uncertainty is:

```
δN_bin = √N_bin = √( L × Σ_{i ∈ bin} kWeight_i )
```

These error bars are shown in `figs/summer26/mjj_{W,Z,H}_Run3.pdf`, produced by
`summer26/plots.py`.

---

## NNj purity study — optimal cut on jet-0 score

Script: `summer26/plots.py` (section starting at line 971).  
Selection: SR defined by `|ΔY_jj| > 3`, weights scaled to L = 300 fb⁻¹ (Run 3).  
Generators averaged (Herwig + MG5+Pythia weighted equally per class).

Scanning NNj_jet0 > threshold to maximise S/(S+B), requiring at least
10 expected signal and 10 expected background events (to avoid the MC-statistics
floor where B drops to zero in the tail).

| Channel | Selection | S | B | S/(S+B) | S/√B |
|---------|-----------|--:|--:|--------:|-----:|
| **H** | SR only (`\|dYjj\|>3`) | 174 000 | 85 000 | 0.672 | 597 |
| **H** | **optimal NNj_jet0 > 0.936** | **1 095** | **40** | **0.965** | **174** |
| **H** | **optimal mjj > 1988 GeV** | **2 130** | **199** | **0.915** | **151** |
| **W** | SR only | 415 000 | 24 700 000 | 0.017 | 83 |
| **W** | **optimal NNj_jet0 > 0.974** | **221** | **247** | **0.472** | **14** |
| **W** | **optimal mjj > 1997 GeV** | **23 118** | **102 334** | **0.184** | **72** |
| **Z** | SR only | 67 500 | 4 800 000 | 0.014 | 31 |
| **Z** | **optimal NNj_jet0 > 0.975** | **33** | **51** | **0.392** | **4.6** |
| **Z** | **optimal mjj > 1999 GeV** | **4 056** | **18 449** | **0.180** | **30** |

**Key observations:**

- **H channel** starts at S/(S+B) = 0.67 already (QCD H+jj cross-section is small).
  With score > 0.936 it reaches **0.965**, keeping 0.6% of signal and 0.05% of
  background. S/√B peaks at moderate cuts (score > 0.5, S/√B = 642) rather than
  at the purity-optimal threshold.
- **W and Z channels** start background-dominated (purity ~1–2%). The NN can push
  purity to **0.47** (W) and **0.39** (Z), but only by retaining a tiny signal tail
  (~0.05% of background remains). S/√B peaks at moderate cuts (score > 0.7) rather
  than at the purity-optimal threshold.
- The optimal thresholds for W/Z (~0.974–0.975) sit near the MC-statistics floor;
  the few tens of background events surviving there carry large MC uncertainties.

---

## SR_Run3 optimization — per-channel (mjj, |ΔYjj|, NNj_jet0) cuts

Script: `summer26/plots.py` (SR optimization section, just before the
`jetSPVA_SR_Run3` plots). Scan performed on **Herwig samples only** at
L = 300 fb⁻¹; the resulting cuts are also applied to the MG5+Pythia figures.

Grid: mjj threshold 0–4000 GeV in 10 GeV steps, |ΔYjj| threshold 0–6 in 0.1
steps, NNj_jet0 threshold 0–1 in 0.01 steps.
Objective: maximize S/(S+B), subject to **S > 1000** expected signal events and
**≥ 10 raw MC background events** (guard against the B→0 MC-statistics floor).
The table is regenerated on every `plots.py` run and written to
`figs/summer26/SR_optimization.txt`.

| Channel | mjj cut (GeV) | \|ΔYjj\| cut | NN cut | S | B | S/(S+B) | raw S | raw B |
|---------|--------------:|------------:|-------:|--:|--:|--------:|------:|------:|
| **H** | > 1250 | > 4.9 | > 0.90 | 1 471 | 15 | **0.990** | 1 659 | 10 |
| **W** | > 2050 | > 5.0 | > 0.86 | 4 632 | 5 192 | **0.472** | 2 059 | 10 |
| **Z** | > 1320 | > 4.2 | > 0.92 | 1 275 | 3 158 | **0.288** | 3 799 | 31 |

For reference (Herwig): the previous common SR (`|ΔYjj| > 3`, `mjj > 2` TeV)
gave H 0.937, W 0.194, Z 0.172; the 2-variable (mjj, |ΔYjj|) optimum gave
H 0.975, W 0.449, Z 0.239. Adding the NN cut both raises purity and relaxes
the mjj cut (more signal kept).

**Note:** the H and W optima sit exactly at the raw-B ≥ 10 guard — the surviving
background estimates rest on only 10 MC events each and carry a ~30%
MC-statistical uncertainty.

The `figs/summer26/jetSPVA_SR_Run3_{H,W,Z}_{Herwig,MG5Pythia}.pdf` figures use
these per-channel optimized cuts (shown in each figure title). The companion
`jetSPVA_jet1_SR_Run3_*.pdf` figures show |θs| of the **sub-leading** jet in the
same SR — largely uncorrelated with the jet-0 |θs| (the NNj cut acts on
jet 0 only).

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
    │   ├── save_scores.py  # writes NNj_jet0/NNj_jet1 sigmoid scores into all 12 h5 files
    │   ├── style.py        # matplotlib style helpers
    │   ├── train_gpu.sh    # LCG_106_cuda wrapper
    │   └── train_gpu.sub   # HTCondor submission
    ├── NNjB/
    │   ├── dataset.py      # jet-level scalars (7: |η|,m,pT + boson pT,η,φ,M) + constituents (3: Δη,ΔΦ,pT/pT_jet)
    │   ├── model.py        # JetNN: phi(3→64→64) + masked pool + rho(71→128→64→1)
    │   ├── train.py        # training loop
    │   ├── inference.py    # ROC on Herwig test + MG5 transfer test
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
    └── <NN>_H/, <NN>_W/    # per-process variants generated by make_variants.py
                            # (PROCESS="H"/"W", QCD{H,W}jj + VBF{H,W} samples)
```
