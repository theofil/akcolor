# Analysis workflow — VBF H→inv vs QCD dijet

Classification of VBF H→invisible signal against QCD dijet background using
jet pull vectors and constituent information.

---

## Environment

One environment covers everything (compilation, ROOT macros, Python training,
inference, and plotting):

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_106_cuda/x86_64-el9-gcc11-opt/setup.sh
```

This is also saved in `init.txt` for convenience:

```bash
source init.txt
```

---

## Step 1 — Build `makefriends`

```bash
make
```

The Makefile detects `root-config` and `fastjet-config` automatically. If
FastJet is not on `PATH` it falls back to the CMS CVMFS installation at
`/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fastjet/3.4.1-...`.

---

## Step 2 — Run event reconstruction (`makefriends.cpp`)

`makefriends` reads a generator-level ROOT file, clusters particles with
anti-kT (R=0.4) via FastJet, computes the **pull vector** for each jet, and
writes a "friend" ROOT tree.

### Quick test (100 events per sample, foreground)

```bash
./run_HtoInv.sh --goFast
# or a custom number of events:
./run_HtoInv.sh --goFast 500
```

### Full production (20 parallel jobs × 50k events, then hadd)

```bash
./run_HtoInv.sh
```

The script:
1. Creates symlinks to the Herwig 13.6 input ROOT files (on apapaefs EOS)
2. Runs `makefriends` in parallel (20 jobs × 50k events each)
3. Merges with `hadd`
4. Moves results to `friends/` and cleans up intermediate files

**Input files** (symlinked automatically):

| Symlink | Sample |
|---------|--------|
| `VBFHtoInv.root` | VBF H→inv, xs = 3.901 pb |
| `QCDHtoInv.root` | ggH+jj (QCD background), xs = 2.26114 pb |

**Output**: `friends/VBFHtoInv.friend.root`, `friends/QCDHtoInv.friend.root`

### Manual single-file invocation

```bash
./makefriends VBFHtoInv.root --genWeight evweight --xs 3.901 --totEve 1000 --output VBFHtoInv.friend.root
```

### Friend ROOT tree structure

Per-jet branches are written with a slot-index suffix (`0` or `1`).
**Slot 0/1 do not correspond to leading/subleading** — they are randomised
(see below).

| Branch | Type | Description |
|--------|------|-------------|
| `jetPt{0,1}`, `jetEta{0,1}`, `jetPhi{0,1}`, `jetM{0,1}` | `Float_t` | Jet 4-vector for slot 0 and slot 1 |
| `jetNC{0,1}` | `Int_t` | Number of constituents |
| `jetPVM{0,1}` | `Float_t` | Pull vector magnitude \|t⃗\| |
| `jetPVA{0,1}` | `Float_t` | Pull vector angle |
| `jetSPVA{0,1}` | `Float_t` | Signed pull vector angle θ_s |
| `jcsPt{0,1}[80]`, `jcsDEta{0,1}[80]`, `jcsDPhi{0,1}[80]`, `jcsM{0,1}[80]`, `jcsW{0,1}[80]` | `Float_t[80]` | Per-constituent features (up to 80, zero-padded) |
| `mjj`, `dYjj`, `dPhijj`, `ptjj` | `Float_t` | Dijet event-level variables |
| `kWeight` | `Float_t` | Event weight (gen weight × xs / total events) |
| `nJets` | `Int_t` | Number of reconstructed jets stored |
| `leadJetIndex` | `Int_t` | Slot index (0 or 1) that holds the true pT-leading jet |

### Jet slot randomisation

The two leading jets are written into slots with a **per-event swap** baked in
at fill time ([makefriends.cpp:337-364](makefriends.cpp)):

```
Line 337: bool swap = (iev % 2 == 0);   // flip for even event numbers
Line 338: o_leadJetIndex = swap ? 1 : 0; // records which slot is truly leading
Line 364: int si = swap ? (1 - i) : i;  // physical jet i → output slot si
```

In 50 % of events the pT-leading jet lands in slot 0; in the other 50 % it
lands in slot 1.  The variable `leadJetIndex` records the true leading slot but
`dataset.py` deliberately ignores it and always reads slot 0 — so the network
sees the leading and subleading jet with equal probability and cannot learn
"leading = signal-like" as a shortcut.

The randomisation is **upstream** in the C++ friend-tree maker; the HDF5 files
and all downstream Python code already have it baked in.

---

## Step 3 — Convert to HDF5 (`root_to_hdf5.py`)

```bash
python root_to_hdf5.py
```

Reads `friends/*.friend.root` with `uproot`, stacks jet pairs into `(N, 2)`
arrays and constituents into `(N, 2, 80)` arrays, and writes compressed HDF5.

**Output**: `friends/VBFHtoInv.h5`, `friends/QCDHtoInv.h5`

### HDF5 dataset layout

Each file contains the following datasets (gzip-compressed, level 4):

| Dataset | Shape | dtype | Description |
|---------|-------|-------|-------------|
| `jetPt` | (N, 2) | float32 | Jet pT for slots 0 and 1 |
| `jetEta` | (N, 2) | float32 | Jet η |
| `jetPhi` | (N, 2) | float32 | Jet φ |
| `jetM` | (N, 2) | float32 | Jet mass |
| `jetSPVA` | (N, 2) | float32 | Signed pull vector angle θ_s |
| `jetPVA` | (N, 2) | float32 | Pull vector angle |
| `jetPVM` | (N, 2) | float32 | Pull vector magnitude \|t⃗\| |
| `jetNC` | (N, 2) | float32 | Constituent count |
| `jcsPt` | (N, 2, 80) | float32 | Constituent pT (zero = absent) |
| `jcsDEta` | (N, 2, 80) | float32 | Constituent Δη (sign-flipped for backward jets) |
| `jcsDPhi` | (N, 2, 80) | float32 | Constituent ΔΦ |
| `jcsM` | (N, 2, 80) | float32 | Constituent mass |
| `jcsW` | (N, 2, 80) | float32 | Pull-vector weight |
| `kWeight` | (N,) | float32 | Event weight |
| `mjj` | (N,) | float32 | Dijet invariant mass |
| `dYjj` | (N,) | float32 | Dijet rapidity separation |
| `dPhijj` | (N,) | float32 | Dijet Δφ |
| `ptjj` | (N,) | float32 | Dijet pT |
| `leadJetIndex` | (N,) | int32 | Slot (0 or 1) of the true pT-leading jet |

The slot ordering in axis-1 of the `(N, 2, …)` arrays mirrors the randomised
slot ordering in the friend ROOT trees.  Training reads only `[:, 0, ...]`
(slot 0), which is the true leading jet in half the events and the subleading
jet in the other half — by design.

---

## Step 4 — Diagnostic plots (`plots.py`)

```bash
python plots.py
```

**Input**: all `friends/*friend.root` files. The six named samples it specifically
looks for are:

| File | Sample |
|------|--------|
| `friends/QCDHtoInv.friend.root` | QCD H |
| `friends/QCDZtoInv.friend.root` | QCD Z |
| `friends/QCDWtoInv.friend.root` | QCD W |
| `friends/VBFHtoInv.friend.root` | VBF H |
| `friends/VBFZtoInv.friend.root` | VBF Z |
| `friends/VBFWtoInv.friend.root` | VBF W |

**Output directory**: `~/www/files/akcolor/`

| Output file | Content |
|-------------|---------|
| `{hname}_{stem}.pdf` | One PDF per TH1F / TH2F / TH2Poly histogram found in each friend file (except `hLeadJetSPVA`, `hLeadJetSPVM`, and `jetSPVA_*` which are handled separately below) |
| `jetSPVA.pdf` | Normalized distributions of the leading-jet `\|θ_s\|` (SPVA angle) for all 6 samples overlaid |
| `jetSPVA_ratio.pdf` | VBF / QCD ratio of `\|θ_s\|` distributions per process (H, Z, W) with RMS / ⟨S²⟩ / KL stats in legend |
| `jetSPVM.pdf` | Normalized distributions of the leading-jet PVM magnitude `\|t⃗\|` for all 6 samples (log-y scale) |
| `jetSPVM_ratio.pdf` | VBF / QCD ratio of `\|t⃗\|` distributions per process |
| `jetNC.pdf` | Normalized distributions of jet number of constituents N_c |
| `jetNC_ratio.pdf` | VBF / QCD ratio of N_c distributions per process |
| `jet_roc.pdf` | ROC curves for 6 variables (`\|θ_s\|`, `\|t⃗\|`, p_T, `\|η\|`, m, N_c) × 3 processes, with AUC in legend |
| `{hname}_ratio2D_{proc}.pdf` | VBF / QCD ratio TH2Poly map for each 2D poly histogram, per process (H, Z, W), with RMS / ⟨S²⟩ / KL stats |
| `index.php` | Auto-generated PHP page listing all PDFs in the output directory |

View plots at:

```
https://theofil.web.cern.ch/theofil/files/akcolor/summary.php
```

---

## Step 5 — Train the neural networks

Both NNs use the same LCG environment. Training was submitted to a GPU node
via HTCondor (`train_gpu.sub`). To re-run on a GPU node:

```bash
condor_submit NNkin/train_gpu.sub
condor_submit NNpull/train_gpu.sub
```

Or run interactively (falls back to CPU):

```bash
cd NNkin  && python train.py
cd NNpull && python train.py
```

### NNkin (`NNkin/`)

Kinematic-only MLP. 12 input features: `dPhijj, dYjj, mjj, ptjj` plus
`|eta|, m, phi, pt` for each of the two leading jets.
Architecture: 12 → 128 → 64 → 32 → 1 (BatchNorm + Dropout at each layer).

### NNpull (`NNpull/`)

DeepSets-style constituent network. Per-jet scalars (6 features: `|eta|, m,
NC, |t⃗|, pt, θ_s`) plus up to 80 constituents × 5 features each
(`dEta, dPhi, m, pt, w`). A shared MLP `phi` processes each constituent;
masked-sum pooling aggregates them; a second MLP `rho` combines with jet
scalars to produce the logit.

**Outputs saved by training** (both NNs):

| File | Contents |
|------|----------|
| `best_model.pt` | Model weights at best validation loss |
| `scaler.pkl` | Fitted `StandardScaler` (or dict of two scalers for NNpull) |
| `split_indices.npz` | `train_idx`, `val_idx`, `test_idx` |
| `loss_curve.pdf` | Train/val loss vs epoch |

---

## Step 6 — Run inference

```bash
source init.txt

cd NNpull && python inference.py   # → roc_nn.pdf
cd NNkin  && python inference.py   # → roc_nnkin.pdf
```

Each script loads `best_model.pt` and `split_indices.npz`, runs the model on
the held-out test set, and plots the NN ROC curve overlaid with
single-variable baselines.

---

## Step 7 — Head-to-head comparison

```bash
cd NNkin && python compare_roc.py   # → compare_roc.pdf
```

Loads both trained models, evaluates on their respective test partitions, and
plots NNpull vs NNkin vs mjj vs dYjj on a single ROC canvas.

---

## Step 8 — Publish plots

Copy output PDFs to the web folder and view via `summary.php`:

```bash
cp NNpull/roc_nn.pdf NNkin/roc_nnkin.pdf NNkin/compare_roc.pdf ~/www/files/akcolor/
```

Then open:
```
https://theofil.web.cern.ch/theofil/files/akcolor/summary.php
```

---

## File map

```
makefriends/C/
├── init.txt              # source this to set up the environment
├── Makefile              # builds makefriends binary
├── makefriends.cpp       # event reconstruction + pull vector (ROOT + FastJet)
├── run_HtoInv.sh         # full production script (parallel jobs + hadd)
├── root_to_hdf5.py       # ROOT → HDF5 converter
├── plots.py              # diagnostic plots → ~/www/files/akcolor/
├── friends/              # friend ROOT files and HDF5 datasets
├── NNkin/
│   ├── dataset.py        # loads HDF5, builds 12-feature kinematic matrix
│   ├── model.py          # KinNN: 4-layer MLP
│   ├── train.py          # training loop, early stopping, loss curve
│   ├── inference.py      # ROC on test set vs single-variable baselines
│   ├── compare_roc.py    # NNkin vs NNpull head-to-head ROC
│   ├── style.py          # matplotlib style helpers
│   └── train_gpu.sub     # HTCondor submission file
└── NNpull/
    ├── dataset.py        # loads HDF5, builds jet scalars + constituent arrays
    ├── model.py          # JetNN: DeepSets (phi MLP + masked pool + rho MLP)
    ├── train.py          # training loop
    ├── inference.py      # ROC on test set vs single-variable baselines
    ├── style.py          # matplotlib style helpers
    └── train_gpu.sub     # HTCondor submission file
```
