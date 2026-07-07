# Analysis workflow — VBF H→inv vs QCD dijet

Classification of VBF H→invisible signal against QCD dijet background using
jet pull vectors and constituent information.

For analysis steps (plots, NN training, inference) see
[summer26.md](summer26.md). Conclusions from the earlier summer25 campaign
are preserved in [studies/](studies/).

---

## Available networks

All networks live in `summer26/<name>/` (Z training), with `_H`/`_W` variant
directories for the H and W processes. All are trained on Herwig samples
(`QCD{Z,W,H}jj_herwig + VBF{Z,W,H}_herwig`); the MG5+Pythia samples are used
only as a generator-transfer test at inference.

| Network | Input features |
|---------|----------------|
| **NNkin** | 12 kinematic scalars: `dPhijj, dYjj, mjj, ptjj` + `eta, m, phi, pt` for each of the two leading jets |
| **NNj** | Leading jet: \|η\|, m, pT + up to 80 constituents × (Δη, ΔΦ, pT/pT_jet) |
| **NNZj** | NNj inputs + generator-boson 4-vector (`bosonPt, bosonEta, bosonPhi, bosonM`) |
| **NNjj** | Both jets: η (signed), m, pT + up to 80 constituents × (Δη_raw, ΔΦ, pT/pT_jet) per jet |
| **NNjjZ** | NNjj inputs + generator-boson 4-vector (`bosonPt, bosonEta, bosonPhi, bosonM`) |

NNkin is a plain MLP; the others are DeepSets over jet constituents. See
[summer26.md](summer26.md) (Step 5) for the exact architectures, training
recipe, and per-process condor submission (`train.txt`, `infer.txt`).

---

## Environment

**Must be sourced before every session** — compilation, production scripts,
Python training, inference, and plotting all rely on it:

```bash
source init.txt
```

which is equivalent to:

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_106_cuda/x86_64-el9-gcc11-opt/setup.sh
```

The entire toolchain (`root-config`, `hadd`, `python3`, `uproot`, `h5py`,
`pytorch`) is taken from this single LCG view — no mixing with system packages.

---

## Step 1 — Build `makefriends`

```bash
make
```

The Makefile resolves `root-config` and `fastjet-config` from `PATH` (set by
`init.txt`), so the binary is compiled and linked against the LCG ROOT and
FastJet libraries. If `fastjet-config` is not on `PATH` it falls back to the
CMS CVMFS installation at
`/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fastjet/3.4.1-...`.

---

## Step 2 — Run event reconstruction (`makefriends.cpp`)

`makefriends` reads a generator-level ROOT file, clusters particles with
anti-kT (R=0.4) via FastJet, computes the **pull vector** for each jet, and
writes a "friend" ROOT tree.

> **Prerequisite**: `source init.txt` must be active — the script uses the
> LCG `hadd`, `python3`, and `makefriends` (compiled against LCG ROOT) without
> any internal environment setup.

### summer26 (MadGraph5+Pythia samples, campaign-00022)

#### Quick test (100 events per sample, foreground)

```bash
./run.summer26.sh --goFast
# or a custom number of events:
./run.summer26.sh --goFast 500
```

#### Full production (20 parallel jobs × 50k events, then hadd)

```bash
./run.summer26.sh
```

**Input files** (symlinked automatically from campaign-00022 EOS) — each
process is produced with two generators, Herwig and MadGraph5+Pythia:

| File | Sample | Cross section | `--boson` |
|------|--------|---------------|-----------|
| `VBFH_herwig.root` | VBF H→inv | 2.97189 pb | H |
| `VBFH_mg5_pythia.root` | VBF H→inv | 2.88990 pb | H |
| `VBFW_herwig.root` | VBF W→inv | 7.494 pb | W |
| `VBFW_mg5_pythia.root` | VBF W→inv | 7.20333 pb | W |
| `VBFZ_herwig.root` | VBF Z→inv | 1.099 pb | Z |
| `VBFZ_mg5_pythia.root` | VBF Z→inv | 1.08581 pb | Z |
| `QCDHjj_herwig.root` | QCD H+jj (ggH) | 5.06051 pb | H |
| `QCDHjj_mg5_pythia.root` | QCD H+jj (ggH) | 4.982 pb | H |
| `QCDWjj_herwig.root` | QCD W+jj | 1733.3 pb | W |
| `QCDWjj_mg5_pythia.root` | QCD W+jj | 1646.805 pb | W |
| `QCDZjj_herwig.root` | QCD Z+jj | 340.2 pb | Z |
| `QCDZjj_mg5_pythia.root` | QCD Z+jj | 316.363 pb | Z |

**Output**: `friends/summer26/<sample>.{friend.root,h5}` for all 12 samples.

See [summer26.md](summer26.md) for sample details.

---

### Manual single-file invocation

```bash
./makefriends VBFH_mg5_pythia.root --genWeight evweight --xs 2.8899 --boson H --totEve 1000 --output VBFH_mg5_pythia.friend.root
```

`--boson H|W|Z` enables the generator-boson branches (see below); without it
they are written as −99.

### Friend ROOT tree structure

Per-jet branches are written with a slot-index suffix (`0`, `1` or `2`).
**Slots 0/1 do not correspond to leading/subleading** — they are randomised
(see below); slot 2 always holds the pT-ordered third jet (zero-padded when
fewer than 3 jets).

| Branch | Type | Description |
|--------|------|-------------|
| `jetPt{0,1,2}`, `jetEta{0,1,2}`, `jetPhi{0,1,2}`, `jetM{0,1,2}` | `Float_t` | Jet 4-vector for slots 0–2 |
| `jetNC{0,1,2}` | `Int_t` | Number of constituents |
| `jetPVM{0,1,2}` | `Float_t` | Pull vector magnitude \|t⃗\| |
| `jetPVA{0,1,2}` | `Float_t` | Pull vector angle |
| `jetSPVA{0,1,2}` | `Float_t` | Signed pull vector angle θ_s |
| `jcsPt{0,1,2}[80]`, `jcsDEta{0,1,2}[80]`, `jcsDEtaRaw{0,1,2}[80]`, `jcsDPhi{0,1,2}[80]`, `jcsM{0,1,2}[80]`, `jcsW{0,1,2}[80]` | `Float_t[80]` | Per-constituent features (up to 80, zero-padded); `jcsDEta` is sign-flipped for backward jets, `jcsDEtaRaw` is not |
| `mjj`, `dYjj`, `dPhijj`, `ptjj` | `Float_t` | Dijet event-level variables |
| `bosonM`, `bosonY`, `bosonEta`, `bosonPt`, `bosonPhi` | `Float_t` | Generator boson 4-vector (M, rapidity, η, pT, φ); −99 if `--boson` not given or no candidate found |
| `kWeight` | `Float_t` | Event weight (gen weight × xs / total events) |
| `nJets` | `Int_t` | Number of reconstructed jets stored (≤ 3) |
| `leadJetIndex` | `Int_t` | Slot index (0 or 1) that holds the true pT-leading jet |

The boson is identified at generator level depending on `--boson`: **H** =
first stable pid 25; **W** = flavour-matched (τ, ν_τ) pair with invariant mass
closest to the W pole; **Z** = (ν_e, ν̄_e) pair closest to the Z pole.

### Jet slot randomisation

The two leading jets are written into slots 0/1 with a **per-event swap**
baked in at fill time; slot 2 is exempt
([makefriends.cpp:445-475](makefriends.cpp)):

```
Line 445: bool swap = (iev % 2 == 0);             // flip for even event numbers
Line 446: o_leadJetIndex = swap ? 1 : 0;          // records which slot is truly leading
Line 475: int si = (i < 2) ? (swap ? 1 - i : i) : i;  // jets 0/1 swapped, jet 2 stays put
```

In 50 % of events the pT-leading jet lands in slot 0; in the other 50 % it
lands in slot 1.  The variable `leadJetIndex` records the true leading slot but
`dataset.py` deliberately ignores it and always reads slot 0 — so the network
sees the leading and subleading jet with equal probability and cannot learn
"leading = signal-like" as a shortcut.

---

## Step 3 — Convert to HDF5 (`root_to_hdf5.py`)

> **Note**: `run.summer26.sh` calls this automatically at the end of the
> production run — no manual step needed for the standard workflow.

To run manually on any folder:

```bash
python root_to_hdf5.py <folder>
```

Pass the folder that contains the `*.friend.root` files.  The script searches
for all matching files in that folder and writes the `.h5` output next to each
input file.

```bash
python root_to_hdf5.py friends/summer26
```

### HDF5 dataset layout

Each file contains the following datasets (gzip-compressed, level 4):

| Dataset | Shape | dtype | Description |
|---------|-------|-------|-------------|
| `jetPt` | (N, 3) | float32 | Jet pT for slots 0–2 |
| `jetEta` | (N, 3) | float32 | Jet η |
| `jetPhi` | (N, 3) | float32 | Jet φ |
| `jetM` | (N, 3) | float32 | Jet mass |
| `jetSPVA` | (N, 3) | float32 | Signed pull vector angle θ_s |
| `jetPVA` | (N, 3) | float32 | Pull vector angle |
| `jetPVM` | (N, 3) | float32 | Pull vector magnitude \|t⃗\| |
| `jetNC` | (N, 3) | int32 | Constituent count |
| `jcsPt` | (N, 3, 80) | float32 | Constituent pT (zero = absent) |
| `jcsDEta` | (N, 3, 80) | float32 | Constituent Δη (sign-flipped for backward jets) |
| `jcsDEtaRaw` | (N, 3, 80) | float32 | Constituent Δη (no sign flip) |
| `jcsDPhi` | (N, 3, 80) | float32 | Constituent ΔΦ |
| `jcsM` | (N, 3, 80) | float32 | Constituent mass |
| `jcsW` | (N, 3, 80) | float32 | Pull-vector weight |
| `kWeight` | (N,) | float32 | Event weight |
| `mjj` | (N,) | float32 | Dijet invariant mass |
| `dYjj` | (N,) | float32 | Dijet rapidity separation |
| `dPhijj` | (N,) | float32 | Dijet Δφ |
| `ptjj` | (N,) | float32 | Dijet pT |
| `bosonM`, `bosonY`, `bosonEta`, `bosonPt`, `bosonPhi` | (N,) | float32 | Generator boson M, rapidity, η, pT, φ (−99 if absent) |
| `leadJetIndex` | (N,) | int32 | Slot (0 or 1) of the true pT-leading jet |
