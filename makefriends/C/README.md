# Analysis workflow — VBF H→inv vs QCD dijet

Classification of VBF H→invisible signal against QCD dijet background using
jet pull vectors and constituent information.

For analysis steps (plots, NN training, inference) see
[summer25.md](summer25.md) and [summer26.md](summer26.md).

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

### summer25 (Herwig 13.6 samples)

#### Quick test (100 events per sample, foreground)

```bash
./run.summer25.sh --goFast
# or a custom number of events:
./run.summer25.sh --goFast 500
```

#### Full production (20 parallel jobs × 50k events, then hadd)

```bash
./run.summer25.sh
```

**Input files** (symlinked automatically from apapaefs EOS, Herwig 13.6):

| Symlink | Sample | Cross section |
|---------|--------|---------------|
| `VBFHtoInv.root` | VBF H→inv | 3.901 pb |
| `QCDHtoInv.root` | ggH+jj (QCD background) | 2.26114 pb |

**Output**: `friends/summer25/VBFHtoInv.{friend.root,h5}`, `friends/summer25/QCDHtoInv.{friend.root,h5}`

---

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

**Input files** (symlinked automatically from campaign-00022 EOS):

| File | Sample | Cross section |
|------|--------|---------------|
| `VBFH_mg5_pythia.root` | VBF H→inv | 2.889 pb |
| `QCDHjj_mg5_pythia.root` | QCD H+jj (ggH) | 4.967 pb |
| `VBFZ_mg5_pythia.root` *(commented out)* | VBF Z→inv | 1.084 pb |
| `VBFW_mg5_pythia.root` *(commented out)* | VBF W→inv | 7.186 pb |
| `QCDZjj_mg5_pythia.root` *(commented out)* | QCD Z+jj | 315.3 pb |
| `QCDWjj_mg5_pythia.root` *(commented out)* | QCD W+jj | 1642 pb |

**Output**: `friends/summer26/VBFH_mg5_pythia.{friend.root,h5}`, `friends/summer26/QCDHjj_mg5_pythia.{friend.root,h5}`

See [summer26.md](summer26.md) for details on enabling Z/W samples.

---

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

---

## Step 3 — Convert to HDF5 (`root_to_hdf5.py`)

> **Note**: `run.summer25.sh` and `run.summer26.sh` call this
> automatically at the end of the production run — no manual step needed for
> the standard workflow.

To run manually on any folder:

```bash
python root_to_hdf5.py <folder>
```

Pass the folder that contains the `*.friend.root` files.  The script searches
for all matching files in that folder and writes the `.h5` output next to each
input file.

```bash
# summer25
python root_to_hdf5.py friends/summer25

# summer26
python root_to_hdf5.py friends/summer26
```

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
