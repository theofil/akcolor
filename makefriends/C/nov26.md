# nov26 campaign: inclusive Z-channel friends with PF-type constituents

The nov26 friends cover the neutral (Z) channel first. They differ from summer26 in five ways:

- **No preselection.** Every input event is written. The loose dijet gate is stored as a flag (`gate`), not applied as a cut.
- **Two jets only.** At most the two leading jets are stored, with their constituents. There is no 3rd jet.
- **CMS-PF-like particle cuts.** Particles need pT > 0.5 GeV and |η| < 4.7.
- **Particle type per constituent.** Each constituent carries a CMS-PF-style type and its charge.
- **Beam connection per jet.** Four new scalars (`bc_in`, `bc_out`, `bc_met`, `bc_dphi`) measure the activity between each jet and its own beam.

Code: [makefriends_nov26.cpp](makefriends_nov26.cpp) (a fork of [makefriends.cpp](makefriends.cpp); `diff makefriends.cpp makefriends_nov26.cpp` lists every change), [run.nov26.sh](run.nov26.sh), [root_to_hdf5_nov26.py](root_to_hdf5_nov26.py). Plots: [nov26/bc_display.py](nov26/bc_display.py) and [nov26/bc_plots.py](nov26/bc_plots.py). Jet networks NNjPF and ParTj, with five-fold cross-fitting: [nov26/training.md](nov26/training.md).
None of the summer26 code or outputs are touched.

**Status (2026-11-17):** the code, including bc, is validated, and the full production is done: `friends/nov26/{VBFZ,QCDZjj}_{herwig,mg5_pythia}.{friend.root,h5}` (see Production). The jet networks NNjPF and ParTj are trained with five-fold cross-fitting, and their out-of-fold scores are in the h5 files ([nov26/training.md](nov26/training.md)).

---

## Samples

These are the same inputs and cross sections as the summer26 Z channel. All four use `--boson Z` (Z → ν_e ν̄_e) and have 1M events each.

| Sample | Input | σ (pb) |
|---|---|---|
| `VBFZ_herwig` | `campaign-00023/merged/VBFZ_herwig.root` | 1.099 |
| `VBFZ_mg5_pythia` | `pythia-dipole-recoil/merged/VBFZ/VBFZ_mg5_pythia_dipolerecoil.root` | 1.08581 |
| `QCDZjj_herwig` | `campaign-00023/merged/QCDZjj_herwig.root` | 340.2 |
| `QCDZjj_mg5_pythia` | `pythia-dipole-recoil/merged/QCDZjj/QCDZjj_mg5_pythia_dipolerecoil.root` | 316.363 |

`campaign-00023` lives under `/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/`. `pythia-dipole-recoil` lives under `/eos/user/a/apapaefs/Projects/Pull/`.

---

## Definitions

**Particles used for clustering**
- Final-state particles from the input `Data` tree.
- Excluded: the Higgs (|pid| = 25), neutrinos (|pid| = 12, 14, 16) and τ (|pid| = 15).
- Kept only if pT > 0.5 GeV (`JC_PT_MIN`) and |η| < 4.7 (`JC_ETA_MAX`).
- The η cut has practically no effect on jets with |η| < 3.

**Jets**
- anti-kT with R = 0.4. A jet is kept if pT > 30 GeV (`JET_PT_MIN`) and |η| < 3.0 (`JET_ETA_MAX`).
- The two leading jets are stored. For each, up to 80 constituents (`NC_MAX`) are stored, pT-ordered and zero-padded.

**Dijet gate** (`gate`, a boolean). An event passes when all three hold:
1. at least 2 jets,
2. mjj > 0 (`MJJ_MIN = 0`, so no real mass cut),
3. η_j1 · η_j2 < 0, i.e. the two leading jets are in opposite hemispheres.

`gate` is only a flag: events that fail it are still written.

**Weight.** `kWeight = evweight × σ / ΣW`, where ΣW is the sum of `evweight` over the whole input file.
- `kWeight` is a cross section per event, in pb.
- Σ `kWeight` over a complete sample equals σ.
- The expected yield is Σ `kWeight` × L, with L in pb⁻¹.

**Slot randomisation** (same trick as summer26). In 2-jet events, the leading jet is written to slot 1 for even input entries and to slot 0 for odd ones. `leadJetIndex` records which slot holds the leading jet. A single jet always sits in slot 0 (`leadJetIndex = 0`). Events with no jets have `leadJetIndex = −1`.

**Particle type** (`jcsPFType`). The codes follow CMS `reco::PFCandidate::ParticleType`:

| Code | PF type | Generator particles |
|---|---|---|
| 1 | charged hadron (h) | π±, K±, p, p̄ |
| 2 | electron (e) | e± |
| 3 | muon (μ) | μ± |
| 4 | photon (γ) | γ |
| 5 | neutral hadron (h0) | K0L, n, n̄ |
| 0 | — | empty constituent slot |

- The rule: |pid| 11 → e, 13 → μ, 22 → γ; any other charged particle → h; any other neutral particle → h0.
- The charge (`jcsCharge`, −1/0/+1) comes from row 5 of the input `objects` array. Row 4 holds the PDG id.
- Herwig and Pythia produce the same stable species. There are no K0S, Λ or Σ in the record.

**Beam connection** (`jetBCin`, `jetBCout`, `jetBCmet`, `jetBCdphi`). These measure the activity between a jet and its own beam; they are computed by `beamConnection()` in `makefriends_nov26.cpp`.
- y_j and φ_j are the jet's rapidity and azimuth, and s = sign(y_j), with s = +1 for y_j = 0.
- The particles are the clustering particles defined above.
- The region R holds the particles with s·(y_i − y_j) > 0.4 and |y_i| < 4.7. It starts 0.4 in rapidity beyond the jet axis and runs to the acceptance edge on the jet's own side: toward +y for a forward jet, toward −y for a negative-rapidity jet.

| Branch | Definition |
|---|---|
| `jetBCin` | Σ pT_i over R with \|Δφ(i, j)\| < 0.4, a strip of half-width 0.4 around the jet's φ |
| `jetBCout` | Σ pT_i over R with \|Δφ\| ≥ 0.4, so bc_in + bc_out = Σ pT over R |
| `jetBCmet` | pT of the vector sum of the particle pT over R |
| `jetBCdphi` | Phi_mpi_pi(φ_j − φ of that vector sum), in [−π, π]; −99 if R is empty |

- Every Δφ uses ROOT's `TVector2::Phi_mpi_pi`, so the strip wraps correctly across φ = 0/2π.
- The original spec wrote the region as y_i − sign(y_j)·y_j > 0.4. That equals y_i > \|y_j\| + 0.4, which for a negative-rapidity jet points at the opposite beam. The mirrored form above was chosen on 2026-11-17.
- Example events: `figs/nov26/bc_regions.pdf`, made by `nov26/bc_display.py`.

---

## Friend tree `events`

There is one entry per input event, in input order: entry i is input entry i, and `evIdx` stores i. This holds because `run.nov26.sh` merges the jobs in numeric order. Branches with the suffix 0/1 are slots, not pT order.

| Branch | Type | Description |
|---|---|---|
| `evIdx` | `Long64_t` | input `Data` entry number |
| `nJets` | `Int_t` | jets stored (≤ 2) |
| `nJetsAll` | `Int_t` | all jets with pT > 30 GeV, \|η\| < 3 |
| `gate` | `Bool_t` | dijet gate (see Definitions) |
| `jetPt{0,1}`, `jetEta{0,1}`, `jetPhi{0,1}`, `jetM{0,1}` | `Float_t` | jet 4-vector (rounded to 0.01 GeV / 0.001) |
| `jetSPVA{0,1}`, `jetPVA{0,1}`, `jetPVM{0,1}` | `Float_t` | pull vector: signed angle, angle, magnitude (−99 if undefined) |
| `jetNC{0,1}` | `Int_t` | number of constituents (all of them, not capped at 80) |
| `jetBCin{0,1}`, `jetBCout{0,1}`, `jetBCmet{0,1}`, `jetBCdphi{0,1}` | `Float_t` | beam connection (see Definitions). With an empty region, in/out/met are 0 and dphi is −99 |
| `jcsPt{0,1}[80]`, `jcsDEta{0,1}[80]`, `jcsDEtaRaw{0,1}[80]`, `jcsDPhi{0,1}[80]`, `jcsM{0,1}[80]`, `jcsW{0,1}[80]` | `Float_t[80]` | as in summer26 (`jcsDEta` is sign-flipped for backward jets, `jcsDEtaRaw` is not; `jcsW` = pull weight) |
| `jcsPFType{0,1}[80]` | `Char_t[80]` | PF type code (table above) |
| `jcsCharge{0,1}[80]` | `Char_t[80]` | −1/0/+1 |
| `mjj`, `dYjj`, `dPhijj`, `ptjj` | `Float_t` | from the two leading jets; −99 if `nJetsAll` < 2 |
| `bosonM`, `bosonY`, `bosonEta`, `bosonPt`, `bosonPhi` | `Float_t` | generator Z (ν_e ν̄_e pair closest to the Z pole); −99 if not found |
| `kWeight` | `Float_t` | pb per event (see Definitions) |
| `leadJetIndex` | `Int_t` | slot of the pT-leading jet (−1 without jets) |

**Empty slots** are all zero, including the PF type. A slot is valid when its index is below `nJets` (equivalently `jetPt > 0`). A constituent is valid when `jcsPt > 0` (equivalently `jcsPFType > 0`).

**Histograms in the file**
- `hSumW`, `hNJets` and `hJetPt1` are filled for all events.
- `jetShape`, `jetShapeRaw`, `hLeadJetSPVA` and `hLeadJetSPVM` are filled for gate events. In summer26 they used the mjj > 400 / |ΔY| > 2.5 selection instead.

## HDF5 layout

[root_to_hdf5_nov26.py](root_to_hdf5_nov26.py) writes `<sample>.h5` next to each friend file, compressed with gzip level 4. The dataset names are the same as in summer26, but there are 2 slots instead of 3.

| Shape | dtype | Datasets |
|---|---|---|
| (N, 2) | float32 | `jetPt, jetEta, jetPhi, jetM, jetSPVA, jetPVA, jetPVM, jetBCin, jetBCout, jetBCmet, jetBCdphi` |
| (N, 2) | int32 | `jetNC` |
| (N, 2, 80) | float32 | `jcsPt, jcsDEta, jcsDEtaRaw, jcsDPhi, jcsM, jcsW` |
| (N, 2, 80) | int8 | `jcsPFType, jcsCharge` |
| (N,) | float32 | `kWeight, mjj, dYjj, dPhijj, ptjj, bosonM, bosonY, bosonEta, bosonPt, bosonPhi` |
| (N,) | int32 / int64 / bool | `nJets, nJetsAll, leadJetIndex` / `evIdx` / `gate` |
| (N, 2) | float32 | `NNjPF`, `ParTj`: out-of-fold jet-tagger score of the jet in each slot, −1 for an empty slot (written by `nov26/save_scores.py`; see [nov26/training.md](nov26/training.md)) |
| (N,) | int8 | `subfold`: the event's cross-fitting sub-fold (0–9) |

**For downstream code:** these files hold every event, so a loader has to select `gate` (and any mjj/ΔY cut) itself. The summer26 loaders divide `jcsPt` by `jetPt`, which gives NaN for an empty slot. Rerunning `root_to_hdf5_nov26.py` drops the `NNjPF`, `ParTj` and `subfold` datasets; rerun `nov26/save_scores.py` afterwards (no retraining needed).

---

## How to run

```bash
source init.source
make makefriends_nov26           # a bare `make` still builds only the summer26 binary
./run.nov26.sh --goFast 2000     # first 2000 events per sample -> fast_tmp_nov26/ (friend.root + h5)
./run.nov26.sh                   # full production: 20 jobs x 50k events per sample -> friends/nov26/
python3 nov26/bc_display.py      # two example events with the bc regions -> figs/nov26/bc_regions.{pdf,png}
python3 nov26/bc_plots.py        # bc distributions of gate events -> figs/nov26/bc_distributions.{pdf,png}
```

- The full production runs locally on lxplus, not on condor. It ends by converting the four friend files to h5.
- Single file: `./makefriends_nov26 <input.root> --genWeight evweight --xs <σ> --boson Z [--totEve N --skip N] --output <out>.friend.root`.
- Manual conversion: `python3 root_to_hdf5_nov26.py <file.friend.root ...>`, or pass a folder to convert every `*.friend.root` in it.
- Outputs are deleted before they are rewritten. On EOS, overwriting a file in place leaves hidden `.sys.v#` version directories behind.

---

## Validation (2026-11-17)

**1. Input charge.** On the first 2000 events of all four inputs (2.85M particles), `objects[5]` equals the PDG charge for every species. There were zero mismatches and no unexpected species.

**2. Regression against summer26.** A scratch build of `makefriends_nov26.cpp` with the summer26 particle cuts (pT > 0, |η| < 20) was run next to the summer26 `makefriends` on the first 3000 events of `VBFZ_herwig` and `QCDZjj_mg5_pythia`.
- The nov26 events with `gate && mjj > 400 && |dYjj| > 2.5` (621 and 92 events) equal the summer26 output exactly, in order.
- The comparison covered slots 0/1 (jet scalars, pull vector, all six jcs arrays), the dijet variables, `kWeight`, the boson variables and `leadJetIndex`.
- Also equal: summer26 `nJets` = min(`nJetsAll`, 3), and the inclusive histograms.
- In total, 43 of 43 comparisons agree for each sample, both before and after bc was added.

**3. Smoke test** (`./run.nov26.sh --goFast 2000`). Every check passed on all four samples:
- **Event bookkeeping**
  - `evIdx` runs 0..1999 and `nJets` = min(`nJetsAll`, 2).
  - Empty slots are all zero.
  - The `leadJetIndex` rule holds, the leading slot has the larger pT, and both slot parities occur.
  - `gate` matches its definition, and the dijet variables are −99 exactly when `nJetsAll` < 2.
- **Constituents**
  - `jcsPt > 0` exactly when `jcsPFType > 0`, and every stored constituent has pT > 0.5 GeV.
  - The number of filled slots equals min(`jetNC`, 80), and constituents are pT-ordered.
  - Charge matches type: h, e and μ have |q| = 1; γ and h0 have q = 0.
- **Beam connection**
  - bc_in, bc_out and bc_met are ≥ 0, and bc_met ≤ bc_in + bc_out.
  - bc_dphi is in [−π, π], and is −99 exactly when the region is empty (0.3–0.6% of jets).
- **h5:** every dataset equals its ROOT branch (shapes and dtypes as listed).

| Sample (first 2000 events) | 0 / 1 / 2 jets stored | ≥3 jets | gate |
|---|---|---|---|
| `VBFZ_herwig` | 145 / 682 / 1173 | 196 | 790 |
| `VBFZ_mg5_pythia` | 161 / 616 / 1223 | 295 | 765 |
| `QCDZjj_herwig` | 441 / 854 / 705 | 177 | 290 |
| `QCDZjj_mg5_pythia` | 378 / 707 / 915 | 313 | 371 |

PF composition per stored jet: the pT-weighted mean fraction, and the share of constituents.

| Sample | h | e | μ | γ | h0 | share h / γ / h0 |
|---|---|---|---|---|---|---|
| `VBFZ_herwig` | 61.6 % | 0.3 % | 0.3 % | 27.4 % | 10.3 % | 53.3 / 40.3 / 5.8 % |
| `VBFZ_mg5_pythia` | 61.4 % | 0.4 % | 0.2 % | 28.0 % | 10.0 % | 52.8 / 40.3 / 6.3 % |
| `QCDZjj_herwig` | 61.9 % | 0.4 % | 0.3 % | 26.0 % | 11.4 % | 54.8 / 37.9 / 6.7 % |
| `QCDZjj_mg5_pythia` | 62.4 % | 0.2 % | 0.2 % | 26.1 % | 11.1 % | 54.2 / 38.2 / 7.2 % |

**4. Independent bc cross-check.** `nov26/bc_display.py` re-clusters the events with the same FastJet library (through PyROOT) and recomputes bc with numpy, separately from the C++ code.
- It was run on every stored jet of the first 2000 events of all four samples: 10 891 jets.
- There were 0 mismatches. The worst relative deviation is 6×10⁻⁸, which is float32 storage rounding.
- The re-clustered jets also reproduce the stored pT, η and φ.

**5. Example events.** `figs/nov26/bc_regions.pdf` shows the bc regions for two VBFZ_herwig events.
- A forward jet (event 16, y_j = +2.21): R runs from y = 2.61 to 4.7.
- A negative-rapidity jet (event 67, y_j = −2.42, φ_j = 6.10): R runs from y = −2.82 to −4.7, and its |Δφ| < 0.4 strip wraps across φ = 0/2π.
- In both panels the Python and C++ values agree.

**Size.** With bc included, the tree takes 483–700 bytes per event, and the h5 about the same.

---

## Production (2026-11-17)

`./run.nov26.sh` ran locally on lxplus with 20 jobs per sample. The four friend files took 19 minutes and the h5 conversion another 10. The log is `friends/nov26/run.nov26.log`.

| Sample | ≥ 3 jets | gate events (share of events / of σ) | gate && mjj > 400 && \|dYjj\| > 2.5 vs summer26 | friend.root / h5 |
|---|---|---|---|---|
| `VBFZ_herwig` | 10.3 % | 375 453 (37.5 % / 37.6 %) | 203 331 vs 205 250 (−0.9 %) | 0.65 / 0.61 GB |
| `VBFZ_mg5_pythia` | 14.9 % | 379 898 (38.0 % / 38.0 %) | 205 548 vs 207 827 (−1.1 %) | 0.69 / 0.65 GB |
| `QCDZjj_herwig` | 7.9 % | 144 200 (14.4 % / 14.4 %) | 21 962 vs 22 821 (−3.8 %) | 0.47 / 0.43 GB |
| `QCDZjj_mg5_pythia` | 14.9 % | 192 208 (19.2 % / 19.2 %) | 29 266 vs 30 452 (−3.9 %) | 0.56 / 0.51 GB |

- Each file has 1 000 000 entries, one per input event, and every smoke-test check also passes on the full files.
- **Weight closure.** Σ `kWeight` equals σ to better than 10⁻⁸ in every file, as it must when every event is kept.
- **Comparison with summer26.** With summer26's cuts applied, 0.9–1.1 % fewer VBF events and 3.8–3.9 % fewer QCD events pass than in summer26. This is the expected direction: dropping particles below 0.5 GeV lowers jet pT slightly, so fewer jets pass 30 GeV and fewer pairs pass mjj > 400.
- **Empty bc region:** 0.3–0.6 % of stored jets.

**bc distributions: a first look** (`figs/nov26/bc_distributions.pdf`, made by `nov26/bc_plots.py`). These are gate events, both stored jets, kWeight-weighted and normalised to unit area.
- bc_out is softer for VBF Z than for QCD Z+jj. VBF peaks near 5 GeV and QCD near 15 GeV, with a longer QCD tail.
- |bc_dphi| peaks at 0 for VBF Z, so the beam-side activity tends to line up with the jet. For QCD it is nearly flat.
- bc_in and bc_met differ less between the processes.
- Herwig and Pythia agree closely for VBF Z. For QCD Z+jj, Pythia has more beam-side activity than Herwig.
