# fall26: friends from the Fall2026 samples (Z channel first)

The Fall2026 samples are a new, common production by the collaborator. This campaign runs them through the nov26 pipeline: inclusive friends and h5 files.

**Status: paused (2026-11-18).** The friends and h5 files are produced and validated (below). No networks are trained on them; the nov26 network code reads only the nov26 friends.

- **Scope:** the Z channel first (VBFZ and QCDZjj, Herwig and MG5+Pythia). The code already handles W and H.
- **Physics definitions:** the same as nov26 (particles, jets, gate, PF type, bc; see [nov26.md](nov26.md)). The changes are in how the inputs are read and normalised.
- **Untouched code:** none of the summer26 or nov26 friends code or outputs is changed.

Code: [makefriends_fall26.cpp](makefriends_fall26.cpp) (a fork of [makefriends_nov26.cpp](makefriends_nov26.cpp); `diff makefriends_nov26.cpp makefriends_fall26.cpp` lists every change), [fall26/manifest.py](fall26/manifest.py), [run.fall26.sh](run.fall26.sh), [root_to_hdf5_fall26.py](root_to_hdf5_fall26.py), [fall26/validate.py](fall26/validate.py).

---

## Samples

Location: `/eos/home-a/apapaefs/Projects/PullPheno/Fall2026/` (read in place; its `README.txt` describes the package).

- **Collisions and generators:** pp at 13.6 TeV. Herwig 7.3.0 (campaign B0 configuration) and MG5_aMC 3.5.15 + Pythia 8.317 (Monash, `SpaceShower:dipoleRecoil = on`).
- **Hard process:** CT14lo, α_s(M_Z) = 0.118, five flavours, fixed scales μ_R = μ_F = m_V.
- **Generation cuts, the same for both generators:** hard partons with pT > 20 GeV, |y| < 5, mjj > 20 GeV; boson mass within m_V ± 15 Γ_V.
- **Decays:** Z → ν_e ν̄_e only (no factor of three). W → τν_τ with a stable τ, W+ and W− as separate components. H stable.
- **Package:** 12 datasets of 1M events. The copy was verified (`verification.json`: `passed: true`).

| Sample | Input chunks | σ (pb) | Integration error (pb) | N_eff / N |
|---|---|---|---|---|
| `VBFZ_herwig` | 20 × 50 000 | 1.070947 | 0.0018 | 0.364 |
| `VBFZ_mg5_pythia` (chain `pythia-on`) | 100 × 10 000 | 1.072718 | 0.0011 | 1.000 |
| `QCDZjj_herwig` | 20 × 50 000 | 332.3014 | 0.26 | 0.894 |
| `QCDZjj_mg5_pythia` | 100 × 10 000 | 332.3232 | 0.40 | 1.000 |

N_eff = (Σw)²/Σw² is the effective number of events of a weighted sample. A sample with N_eff/N = 0.364 has the statistical power of 364 000 unweighted events.

**Weights.** For each event, `kWeight = evweight × pb_per_raw_weight`, in pb.
- The factor comes from `manifest.json`, per chunk. It equals σ divided by the sum of `evweight` over the whole component, so it is the same for every chunk of a component.
- Σ kWeight over a complete sample is σ.
- Pythia weights are constant.
- Herwig weights vary. They are never negative, but they are broad: in VBFZ they range from 1 to 150. The package README says they must never be replaced by unit weights.

**Compared with the nov26 inputs** (campaign-00023 Herwig and the pythia-dipole-recoil MG5+Pythia files):

| | nov26 VBFZ | Fall2026 VBFZ | nov26 QCDZjj | Fall2026 QCDZjj |
|---|---|---|---|---|
| σ Herwig (pb) | 1.099 | 1.0709 | 340.2 | 332.30 |
| σ MG5+Pythia (pb) | 1.0858 | 1.0727 | 316.36 | 332.32 |

- The two generators now agree on the Z cross sections within 0.2 % (and on all six processes within 0.4 %).
- σ_VBF/σ_QCD without selection is 1 : 310.3 (Herwig) and 1 : 309.8 (Pythia). In the nov26 inputs it was 1 : 310 and 1 : 291.

**Chunks and their reading.** Each dataset is a directory of ROOT chunks: tree `Data`, branches `numparticles`, `objects[8][10000]`, `evweight`.
- Herwig stores `objects` as **float32** and Pythia as float64.
- `makefriends.cpp` and `makefriends_nov26.cpp` bind a `Double_t` buffer. On a float32 chunk they run without any error and write garbage: jet pT of about 10²⁴ GeV, hundreds of "jets", no boson. **Never run the old binaries on Fall2026 Herwig files.**
- `makefriends_fall26` checks the stored type and reads either.

---

## What the friends contain

The nov26 format: [nov26.md](nov26.md) has the branch tables, the gate, the PF types and bc. The particle definition is the nov26 one: pT > 0.5 GeV and |η| < 4.7, excluding H, ν and τ. Both generators store every final-state particle above pT 0.1 GeV, so the definition is the same for both.

Changes against nov26:

| Branch | Type | Content |
|---|---|---|
| `kWeight` | `Float_t` | `evweight × pb_per_raw_weight` (pb per event); nov26 used σ / Σ evweight of the input file |
| `evIdx` | `Long64_t` | row of the merged sample = the chunk's offset + `chunkEntry` (nov26: the input entry) |
| `chunk` | `Int_t` | production index of the input chunk (manifest `index`), e.g. 4 for `VBFZ-000004.root` |
| `chunkEntry` | `Int_t` | entry in that chunk |
| `bosonCharge` | `Char_t` | charge of the reconstructed W (from the τ); 0 for H and Z |

**Chunk order.**
- Chunks are merged in the order given by `fall26/manifest.py`: sorted by (component, index), with the offsets of the chunks before.
- Every chunk has an even number of events, so the nov26 jet-slot randomisation (by the parity of the entry) is unchanged. The code applies it to `evIdx`.

**h5 layout.** The nov26 layout ([nov26.md](nov26.md)), plus `chunk`, `chunkEntry` (int32) and `bosonCharge` (int8), each (N,).

---

## Production

```bash
source init.source
make makefriends_fall26
./run.fall26.sh                  # PROCESSES="VBFZ QCDZjj" and NPROC=40 by default
./run.fall26.sh --goFast 2000    # quick check: first chunk of each sample -> fast_tmp_fall26/
```

- **Run:** 2026-11-17, on lxplus, about 19 minutes in all.
  - 240 chunk jobs, 40 at a time: 2.5 minutes.
  - The four merges: 8 minutes.
  - The h5 conversion: 8 minutes.
- **Outputs:** `friends/fall26/{VBFZ,QCDZjj}_{herwig,mg5_pythia}.{friend.root,h5}`.
  - friend.root files: 477–686 MB each; h5 files: 431–638 MB each.
  - Logs: `friends/fall26/run.fall26.log`, and the per-chunk logs in `run.fall26.chunks.log`.
- **Driver steps (`run.fall26.sh`):** one `makefriends_fall26` process per chunk, then a merge per sample in manifest order, then the h5 conversion. Outputs are deleted before they are rewritten (the EOS `.sys.v#` issue).

---

## Validation

**1. The float/double reader.** Pythia chunks are float64, which `makefriends_nov26` reads correctly. So the two binaries were compared on the first 1000 events of:
- a Pythia chunk;
- a float64 copy of a Herwig chunk (float32 → float64 is exact).

All 54 branches they share are identical, apart from `kWeight`, which is normalised differently. `kWeight` equals float32(`evweight` × `pb_per_raw_weight`), and `evIdx` and `chunkEntry` count 0, 1, 2, … .

**2. The full samples** (`python3 fall26/validate.py --inputs`):

| Sample | Rows | σ (pb) | Σ kWeight / σ − 1 | median bosonM (GeV) | events without a boson |
|---|---|---|---|---|---|
| `VBFZ_herwig` | 1 000 000 | 1.07095 | +1.1e-08 | 91.22 | 0 |
| `VBFZ_mg5_pythia` | 1 000 000 | 1.07272 | +4.6e-09 | 91.22 | 2 |
| `QCDZjj_herwig` | 1 000 000 | 332.301 | −3.4e-08 | 91.19 | 0 |
| `QCDZjj_mg5_pythia` | 1 000 000 | 332.323 | +1.6e-09 | 91.19 | 12 |

- **Order:** `evIdx` counts 0, 1, 2, …, and (`chunk`, `chunkEntry`) follow the manifest chunk order.
- **Weights:** `kWeight` equals `evweight × pb_per_raw_weight` read from every input chunk, event by event, for all 4M events.
- **Missing bosons:** in 14 Pythia events only one of the two Z-decay neutrinos is in the record. Pythia stores particles with pT > 0.1 GeV only, so a neutrino emitted along the beam is lost. Herwig stores every particle and has no such event. These events stay in the friends, with the boson variables at −99.

**3. Cross sections, events and N_eff at each stage.** The stages are all events, the gate, and the VBF preselection (gate, mjj > 400 GeV, |ΔY_jj| > 2.5). σ is in pb.

| Sample | σ all | N_eff all | σ gate | N gate | N_eff gate | σ VBF presel. | N VBF presel. | N_eff VBF presel. |
|---|---|---|---|---|---|---|---|---|
| `VBFZ_herwig` | 1.071 | 364 329 | 0.4044 | 342 390 | 81 523 | 0.2200 | 168 269 | 41 223 |
| `VBFZ_mg5_pythia` | 1.073 | 999 995 | 0.3985 | 371 503 | 371 501 | 0.2175 | 202 719 | 202 718 |
| `QCDZjj_herwig` | 332.3 | 894 010 | 47.78 | 148 612 | 142 518 | 7.325 | 23 116 | 22 675 |
| `QCDZjj_mg5_pythia` | 332.3 | 999 993 | 43.15 | 129 847 | 129 846 | 6.591 | 19 831 | 19 831 |

**The large Herwig VBF weights sit where the analysis looks.**
- Over all events, VBFZ Herwig keeps N_eff/N = 0.36.
- For gate events and VBF-preselected events this drops to 0.24: the 168 269 preselected events are worth 41 223 unweighted ones.

**σ_VBF / σ_QCD:**

| Generator | Campaign | All events | Gate | VBF presel. |
|---|---|---|---|---|
| Herwig | fall26 | 1 : 310 | 1 : 118 | 1 : 33 |
| Herwig | nov26 | 1 : 310 | 1 : 119 | 1 : 33 |
| MG5+Pythia | fall26 | 1 : 310 | 1 : 108 | 1 : 30 |
| MG5+Pythia | nov26 | 1 : 291 | 1 : 147 | 1 : 41 |

- **Herwig:** the ratios are the same as in nov26 at every stage.
- **MG5+Pythia:** the QCD Z+jj gate fraction changed the most.
  - nov26: 19.2 % of the QCD cross section passed the gate, against 14.4 % for Herwig.
  - fall26: 13.0 % passes, against 14.4 % for Herwig.
  - The VBF gate fractions agree between the generators in both campaigns (37–38 %).
