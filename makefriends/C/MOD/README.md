# MOD friends: Zenodo CMS Open Data → friends/*.h5

Friends-building pipeline for the public [Zenodo 3340205](https://zenodo.org/records/3340205)
family ("CMS 2011A Open Data | Jet Primary Dataset | pT > 375 GeV | MOD HDF5 Format",
Komiske/Mastandrea/Metodiev/Naik/Thaler, MIT — arXiv:1908.08542), reproducing the
same physics computations and output schema as `../makefriends.cpp` +
`../root_to_hdf5.py` so `summer26/NN*/dataset.py` can read `MOD/friends/*.h5` with
no code changes.

This is **not** `makefriends.cpp` pointed at a new file — the input is structurally
different (HDF5, pre-formed jets, no full event particle record), so
`make_mod_friends.py` is a standalone Python/h5py converter instead. See
"Showstoppers" below for why, and the module docstring in `make_mod_friends.py`
for the exact per-branch derivation.

## Data sources

| Sample | Zenodo record | Content |
|---|---|---|
| `CMS2011A_data` | [3340205](https://zenodo.org/records/3340205) | Real 2011A collision data, Jet300 HLT, pT>375 GeV |
| `QCD_170-300` … `QCD_1800-inf` | 3341500, 3341498, 3341419, 3364139, 3341413, 3341502, 3341770, 3341772 | Pythia6 QCD dijet MC, binned in hard-parton p̂T |

Each SIM record contains two file kinds per index: `GEN*_compressed.h5` (truth-only,
no PF candidates — **not used**) and `SIM*_compressed.h5` (detector-level, has
`pfcs`/`quality`/`npv` — **this is what's converted**). `download_zenodo.py` only
fetches the `SIM*`/`CMS_Jet300*` files.

## Showstoppers / missing info found

1. **No H/W/Z boson truth anywhere in this dataset family.** It's a generic QCD
   jet-substructure sample. `bosonM/Y/Eta/Pt/Phi` are always `-99` (same sentinel
   `makefriends.cpp` uses when `--boson` is omitted). SIM carries `hard_pid`
   (quark/gluon flavor of the hard parton) but that's unused here.
2. **Only the two hardest jets' own PF candidates are recorded per event** —
   nothing else in the event. So there's no re-clustering step (unlike
   `makefriends.cpp`, which clusters from the full particle list); the two
   pre-formed jets are taken as-is. `NJETMAX` is 2, never 3 — there is no "slot 2"
   in this output at all (per instruction: only save what's actually there).
3. **Most events only have 1 jet stored, not 2.** Empirically (first DATA file):
   64% of events have exactly 1 jet row, 36% have 2. Only the 2-jet events can
   form a dijet pair and are kept — matching the spirit of `makefriends.cpp`'s own
   "≥2 jets" requirement. For low-p̂T SIM bins the 2-jet fraction is far smaller
   (~1% for 170-300 GeV) since a jet rarely fluctuates up to the 375 GeV threshold
   at all in that bin, let alone twice.
4. **No per-event MC weight in DATA** in the conceptual sense (it's real collision
   data) — but empirically `jets_f` *does* carry a `weight` column even for DATA
   (constant across a file, ~4.4e-7), contrary to what the Zenodo docs implied.
   `kWeight` = this `weight` column directly, for both DATA and SIM — no
   `--xs`/`sumWtot` rescaling needed.
5. **Dropped the VBF-style opposite-eta-hemisphere cut** that `makefriends.cpp`
   applies (`jets[0].eta()*jets[1].eta()<0`) — confirmed with same-sign eta pairs
   being common in this generic QCD sample; that cut was specific to the VBF
   forward-jet topology and has no motivation here (confirmed with user).
   Kept selection: the dataset's own recommended `|jet_eta|<1.9` and
   `quality>=2` ("medium"), both jets of a pair must pass.
6. **eta vs rapidity**: PF candidates only carry rapidity `y`, not pseudorapidity.
   `jetPVA/PVM/SPVA` (`fillPV` in `makefriends.cpp`) are rapidity-based already, so
   no conversion is needed there. `jcsDEta`/`jcsDEtaRaw` are pseudorapidity-based
   in the original code; the converter reconstructs each constituent's eta from
   its own `(pt, y, m)` (`pz = sqrt(m²+pt²)·sinh(y)`, `eta = asinh(pz/pt)`) to stay
   consistent with that convention.

## Files

- `download_zenodo.py [sample ...]` — idempotent downloader (Zenodo API file list,
  skips files already present with matching size, skips `GEN*` files).
- `make_mod_friends.py [sample ...]` — converter; reads `raw/<subdir>/*.h5`,
  writes `friends/<sample>.h5`. Streams one input file at a time into resizable
  HDF5 datasets (`append_h5`) rather than buffering a whole sample's events in
  Python lists — the big pT-hat bins (e.g. 800-1000 GeV: 3.9M events) OOM-killed
  the first version of this script against this container's ~35GB memory cgroup
  limit before that fix.
- `plot_observables.py [sample]` — mjj, jet pT, jet eta, theta_21 (pull angle),
  and leading-jet signed pull-vector angle, from a friends file. Writes to `figs/`.

## Output schema

Identical dataset names/dtypes to `root_to_hdf5.py`, minus the jet-slot-2 axis:
`jetPt/Eta/Phi/M/SPVA/PVA/PVM/NC` shape `[N,2]`; `jcsPt/DEta/DEtaRaw/DPhi/M/W`
shape `[N,2,80]`; `kWeight/mjj/dYjj/dPhijj/ptjj/bosonM/Y/Eta/Pt/Phi/leadJetIndex`
shape `[N]`. Verified byte-for-byte-compatible key set/dtypes against an existing
`../friends/summer26/*.h5` file (differing only in the 2-vs-3 jet axis).
