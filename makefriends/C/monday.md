# monday.md — jet2 + boson 4-vector extension of summer26 friends/h5

**Goal**: extend the summer26 friend trees and h5 files with
(a) the **3rd pt-ordered jet ("jet2")** — never randomized, always softer than jet0/jet1 —
including its constituent arrays, and
(b) the generator-level **W/Z/h boson 4-vector** (`M, Y, PT, Phi`) as scalar event variables.

Written 2026-07-05 (planning session, no code changed). Execute the steps below.

---

## How the bosons appear in the samples (verified 2026-07-05)

Facts established by directly inspecting the campaign ROOT files (2000–3000 events per sample,
scripts kept in the planning-session scratchpad; trivial to re-derive with uproot).

The `Data` tree stores **final-state particles only** in `objects[8][10000]`:
columns 0=E, 1=px, 2=py, 3=pz, 4=pid, 5=charge, 6–7=unused(0). No status codes,
no intermediate resonances.

| Channel | What's in the record | Correct 4-vector recipe |
|---------|---------------------|------------------------|
| **h** (VBFH, QCDHjj) | Higgs is **kept stable** (models H→inv): exactly one pid==25 per event, m = 125.00 exactly | Take the pid==25 particle directly |
| **W** (VBFW, QCDWjj) | **No pid ±24.** W→τν_τ with the **τ kept stable** (m=1.78; τ and ν's are already excluded from jet clustering in `makefriends.cpp`) | p4(τ) + p4(ν_τ), flavor-matched: pid 15 ↔ −16, pid −15 ↔ 16 |
| **Z** (VBFZ, QCDZjj) | **No pid 23.** Z→ν_e ν̄_e, **always electron flavor** (pid ±12) | p4(ν_e) + p4(ν̄_e), i.e. a (12, −12) pair |

Ambiguity resolution (extra ν's from hadron decays in ~13% of Z events; a 2nd shower τ in
~1.2% of MG5 W events / ~0.03% Herwig, never 0 τ): among flavor-consistent candidates pick the
pair with invariant mass **closest to the pole mass** (m_W = 80.379, m_Z = 91.1876 GeV).
Verified peaks: W median 80.4 GeV; Z median 91.2 GeV (Breit-Wigner in VBF, essentially on-shell
in QCDZjj); H = 125 exactly. The mass-closest choice differs from naive hardest-pairing in only
~1.5% (Z) / ≤0.1% (W) of events.

**Do not use** the `theETmiss[4]` branch — its components do not match the ν system.
(`theHiggses`/`numHiggses` confirm 1 Higgs/event but reading pid 25 from `objects` suffices.)

---

## Step 1 — `makefriends.cpp`

- `NJETMAX` 2 → 3. The existing `memset` zero-padding then covers jet2 automatically
  (jet2 slots stay 0 when only 2 jets pass cuts — same convention the datasets already use
  as validity mask `jcsPt > 0`).
- New branches, same style/rounding as existing (pt/m 2dp, eta/phi 3dp):
  `jetPt2, jetEta2, jetPhi2, jetM2, jetSPVA2, jetPVA2, jetPVM2, jetNC2` and
  `jcsPt2, jcsDEta2, jcsDEtaRaw2, jcsDPhi2, jcsM2, jcsW2` (each `[NC_MAX]`).
- Fill loop (currently line ~404): `nmax = min((int)jets.size(), NJETMAX)`; index map
  `si = (i < 2) ? (swap ? 1-i : i) : i` — **randomization stays confined to jets 0/1**;
  jet2 is always `jets[2]` (pt-sorted 3rd jet, hence always softer than both leading jets by
  construction). `leadJetIndex` semantics unchanged. `nJets` becomes 2 or 3 (was always 2).
- Jet2 passes the same jet cuts (pt > 30, |η| < 3) automatically — it comes from the same
  `jets` vector. The dijet event selection itself is unchanged (two leading jets only).
- New CLI arg `--boson H|W|Z` (explicit per sample, **no auto-detect** — avoids mis-tagging
  the rare Z events containing a shower τ pair as W). Absent/empty → boson branches = −99.
- Boson extraction inside the existing particle loop (reuse E/px/py/pz/pid already read):
  - **H**: the pid==25 particle.
  - **W**: collect τ (|pid|=15) and ν_τ (|pid|=16); pick the flavor-consistent
    (15,−16)/(−15,16) pair with m closest to 80.379.
  - **Z**: collect pid 12 and −12; pick the pair with m closest to 91.1876.
  - New branches `bosonM, bosonY, bosonEta, bosonPt, bosonPhi` via TLorentzVector
    (`.M(), .Rapidity(), .Eta(), .Pt(), .Phi()`; M/PT 2dp, Y/Eta/Phi 3dp); −99 sentinel
    if no candidate (essentially never happens).

## Step 2 — `root_to_hdf5.py`

- Stack three columns instead of two for the per-jet and per-constituent variables:
  `np.stack([arrays[f'{var}0'], arrays[f'{var}1'], arrays[f'{var}2']], axis=1)`
  → shapes (N,3) and (N,3,80). Rename `SCALAR_PAIRS` → `SCALAR_JETS` for clarity.
- Add `bosonM, bosonY, bosonPt, bosonPhi` to `SCALAR_EVT`.
- **Backward compatible** (verified): all `summer26/NN*/dataset.py` index `[:, 0]` / `[:, 1]`
  only, so existing trained models keep working; `plots.py` reads friend branches by name —
  unaffected.

## Step 3 — `run.summer26.sh`

- Add `--boson H/W/Z` to all 24 `makefriends` invocations (12 goFast lines + the 12 entries of
  the `samples=` array, which gain a 6th `:boson` field parsed in `run_sample`).
  Boson letter follows the sample name: VBFH/QCDHjj → H, VBFW/QCDWjj → W, VBFZ/QCDZjj → Z.

## Step 4 — Documentation

- Update `summer26.md`: new branches table, the boson-reconstruction recipe above, and a note
  that friends/h5 must be regenerated.
- `makefriends.py` (legacy Python version) stays untouched.

## Step 5 — Verification

1. `make`, then `./run.summer26.sh --goFast 2000` → outputs in `fast_tmp/`.
2. Checks on the fast_tmp friend/h5 files:
   - `bosonM` peaks: H = 125 exactly; W median ≈ 80.4; Z median ≈ 91.2; count of −99 ≈ 0.
   - `jetPt2 <= min(jetPt0, jetPt1)` wherever `jetPt2 > 0`; jet2 slots all-zero when `nJets==2`.
   - jet0/jet1 still randomized (`leadJetIndex` ~50/50); jet2 never swapped.
   - h5 shapes (N,3) / (N,3,80); existing `summer26/NNFullRaw/dataset.py` loads the new h5
     unmodified.
3. Only after the goFast checks pass: full production `./run.summer26.sh`
   (rerun of all 12 samples, ~178 GB inputs) + h5 conversion.
