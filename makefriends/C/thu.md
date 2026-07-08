# thu.md — worklog Wed 2026-07-08 / Thu 2026-07-09

**Thursday update (read this first):** both pending items are done and
committed on `CNN` — Wednesday's work as `21ee021`, the per-column SR
optimization as `f5f52c9`. Results in section 5 at the bottom; sections 1–4
are the Wednesday notes, kept for context.

---

## 1. summer26.md cleanup

Removed the stale "NNj purity study — optimal cut on jet-0 score" section:
its generating code was no longer in `plots.py` (the "line 971" pointer was
dead), no figure used its thresholds, and the 3-D SR_Run3 optimization
supersedes it. HTML regenerated.

## 2. Per-channel event-level NN scores (the main work)

**Convention fixed (user decision):** every net is trained on its own
channel's **Herwig** samples and applied **only to that channel's 4 files**
(Herwig **and** MG5+Pythia — generator transfer is intended; cross-process
inference is not). The old NNj setup violated this: the Z-trained net had
scored all 12 files.

New/rewritten columns in all 12 `friends/summer26/*.h5`:

| Column | Net | Notes |
|--------|-----|-------|
| `NNj_jet0`, `NNj_jet1` | NNj{,_H,_W} | **rewritten** per-channel (H/W changed, Z identical) |
| `NNjB` | NNjB{,_H,_W} | **new**, event-level (leading jet + constituents + boson pT/η/φ) |
| `NNjjBj` | NNjjBj{,_H,_W} | **new**, event-level (jets 0+1 + jet2 + boson) |

Implementation: nine `save_scores.py` (+`.sh`/`.sub`), one per net dir,
identical up to `PROCESS = "Z"/"H"/"W"` (file list derives from it). The
NNjB/NNjjBj ones reuse the dir's own `dataset.load_features` and
`inference.run_nn`, so inference inputs are exactly the trained layout
(`_H` nets keep `bosonM`). `make_variants.py` skips `save_scores*` by design
— the copies are hand-maintained (comment updated in make_variants.py).

Ran as HTCondor GPU jobs in **3 waves** (NNj family → NNjB → NNjjBj) so no
two jobs ever write the same h5 file. Submission recipe unchanged:
`_condor_SCHEDD_HOST=bigbird24.cern.ch condor_submit -append
'requirements = (Machine =!= "b9g57n0009.cern.ch")' <dir>/save_scores.sub`.

### EOS-FUSE gotcha (cost two job retries)

On the GPU nodes, `del f[key]` on an **existing** h5 dataset intermittently
fails with `KeyError: "Couldn't delete link (bad object header version
number)"` — and persists across in-job retries (a retry loop is now in the
scripts anyway). Deletes work fine from lxplus. Remedy that worked: delete
the stale columns interactively from lxplus, resubmit — pure creates never
hit it. Also: after mutating an h5, verify from a **fresh process**
(same-process reopen can show a phantom "bad symbol table node signature").

### Verification (all green)

- All 12 files: 4 columns present, length == `len(kWeight)`, values in [0,1].
- AUC from the **stored scores** on the full MG5+Pythia files reproduces
  `auc_mg5` from each dir's `roc_data_{P}.npz` to 4 decimals, for all 9
  net/channel combos — proves the right model scored the right files:

```
NNj    Z 0.7728  H 0.7930  W 0.7602
NNjB   Z 0.8041  H 0.8009  W 0.7971
NNjjBj Z 0.8558  H 0.8837  W 0.8508
```

## 3. Downstream refresh

`summer26/plots.py` re-run (NNj_jet0-dependent figures + SR optimization).
New SR_Run3 optimum with per-channel NNj scores (Herwig, 300 fb⁻¹):

| Channel | mjj cut | \|ΔYjj\| cut | NN cut | S | B | S/(S+B) |
|---------|--------:|------------:|-------:|--:|--:|--------:|
| **H** | > 870 | > 4.7 | > 0.96 | 4 724 | 15 | **0.997** |
| **W** | > 2210 | > 4.9 | > 0.85 | 4 526 | 5 192 | **0.466** |
| **Z** | > 1320 | > 4.2 | > 0.92 | 1 275 | 3 158 | **0.288** |

H gains massively over the Z-net scores (was 0.990 at 1/3 the signal);
Z unchanged as expected. summer26.md updated (new table, new "Event-level NN
score columns" section, corrected file-map tree), summer26.html regenerated.

## 4. Natural next steps (Wed evening list — both done Thursday)

- ~~Redo the SR optimization / SR figures with `NNjB` or `NNjjBj`~~ → §5
- ~~Commit~~ → `21ee021`

---

## 5. Thursday 2026-07-09 — per-column SR optimization (commit `f5f52c9`)

`summer26/plots.py` changes:

- The four score-distribution figure blocks (normalized / SR / abs / SR-abs,
  same variants as the mjj figures) now loop over a `SCORE_COLS` list —
  `NNj_jet0`, `NNjB`, `NNjjBj` — producing
  `figs/summer26/{col}_{H,W,Z}[_Run3][_abs].pdf` (24 new PDFs; the NNj_jet0
  names are unchanged).
- The SR_Run3 optimization scans each score column **independently** (same
  grid, same S > 1000 and raw-B ≥ 10 constraints); the combined table goes to
  `figs/summer26/SR_optimization.txt`. The NNj_jet0 rows reproduce
  Wednesday's optima exactly (regression check passed).
- The `jetSPVA_SR_Run3_*` figures now use the **NNjjBj**-optimized SR — set
  by `SR_SPVA_COL = 'NNjjBj'` in plots.py (one-line switch to go back).

New optima (Herwig, 300 fb⁻¹, maximize S/(S+B)):

| Score | Channel | mjj cut | \|ΔYjj\| cut | NN cut | S | B | S/(S+B) |
|-------|---------|--------:|------------:|-------:|--:|--:|--------:|
| NNjjBj | **H** | > 1030 | > 5.0 | > 0.99 | 6 360 | 15 | **0.998** |
| NNjjBj | **W** | > 1620 | > 4.2 | > 0.99 | 15 857 | 5 192 | **0.753** |
| NNjjBj | **Z** | > 1760 | > 3.4 | > 0.99 | 2 850 | 1 019 | **0.737** |

vs NNj_jet0: W 0.466 → 0.753 at 3.5× the signal, Z 0.288 → 0.737 at 2.2×;
H already saturated but +35% signal. NNjB is intermediate (W 0.592, Z 0.524).
Full 9-row table in SR_optimization.txt and summer26.md (both updated,
html regenerated). Spot-checked `NNjjBj_W.pdf` and
`jetSPVA_SR_Run3_W_Herwig.pdf` — styling, titles and yields all consistent.

### Caveats (also in summer26.md)

- Nearly every optimum sits **at the raw-B ≥ 10 guard**: the background
  estimates rest on 10 MC events (~30% MC-stat) — visible as spiky QCD
  histograms in the SR figures.
- All three NNjjBj optima select the **last NN-cut scan bin (> 0.99)**:
  grid-limited. A finer threshold scan near 1 (e.g. 0.99–1 in 0.001 steps,
  or scan on −log(1−score)) is needed to map the true optimum — and more
  QCD MC (or a looser working point) to get off the raw-B floor.

### Natural next steps (not started)

- Finer NNjjBj threshold scan near 1 (see caveats) — the purity ceiling is
  currently set by the grid and the MC-stat guard, not by the net.
- Significance-style objective (e.g. S/√B or Asimov Z) as a cross-check —
  S/(S+B) at the raw-B floor rewards whoever reaches B ≈ 10 events first.
- Delete this worklog once read (pattern: monday/tuesday.md).
