# thu.md — state after Wed 2026-07-08 evening: per-channel NN scores in the h5 files

Session summary / notes to resume Thursday. Everything below is **done and
verified but uncommitted** on branch `CNN` (5 modified + 25 new files, incl.
this worklog) — first action Thursday: review `git status` and commit.

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

## 4. Natural next steps (not started)

- Redo the SR optimization / SR figures with `NNjB` or `NNjjBj` instead of
  `NNj_jet0` — the whole point of storing them (NNjjBj is +0.05–0.09 AUC
  over NNj on MG5). Score distributions per channel would come first.
- Commit (see top).
