# Tuesday — evening session: validate trainings, then run inference

## What was done today (2026-07-07)

Network lineup renamed and extended (see `README.md` → *Available networks* and
`summer26.md` Step 5):

| New name | Was | Change |
|----------|-----|--------|
| NNj | NNPrePro | pure rename — trained models kept |
| NNZj | NNpull (replaced) | **new**: NNj inputs + boson 4-vector (pT, η, φ, M) |
| NNjj | NNFullRaw | pure rename — trained models kept |
| NNjjZ | NNRaw (replaced) | **new**: NNjj inputs + boson 4-vector (pT, η, φ, M) |
| NNkin | — | unchanged |

Submitted at 12:29 to HTCondor (**eossubmit pool, schedd bigbird24**):

| Cluster | Job |
|---------|-----|
| 933533 | NNZj train (Z) |
| 933534 | NNZj_H train |
| 933535 | NNZj_W train |
| 933536 | NNjjZ train (Z) |
| 933537 | NNjjZ_H train |
| 933538 | NNjjZ_W train |
| 933539 | NNj save_scores — rewrites `NNj_jet0`/`NNj_jet1` into all 12 h5 files (the 2026-07-06 production regeneration wiped the old score columns; `plots.py` needs them) |

> ⚠️ These jobs live on the **eossubmit** schedds. Run
> `module load lxbatch/eossubmit` first, otherwise `condor_q` shows nothing
> and `condor_submit` of the /eos .sub files fails.

## 1. Check the training jobs finished successfully

For each of `NNZj{,_H,_W}` and `NNjjZ{,_H,_W}` (PROCESS = Z/H/W):

```bash
module load lxbatch/eossubmit
condor_q                       # should be empty once all jobs have left the queue
condor_history -limit 7        # ExitCode 0 for clusters 933533-933539
```

Then per network dir `summer26/<net>/`:

- `logs/train_gpu.<Cluster>.out` ends with `Early stopping at epoch N` (or hit
  MAX_EPOCHS=200) followed by `Best val loss: …` and
  `Saved …/<net>_<P>_loss_curve.pdf`.
- `logs/train_gpu.<Cluster>.err` has no Python `Traceback`.
- Artifacts exist in the dir: `best_model_<P>.pt`, `scaler_<P>.pkl`,
  `split_indices_<P>.npz` (P = Z for base dirs, H/W for variants).
- Loss curve PDF exists in `figs/summer26/` (`NNZj_Z_loss_curve.pdf`, …) and
  looks sane: val loss tracks train loss and converges, no divergence or
  flat-at-0.69 curve (0.693 = ln 2 means the net learned nothing).

For save_scores (933539): `summer26/NNj/logs/save_scores.<Cluster>.out` should
list all 12 files with `Written NNj_jet0, NNj_jet1 -> <file>.h5`, then spot-check:

```bash
python3 -c "import h5py; f=h5py.File('friends/summer26/VBFZ_herwig.h5'); print('NNj_jet0' in f, 'NNj_jet1' in f)"
```

## 2. Troubleshooting at the training level

- Job held / vanished: `condor_q -hold -af HoldReason`, or
  `condor_history <cluster> -af ExitCode HoldReason`.
- Python errors → `.err` log. Most likely candidates for the *new* nets:
  h5 key typos in `BOSON_FEATURES` (datasets are `bosonPt, bosonEta, bosonPhi,
  bosonM`, flat `(N,)`) or a tensor-order mismatch in the `DijetDataset` /
  `model.forward` call chain (order: jets, boson, constituents, masks).
- Out-of-memory: bump `request_memory` in `<net>/train_gpu.sub` (currently 8GB).
- Environment problems: `train_gpu.sh` sources the LCG_106_cuda view — same as
  `source init.txt`.
- After any fix, resubmit **only** the failed dir:
  `module load lxbatch/eossubmit && condor_submit summer26/<net>/train_gpu.sub`.
  (Retraining is deterministic-seeded, SEED=42.)

## 3. Then submit inference (only after all 6 trainings are validated)

```bash
module load lxbatch/eossubmit
condor_submit summer26/NNZj/infer_gpu.sub
condor_submit summer26/NNZj_H/infer_gpu.sub
condor_submit summer26/NNZj_W/infer_gpu.sub
condor_submit summer26/NNjjZ/infer_gpu.sub
condor_submit summer26/NNjjZ_H/infer_gpu.sub
condor_submit summer26/NNjjZ_W/infer_gpu.sub
```

Produces `figs/summer26/roc_nnzj_{Z,H,W}.pdf` and `roc_nnjjz_{Z,H,W}.pdf`;
`index3.php` picks them up automatically from the filenames (new "nnzj" and
"nnjjz" rows in the Inference section).

The renamed nets need **no** re-run: NNj / NNjj / NNkin models and their ROC
PDFs (`roc_nnj_*`, `roc_nnjj_*`, `roc_nnkin_*`) were carried over.

## Results (evening session, 2026-07-07)

All 7 morning jobs (933533–933539) finished cleanly: early stopping at epochs
43–86, no tracebacks, artifacts + loss curves OK, `NNj_jet0/jet1` rewritten in
all 12 h5 files. Inference submitted at 15:58 (clusters 933678–933683), done
by 16:13.

### AUC comparison (train Herwig; AUC on Herwig test split / full MG5+Py)

| Process | jets only | + boson 4-vec | jets only | + boson 4-vec |
|---------|-----------|---------------|-----------|---------------|
|         | **NNj**   | **NNZj**      | **NNjj**  | **NNjjZ**     |
| Z | 0.773 / 0.773 | 0.986 / **0.524** | 0.838 / 0.833 | 0.986 / **0.547** |
| H | 0.852 / 0.793 | 0.861 / 0.801 | 0.927 / 0.869 | 0.935 / 0.876 |
| W | 0.764 / 0.760 | 0.960 / **0.599** | 0.827 / 0.819 | 0.961 / **0.656** |

### ⚠️ Finding: `bosonM` is a generation artifact in the Herwig samples

For Z and W the Herwig AUC explodes while the MG5+Py AUC collapses below the
simple `pT lead` cut baseline. Cause: in the **Herwig QCD samples the boson is
generated exactly on-shell** while the VBF samples carry the full off-shell
spectrum; in MG5+Py both samples share the same Breit-Wigner:

| Sample | bosonM std (Herwig) | bosonM std (MG5+Py) |
|--------|--------------------:|--------------------:|
| VBFZ / QCDZjj | 20.6 / **0.0** (1 unique value!) | 5.2 / 5.2 |
| VBFW / QCDWjj | 16.3 / 3.7 | 5.4 / 5.4 |
| VBFH / QCDHjj | 0.0 / 0.0 | 0.0 / 0.0 |

`|bosonM − mZ|` **alone** gives AUC 0.999 on Herwig and 0.501 on MG5+Py
(W: 0.942 / 0.500) — reproducing the trained nets almost exactly. The Z/W nets
learned the artifact, not physics. H is immune (M ≡ 125 everywhere), which is
why NNZj_H / NNjjZ_H show consistent genuine gains (~+0.01 on both
generators) from boson pT/η.

**Next step:** drop `bosonM` from `BOSON_FEATURES` (keep pT, η, φ) and retrain
the 4 Z/W nets — or regenerate the Herwig QCD samples with off-shell bosons.

## 4. Optional follow-ups

- Re-run `python3 summer26/plots.py` once `NNj_jet0` is back in the h5 files
  (score-distribution + SR-optimization sections were computed from the old
  production).
- Compare AUCs: NNZj vs NNj and NNjjZ vs NNjj quantify what the boson
  4-vector adds, per process and per generator.
- `git status` is left uncommitted on branch `CNN` — the renames/deletions are
  staged, new dirs (`NNZj*`, `NNjjZ*`) are untracked; commit when happy.
