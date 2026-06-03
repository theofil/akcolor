# testAUC — jet-feature dominance experiment

## What this folder is

Both `NNpull` (DeepSets) and `NNtrans` (Transformer) use identical input features and reach
the same AUC, which is suspicious given how different the two architectures are.

**Hypothesis:** the 6 jet-level scalars passed directly to the final MLP classifier —
especially `jetPVM` (pull-vector magnitude) and `jetSPVA` (signed pull-vector angle), which
are hand-computed aggregates over the constituent-level pull information — already contain
most of the discriminating power. The constituent-level encoder (DeepSets sum-pool vs
Transformer attention) is largely irrelevant because the model gets the answer handed to it
as a jet-level scalar.

**Experiment:** train the same two architectures **without** jet-level features (constituent
features only), then compare all four AUCs on one plot.

```
testAUC/
  deepsets/          ← DeepSetsNoJet   (identical phi/rho to NNpull, x_jet never seen)
  transformer/       ← TransNoJet      (identical transformer to NNtrans, x_jet never seen)
  deepsets_nofun/    ← DeepSetsNoFun   (4 constituent features: jcsDEta, jcsDPhi, jcsM, jcsPt — jcsW removed)
  deepsets_noscale/  ← DeepSetsNoScale (3 constituent features: jcsDEta, jcsDPhi, jcsPt/jetPt — fraction only)
  compare_roc.py                    → roc_comparison.pdf
  roc_comparison_DeepSetNoFun.py    → roc_comparison_DeepSetNoFun.pdf
  roc_comparison_DeepSetNoScale.py  → roc_comparison_DeepSetNoScale.pdf
  compare_gpu.sh / compare_gpu.sub
```

**This folder is temporary** — delete it once `roc_comparison.pdf` has been inspected.

---

## How to run

> Must be done from a terminal with a valid Kerberos ticket.

### Step 1 — source the environment

```bash
source /eos/home-t/theofil/work/akcolor/makefriends/C/init.source
```

This loads the LCG_106_cuda software stack **and** `module load lxbatch/eossubmit`,
which is required for HTCondor to accept `/eos` paths in submit files.

### Step 2 — submit the two training jobs

```bash
condor_submit /eos/home-t/theofil/work/akcolor/makefriends/C/testAUC/deepsets/train_gpu.sub
condor_submit /eos/home-t/theofil/work/akcolor/makefriends/C/testAUC/transformer/train_gpu.sub
```

Each job runs `train.py` followed by `inference.py` on a GPU worker (`longlunch` flavour).
Outputs written to the respective sub-folder:
- `best_model.pt`, `scaler.pkl`, `split_indices.npz`, `loss_curve.pdf`
- `scores.npz`, `roc_deepsets_nojet.pdf` / `roc_trans_nojet.pdf`

### Step 3 — monitor

```bash
condor_q
```

### Step 4 — submit the comparison job (after both training jobs finish)

```bash
condor_submit /eos/home-t/theofil/work/akcolor/makefriends/C/testAUC/compare_gpu.sub
```

This runs `compare_roc.py`, which loads all four models (NNpull, NNtrans, DeepSetsNoJet,
TransNoJet) on the same test set and writes `testAUC/roc_comparison.pdf`.

---

## Interpreting the result

| Outcome | Interpretation |
|---|---|
| nojet AUC << full AUC (both models) | Jet features dominate — hypothesis confirmed |
| nojet AUC ≈ full AUC | Constituent features alone are sufficient |
| TransNoJet AUC > DeepSetsNoJet AUC | Attention extracts something extra that sum-pooling misses |

---

## Results

All four models evaluated on the same held-out test set (27 077 events: 17 405 signal, 9 672 background):

| Model | Architecture | Constituent features | Epochs | Best val loss | AUC |
|-------|-------------|----------------------|--------|---------------|-----|
| NNpull | DeepSets | jet scalars (6) + jcsDEta, jcsDPhi, jcsM, jcsPt, jcsW (5×80) | — | — | **0.8132** |
| NNtrans | Transformer | jet scalars (6) + jcsDEta, jcsDPhi, jcsM, jcsPt, jcsW (5×80) | — | — | **0.8130** |
| DeepSets (no jet) | DeepSets | jcsDEta, jcsDPhi, jcsM, jcsPt, jcsW (5×80) | 45 | 0.3801 | 0.8041 |
| Trans (no jet) | Transformer | jcsDEta, jcsDPhi, jcsM, jcsPt, jcsW (5×80) | — | — | 0.8036 |
| DeepSets (no fun) | DeepSets | jcsDEta, jcsDPhi, jcsM, jcsPt (4×80) — jcsW removed | 67 | 0.3806 | 0.8036 |
| DeepSets (no scale) | DeepSets | jcsDEta, jcsDPhi, jcsPt/jetPt (3×80) — jcsM, jcsW removed | 71 | 0.4018 | 0.7810 |

**Jet scalars** (6): `|jetEta|`, `jetM`, `jetNC`, `jetPVM`, `jetPt`, `jetSPVA`

**Constituent features** (per constituent, up to 80 per jet): `jcsDEta`, `jcsDPhi`, `jcsM`, `jcsPt`, `jcsW`

**Interpretation of the constituent-ablation series (DeepSets only, no jet scalars):**

- **Removing jet scalars** (no jet, 5 features → AUC 0.8041 vs 0.8132): costs ~0.9% AUC.
  Most discriminating power already resides in the raw constituents; jet-level aggregates add
  only a marginal gain. The two architectures (DeepSets vs Transformer) are interchangeable.

- **Removing jcsW** (no fun, 4 features → AUC 0.8036): essentially zero effect (−0.05%).
  The pull-vector weight per constituent carries no information beyond what `jcsDEta`, `jcsDPhi`,
  `jcsM`, `jcsPt` already provide.

- **Removing jcsM and replacing jcsPt with jcsPt/jetPt** (no scale, 3 features → AUC 0.7810):
  costs ~2.3% AUC relative to no-jet, a clear degradation. The absolute momentum scale (`jcsPt`)
  and the constituent mass (`jcsM`) together contribute meaningfully — making `jcsPt` dimensionless
  (fraction of jet pT) and dropping `jcsM` removes real discriminating information.

Output files: `roc_comparison.pdf`, `roc_comparison_DeepSetNoFun.pdf`, `roc_comparison_DeepSetNoScale.pdf`

---

## DeepSets in LHC analyses

DeepSets (Zaheer et al. 2017, NeurIPS) is a general framework for learning on unordered sets:
a shared per-element MLP (phi) maps each element to a latent vector, a symmetric aggregation
(sum/mean/max) produces a permutation-invariant representation, and a second MLP (rho) maps
that to the output. In HEP the "elements" are jet constituents, tracks, or event particles —
variable-size collections with no natural ordering — making DeepSets a natural fit.
NNpull follows this paradigm exactly (phi → masked sum → rho + jet scalars).

| Experiment | Analysis / Algorithm | Physics task | Reference |
|---|---|---|---|
| Theory | Deep Sets | General set learning framework | Zaheer et al., NeurIPS 2017 |
| Phenomenology | Energy Flow Networks (EFN/PFN) | Quark/gluon jet tagging (IRC-safe) | arXiv:1810.05165 |
| ATLAS | Impact-parameter flavour tagging | b/c-jet ID from tracks | ATL-PHYS-PUB-2020-014 |
| ATLAS | DIPz pile-up jet rejection | Jet origin regression along beamline | arXiv:2512.10819 |
| CMS | HL-LHC Level-1 Trigger jet tagging | Multi-class jet ID for di-Higgs triggers | arXiv:2509.24371 |
| CMS | Pile-up mitigation | Event-wide PU suppression | arXiv:2503.02860 |
| LHCb | Inclusive flavour tagging | B⁰/B⁰s production flavour from full-event tracks | arXiv:2602.15625 |
| LHCb | Fast inclusive flavour tagging | Real-time triggered analysis variant | arXiv:2404.14145 |
| Phenomenology | Deep Set Autoencoders | Model-agnostic anomaly detection / BSM searches | arXiv:2109.01695 |
