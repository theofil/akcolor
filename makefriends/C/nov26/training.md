# nov26 jet networks: training sample, five-fold cross-fitting and scores

Jet taggers are trained on the nov26 friends to tell a jet from VBF Z production (signal) from a jet in QCD Z+jj production (background). All of them are built from the DeepSets net NNj of summer26. They differ in which beam-connection (bc) values of the jet ([nov26.md](../nov26.md), Definitions) they get on top of NNj's inputs, and two of them, Kj and KjBC, have no constituent branch:

| Net | bc inputs | Purpose |
|---|---|---|
| **NNj** | none | baseline: NNj trained on the nov26 folds |
| **NNjBC** | bc_in, bc_out, bc_met, bc_dphi | NNj + the beam connection |
| NNjBCin, NNjBCout, NNjBCmet, NNjBCdphi | that one value | add-one: the gain of one value on its own, NNjBC*v* − NNj |
| NNjBC_noin, NNjBC_noout, NNjBC_nomet, NNjBC_nodphi | the other three | drop-one: what only that value adds, NNjBC − NNjBC_no*v* |
| **Kj** | none, and no constituents | kinematics-only reference: NNj without the constituent branch, on the jet's \|η\|, m, pT alone |
| **KjBC** | bc_in, bc_out, bc_met, bc_dphi, and no constituents | Kj + the beam connection |

Each net sees one jet at a time. The jet's inputs are |η|, m and pT, plus its bc values; each constituent contributes the signed Δη (`jcsDEta`), ΔΦ and pT/pT_jet (not for Kj and KjBC). The nets are trained on Herwig only; MG5+Pythia is used purely to test how well they transfer to another generator.

Kj, KjBC, NNj and NNjBC form a 2 × 2 of constituents (no / yes) × beam connection (no / yes) on the same jets, folds and recipe. Comparing them at equal background efficiency measures what the constituents and the beam connection add beyond the jet kinematics ([roc.summary.md](../roc.summary.md), section 5).

The bc values are correlated (bc_in + bc_out is the total pT in the region, and bc_met ≤ that total), so a value's add-one gain and its drop-one loss bracket its contribution: the first counts what it shares with the others, the second only what none of the others carries.

This note explains how the training sample is built, why it is weighted, and how five-fold cross-fitting gives every event a score that no network has seen. The code is in this folder: `common.py` (sample, folds, weights, bc inputs, the list of nets), `model.py`, `train_fold.py`, `train_gpu.sub`, `save_scores.py`, `roc_scores.py`, `bootstrap_scores.py` and `gain_curves.py`.

---

## 1. One training example is one jet

The nov26 files keep every event and at most its two leading jets (slots 0 and 1). Each stored jet has one of three **roles**:

| Role | Which jets | Used for |
|---|---|---|
| **A** | gate events, slot 0 | training and validation |
| **B** | gate events, slot 1 | never trained on; only scored |
| **C** | events failing the gate: every stored jet (1 or 2 per event) | training and validation |

Slot 0 of a gate event is the leading jet in half of the events and the subleading one in the other half. This is the same random slot assignment NNj used, so the net cannot learn "the leading jet looks like signal".

The second jet of a gate event (role B) is held back from training altogether, for two reasons. Training on both jets of an event would put two strongly correlated examples (same event, same colour flow) into the training set. Keeping one jet out also gives a jet in every gate event that is guaranteed to play no part in fitting any network.

## 2. How much training data this gives

Counts of Herwig jets over all 1M events of each sample:

| | VBF Z (signal) | QCD Z+jj (background) |
|---|---|---|
| Role A: gate events, slot 0 | 375 453 | 144 200 |
| Role C: events with 1 stored jet | 343 854 | 428 326 |
| Role C: events with 2 stored jets in the same hemisphere (2 jets each) | 2 × 210 410 | 2 × 197 230 |
| **Training population A + C** | **1 140 127** | **966 986** |
| Gain over role A alone | ×3.0 | ×6.7 |
| Gain over the VBF-preselected sample NNj used (gate, mjj > 400, \|ΔY\| > 2.5) | ×5.6 | ×44 |

**The gain is largest where it matters most.** In summer26 the background was the limiting sample: NNj trained on about 16 thousand QCD jets. Each nov26 fold model trains on about 677 thousand QCD jets, roughly 42 times more.

**Weighted events.** Herwig events carry generator weights. A weighted sample of N events is only as informative as N_eff unweighted events, where N_eff = (Σw)² / Σw². If a few events had huge weights, N_eff would be much smaller than N. Here N_eff / N is 97–100 % (1.10M of 1.14M for VBF, 0.96M of 0.97M for QCD), so the weights cost almost no statistics.

Only two jets per event are stored. The 3rd and later jets of gate == 0 events would add another 4–7 %, which is not worth storing a third jet for.

## 3. Population shift: gate jets and gate == 0 jets

Roles A and C come from two different **populations** of jets, and they differ in two ways.

**1. What the jets look like.** Gate == 0 jets are softer, and their |η| differs:

| Mean value | VBF, gate (A) | VBF, gate == 0 (C) | QCD, gate (A) | QCD, gate == 0 (C) |
|---|---|---|---|---|
| pT [GeV] | 135 | 110 | 80 | 70 |
| \|η\| | 1.54 | 1.37 | 1.17 | 1.32 |
| constituents | 15.1 | 14.2 | 14.3 | 13.4 |

**2. The mix of the two classes.** In the training set of fold model 0, 72 % of the gate jets are VBF jets, but only 48 % of the gate == 0 jets are.

**Why the second point is a problem.** A network learns whatever lowers its loss. Suppose it can tell from a jet's kinematics alone which population the jet probably came from. Then "looks like a gate jet" becomes a hint for "VBF", because gate jets are mostly VBF in this mixture. That hint has nothing to do with VBF-versus-QCD physics; it is an accident of how the training sample was assembled. A network that used it would score jets partly by population, and its performance on gate jets could be distorted.

**The fix: balance the classes inside each population.** Every training jet gets the weight

  w = N_P / (2 · N_{P,c}),

where N_P is the number of training jets in its population P (gate or gate == 0) and N_{P,c} is the number of those that belong to its class c (VBF or QCD).

| Fold model 0, training jets | VBF share before | VBF share after weighting | Total weight of the population |
|---|---|---|---|
| Gate jets (A) | 72.3 % | 50.0 % | 363 755 = the number of jets |
| Gate == 0 jets (C) | 48.2 % | 50.0 % | 1 111 218 = the number of jets |

After weighting, both populations are half VBF and half QCD. Knowing which population a jet came from then tells the network nothing about its class, so the shortcut no longer helps. Each population keeps its full weight, so the larger gate == 0 sample still contributes all of its statistics. This is the per-population version of what NNj already did with `pos_weight`, which balanced the classes over the whole sample. The validation loss uses the same weights, computed on the validation jets.

**What the weighting does not change.** The jets of the two populations still look different, and the networks learn from both. Performance therefore matters most on gate jets, and that is reported separately (roles A + B, section 10). If the gate == 0 jets ever turned out to hurt gate-jet performance, the next steps would be to fine-tune on gate jets only, or to give the network a population flag.

## 4. What a fold is, and what cross-fitting does

A network's score on an event it was trained on is too optimistic: to some degree the network has memorised that event. **Honest** performance, and honest scores for an analysis, need events the network never saw. summer26 used a single 70 / 15 / 15 split, so only 15 % of the events had honest scores. For the Z channel that left about 34 raw background events at the tightest working point.

A **fold** is one of K equal, non-overlapping groups of events. **Cross-fitting** trains K networks: network k never sees fold k and is the one that scores it. Every event then gets exactly one honest score.

nov26 uses **10 sub-folds and 5 networks**. Network k:
- **infers** (scores) sub-folds 2k and 2k+1: 20 % of the events;
- **validates** on sub-fold 2k+2 (mod 10): 10 %;
- **trains** on the other 7: 70 %.

```
                 sub-fold:  0  1  2  3  4  5  6  7  8  9
network 0 (fold 0)          I  I  V  T  T  T  T  T  T  T
network 1 (fold 1)          T  T  I  I  V  T  T  T  T  T
network 2 (fold 2)          T  T  T  T  I  I  V  T  T  T
network 3 (fold 3)          T  T  T  T  T  T  I  I  V  T
network 4 (fold 4)          V  T  T  T  T  T  T  T  I  I
                 T = training (70 %), V = validation (10 %), I = inference (20 %)
```

Each column has exactly one I, so every event is scored once.

- **Training** sub-folds are the data the network learns from.
- **Validation** monitors over-fitting. It decides when to stop training (early stopping) and when to lower the learning rate. Because it influences the network, its scores are not honest either, so it is kept separate from inference.
- **Inference** sub-folds are scored by the finished network and by nothing else. These are the scores written into the h5 files.

**Folds are made of events, not jets.** The two jets of an event share their underlying event and colour flow. If one jet were in training and the other in inference, information would leak from training into the "honest" score. Keeping whole events together prevents this. It also means a gate event's held-back jet (role B) is scored by a network that saw neither of the event's jets.

**Why 70 / 10 / 20.** A 10 % validation set still holds about 97 thousand QCD and 114 thousand VBF jets. That measures the validation loss to about 0.3 %, which is plenty for early stopping. Taking 20 % instead would cost each network 17 % of its training data for no benefit.

## 5. What "stratified" means

If events were dropped into folds purely at random, the folds would differ slightly by chance. One fold might, for example, get a few hundred more gate events than another.

A **stratified** split forces every fold to have the same composition as the whole sample, for the properties that matter. Here those are the sample (which also fixes the class and the generator), whether the event passed the gate, and how many jets it stored (0, 1 or 2).

It works like dealing cards:
1. Sort the events into groups by (sample, gate, number of stored jets).
2. Shuffle each group with a fixed seed.
3. Deal its events round-robin into the 10 sub-folds.

Every sub-fold then receives a tenth of every group. On the full files the sub-folds of a group differ by at most one event.

The assignment depends only on the file content and the seed, so all jobs and all nets get the same folds. It is stored in the h5 files as `subfold`.

## 6. Five networks, or one network trained on more events?

Each of the five networks trains on 70 % of the events. Would one network trained on, say, 85 % be better? These are two different questions:

1. How good is one network? More training data helps, with diminishing returns.
2. How well can we measure and use the score on these same events? That needs honest scores, and the more events have them, the better.

| Scheme | Training data per network | Events with an honest score | Training cost |
|---|---|---|---|
| Single 70 / 15 / 15 split (summer26) | 70 % | 15 % | 1 network |
| One network on 85 %, 15 % validation, no test | 85 % | none | 1 network |
| **Five-fold cross-fitting (nov26)** | 70 % | **100 %** | 5 networks |

**What cross-fitting buys.** Each network gets the same share of training data as summer26's single split, but honest scores now cover every event, 6.7 times more than a 15 % test split. With the VBF preselection applied, all 21 962 QCD events get honest scores, instead of about 3 300 (15 %). The statistical error at a tight working point shrinks by about √6.7 ≈ 2.6.

**What a single bigger network would buy.** It would have 21 % more training jets. With about 680 thousand QCD jets per network, 42 times what NNj had, that is expected to be a small effect. It can be checked with a learning curve: train one fold on 25, 50 and 100 % of its training set and see whether the performance is still rising.

**The five networks are not five different taggers.** They have the same architecture, the same recipe, and 70 % training sets that overlap heavily: any two share 4 or 5 of their 7 training sub-folds. Their scores are statistically equivalent. The spread of their individual AUCs (section 10) even measures how much a network varies from one training to the next.

**Scoring new events.** For events none of them trained on (a new MC sample, or data), average the five scores. This ensemble is usually slightly better than any single network. For the nov26 MC events, use the out-of-fold score stored in the h5.

**MG5+Pythia.** No network trains on Pythia, so any of them would give an honest score. Pythia events still get sub-folds by the same rule and are scored by their fold's network. Herwig and Pythia scores then come from the same mixture of networks, so any difference between them is a generator effect and not a network effect.

## 7. Training recipe

All nets share the recipe of NNj:

| | |
|---|---|
| Architecture | DeepSets (`model.py`): φ 3 → 64 → 64 on each constituent, masked sum; ρ (64 + jet inputs) → 128 → 64 → 1 with BatchNorm and dropout 0.3. Kj and KjBC: no φ, the same ρ on the jet inputs alone |
| Jet inputs | \|η\|, m, pT, then the net's bc inputs: 3 for NNj and Kj, 8 for NNjBC and KjBC |
| Parameters | 21 825 (NNj) to 22 465 (NNjBC); 9 217 (Kj) and 9 857 (KjBC) |
| Optimiser | Adam, learning rate 1e-3, batch size 512 |
| Learning-rate schedule | ×0.5 after 5 epochs without improvement |
| Early stopping | 15 epochs without improvement, at most 100 epochs |

**The bc inputs** (`common.load_jets`, column `x_bc`):
- bc_in, bc_out and bc_met enter as log(1 + x / GeV). Their tails reach several hundred GeV.
- bc_dphi enters as cos and sin, two inputs that count as one value in the ablation. The sign is kept: `jcsDPhi` is not mirrored, so the net can relate the side of the beam-side activity to the constituents.
- A jet whose bc region is empty (bc_dphi = −99, 0.3–0.6 % of jets) gets cos = sin = 0; its bc_in, bc_out and bc_met are 0.

Common to all nets:
- Each fold fits its own StandardScaler to its training jets. It covers all jet inputs (bc included) and the constituent features of valid constituents.
- The loss is binary cross-entropy weighted as in section 3.
- The constituent axis is trimmed to the largest constituent count in each batch.
- The nets run as condor GPU jobs, one per net and fold: `source init.source; mkdir -p nov26/logs; condor_submit -name bigbird24.cern.ch nov26/train_gpu.sub` submits all 60 (12 nets × 5 folds). Kj and KjBC were added later, from a copy of the file whose queue line lists only them (`queue 5 net in (Kj KjBC)`); resubmitting the full file queues every net again.
- Kj and KjBC keep only the first constituent slot in memory, and their model ignores it; everything else, including the scaler fit, is unchanged.

Each job writes `<net>/best_model_fold{k}.pt`, `scaler_fold{k}.pkl`, `history_fold{k}.npz` and `scores_fold{k}.npz`, plus a loss curve `figs/nov26/<net>_fold{k}_loss.pdf`; the condor logs go to `nov26/logs/`.

## 8. The scores in the h5 files

After all five jobs of a net finish, `python3 nov26/save_scores.py --net <net>` (or `--net all`) writes into each of the four `friends/nov26/*.h5`:

| Dataset | Shape | Content |
|---|---|---|
| one per net: `NNj`, `NNjBC`, `NNjBC<v>`, `NNjBC_no<v>`, `Kj`, `KjBC` | (N, 2) float32 | out-of-fold score (sigmoid, 0–1) of the jet in each slot; −1 for an empty slot |
| `subfold` | (N,) int8 | the event's sub-fold; network k scored sub-folds 2k and 2k+1 |

Before writing, `save_scores.py` checks that every stored jet was scored exactly once, by the right network.

For example, to get the scores of both jets of the VBF-preselected gate events:

```python
import h5py, numpy as np
with h5py.File('friends/nov26/VBFZ_herwig.h5') as f:
    sel = f['gate'][:] & (f['mjj'][:] > 400) & (np.abs(f['dYjj'][:]) > 2.5)
    s = f['NNjBC'][:][sel]          # (n, 2): both jets of every selected event, all out-of-fold
```

Rerunning `root_to_hdf5_nov26.py` rewrites the h5 files and drops these datasets. Afterwards, rerun `python3 nov26/save_scores.py --net all`; no retraining is needed.

## 9. Inference for the 1st and the 2nd jet of a gate event

A gate event has two stored jets. The **1st jet** is in slot 0 (role A) and the **2nd jet** in slot 1 (role B). Slot 0 holds the pT-leading jet in half of the events and the subleading jet in the other half; `leadJetIndex` records which.

**Both jets are scored by the same network, and that network saw neither of them.** Folds are made of events, so the two jets share their event's sub-fold. The one network whose inference sub-folds include that sub-fold scores both jets, using the same weights and the same fold scaler. The event was in neither its training nor its validation sub-folds.

**Following one event through the five networks.** Take a gate event in sub-fold 2:

| Network | Sub-fold 2 is its… | 1st jet (slot 0, role A) | 2nd jet (slot 1, role B) |
|---|---|---|---|
| 0 | validation | used to decide when to stop training | not used |
| **1** | **inference** | **scored: this is the score stored in the h5** | **scored: this is the score stored in the h5** |
| 2 | training | trained on | not used |
| 3 | training | trained on | not used |
| 4 | training | trained on | not used |

- Only network 1's scores are written to the h5 file.
- No network uses the 2nd jet for anything other than being scored, and only network 1 scores it.
- Validation sub-folds are always even (2k + 2). An event in an odd sub-fold, for example 3, is therefore never used for validation, and its 1st jet is trained on by four networks instead of three.

**For the stored scores, the two jets are on an equal footing.** Both come from a network that never saw the event. Keeping the 2nd jet out of training still adds three things:

1. **One jet per gate event in training.** The training set never holds two correlated jets of the same event (same underlying event, same colour flow).
2. **A built-in check.** The slot assignment is random, so 1st and 2nd jets come from the same distribution. In VBF Z Herwig both slots have a mean pT of 134.7 GeV and a mean |η| of 1.54. Their score distributions and AUCs must therefore agree within statistics, and a difference would point to leakage or to a slot bias. `roc_scores.py` reports role A and role B separately for this comparison.
3. **A jet that no network ever trained on.** If a single network is later trained on all role A and C jets, for example to use on data, the 2nd jets of the gate events are still jets it never saw. Each of them shares an event with a 1st jet that network did train on, though, so for the nov26 MC events the cross-fitted scores remain the cleaner choice.

**Events that fail the gate (role C)** follow the same rule. Their one or two stored jets share the event's sub-fold and are scored by its inference network. The only difference is that a network training on that sub-fold uses all of the event's jets.

**Reading the scores of the leading and subleading jet of gate events:**

```python
import h5py, numpy as np
with h5py.File('friends/nov26/VBFZ_herwig.h5') as f:
    gate, lead, s = f['gate'][:], f['leadJetIndex'][:], f['NNjBC'][:]
rows = np.flatnonzero(gate)              # gate events always store two jets
s_lead = s[rows, lead[rows]]             # score of the pT-leading jet
s_sub = s[rows, 1 - lead[rows]]          # score of the subleading jet
```

Both arrays are out-of-fold scores. They can be combined into an event-level score: a product, a mean, or a small second-level model. A second-level model trained on these scores has to reuse the same sub-folds (the `subfold` column) so that it stays honest too.

## 10. Results (2026-11-18)

All 50 jobs (condor cluster 1648339) finished cleanly, and `save_scores.py --net all` wrote the ten score columns into the four h5 files, after checking that every stored jet was scored exactly once by the right model. The values in this section are unweighted, as before in this note. The event-weighted values, the comparison with summer26 NNj on the same events and the full ablation tables, with uncertainties, are in [roc.summary.md](../roc.summary.md).

### Training

| Net | Stopped at epoch | Best validation loss | Wall time per fold |
|---|---|---|---|
| NNj | 47–85 | 0.6088–0.6113 | 7–23 min |
| NNjBC | 61–88 | 0.5849–0.5875 | 10–14 min |
| NNjBCin / out / met / dphi | 50–99 | 0.6073–0.6102 / 0.5993–0.6023 / 0.5940–0.5968 / 0.6053–0.6079 | 6–19 min |
| NNjBC_noin / noout / nomet / nodphi | 53–99 | 0.5853–0.5881 / 0.5887–0.5915 / 0.5935–0.5968 / 0.5866–0.5895 | 6–19 min |
| Kj | 26–51 | 0.6143–0.6171 | 7–11 min |
| KjBC | 54–88 | 0.5910–0.5939 | 9–12 min |

- The jobs ran on H100 MIG slices and V100s. Each fold model trained on 1.47M jets and was validated on 211k, as in section 2.
- The validation losses already rank the nets: on every fold NNjBC has the lowest loss. Without bc_in the loss is unchanged, and it rises most without bc_met.
- Kj and KjBC were trained later, in condor cluster 1649485 (10 jobs, all exit 0), with the same folds and recipe. KjBC's validation loss (0.5910–0.5939) lies below NNj's (0.6088–0.6113) on every fold. The loss includes the gate == 0 jets, where the beam connection helps most, but KjBC also beats NNj on gate jets (AUC 0.7666 vs 0.7624).

### NNj and NNjBC

| Net | Jets | AUC Herwig | AUC MG5+Pythia | ε_S(1 %) Herwig | ε_S(1 %) MG5+Pythia |
|---|---|---|---|---|---|
| NNj | gate jets (roles A + B) | 0.7624 | 0.7631 | 0.0949 | 0.1017 |
| NNj | 1st jet of gate events (role A) | 0.7622 | 0.7629 | 0.0940 | 0.1003 |
| NNj | 2nd jet of gate events (role B) | 0.7627 | 0.7634 | 0.0959 | 0.1031 |
| NNj | gate == 0 jets (role C) | 0.7175 | 0.7078 | 0.0608 | 0.0586 |
| NNj | VBF-preselected gate jets | 0.6832 | 0.6946 | 0.0464 | 0.0516 |
| NNjBC | gate jets (roles A + B) | 0.7754 | 0.7764 | 0.1022 | 0.1098 |
| NNjBC | 1st jet of gate events (role A) | 0.7752 | 0.7760 | 0.1024 | 0.1095 |
| NNjBC | 2nd jet of gate events (role B) | 0.7755 | 0.7768 | 0.1021 | 0.1102 |
| NNjBC | gate == 0 jets (role C) | 0.7497 | 0.7366 | 0.0949 | 0.0910 |
| NNjBC | VBF-preselected gate jets | 0.6980 | 0.7081 | 0.0579 | 0.0590 |

- **1st vs 2nd jet:** the AUCs agree within 0.0008 ± 0.0008 for both nets and both generators, and ε_S(1 %) within 1.2σ (`bootstrap_scores.py`). There is no sign of leakage or of a slot bias.
- **Spread of the five fold models** on their own Herwig gate jets: NNj AUC 0.7603–0.7648 (rms 0.0017), NNjBC 0.7733–0.7777 (rms 0.0016). For ε_S(1 %), the rms is 0.0022 for NNj and 0.0023 for NNjBC.

### What the beam connection adds

NNjBC − NNj per fold model, on the same jets (mean ± standard error of the five; every fold model gains):

| Jets | AUC, Herwig | AUC, MG5+Pythia | ε_S(1 %), Herwig | ε_S(1 %), MG5+Pythia |
|---|---|---|---|---|
| Gate jets | +0.0129 ± 0.0003 | +0.0133 ± 0.0002 | +0.0073 ± 0.0010 | +0.0087 ± 0.0014 |
| Gate == 0 jets | +0.0322 ± 0.0001 | +0.0287 ± 0.0001 | +0.0341 ± 0.0011 | +0.0325 ± 0.0004 |
| Gate jets of VBF-preselected events | +0.0148 ± 0.0008 | +0.0136 ± 0.0002 | +0.0117 ± 0.0015 | +0.0086 ± 0.0022 |

Background rejection at fixed signal efficiency (`bootstrap_scores.py`), with NNjBC / NNj as the ratio:

| Signal efficiency | 5 % | 10 % | 20 % | 30 % | 50 % | 70 % | 90 % |
|---|---|---|---|---|---|---|---|
| NNj rejection, gate jets, Herwig | 249 | 92.6 | 32.0 | 16.3 | 6.6 | 3.2 | 1.6 |
| NNjBC / NNj, gate jets, Herwig | 1.08 | 1.11 | 1.13 | 1.13 | 1.10 | 1.08 | 1.04 |
| NNjBC / NNj, gate jets, MG5+Pythia | 1.16 | 1.12 | 1.14 | 1.13 | 1.10 | 1.08 | 1.04 |
| NNj rejection, VBF-preselected, Herwig | 93.3 | 37.6 | 14.7 | 8.5 | 4.0 | 2.3 | 1.4 |
| NNjBC / NNj, VBF-preselected, Herwig | 1.32 | 1.18 | 1.18 | 1.13 | 1.08 | 1.05 | 1.03 |
| NNjBC / NNj, VBF-preselected, MG5+Pythia | 1.17 | 1.17 | 1.15 | 1.10 | 1.08 | 1.05 | 1.03 |

- **The gain transfers to the other generator.** The AUC gain on gate jets is the same on MG5+Pythia as on Herwig (+0.013). For scale: the particle type and charge of the discontinued NNjPF were worth about 0.002 (Herwig) and 0.0006 (MG5+Pythia) on the same gate jets (NNjPF 0.7643 and 0.7637).
- **It is largest where NNj is weakest.** For the gate == 0 jets, the AUC rises by 0.03 and ε_S(1 %) by 0.034, which is more than half again. For gate jets, NNjBC rejects 4–16 % more background at fixed signal efficiency (10–14 % between 10 % and 50 % signal efficiency), and for VBF-preselected jets up to 32 % more at 5 % signal efficiency on Herwig.
- **Which bc values carry the gain:** add-one and drop-one results per value are in [roc.summary.md](../roc.summary.md), section 4. In short, bc_met and bc_dphi carry most of it, bc_out adds a Herwig-only gain on its own, and bc_in adds nothing.

ROC curves: `figs/nov26/roc_nets.pdf` (unweighted, NNj and NNjBC, three jet categories with the rejection ratio), `figs/roc_summary.pdf` (event-weighted, with summer26 NNj), `figs/nov26/bc_ablation.pdf` (the ablation) and `figs/nov26/gain_curves.pdf` (Kj, KjBC, NNj, NNjBC at equal background efficiency).

### Kj and KjBC: the kinematics-only references (2026-11-18)

Unweighted, as above. The event-weighted values, the AUC at fixed jet kinematics and the comparison at equal background efficiency are in [roc.summary.md](../roc.summary.md), section 5.

| Net | Jets | AUC Herwig | AUC MG5+Pythia | ε_S(1 %) Herwig | ε_S(1 %) MG5+Pythia |
|---|---|---|---|---|---|
| Kj | gate jets (roles A + B) | 0.7528 | 0.7556 | 0.0839 | 0.0929 |
| Kj | gate == 0 jets (role C) | 0.7105 | 0.7009 | 0.0573 | 0.0552 |
| Kj | VBF-preselected gate jets | 0.6630 | 0.6773 | 0.0397 | 0.0452 |
| KjBC | gate jets (roles A + B) | 0.7666 | 0.7697 | 0.0935 | 0.1043 |
| KjBC | gate == 0 jets (role C) | 0.7429 | 0.7303 | 0.0872 | 0.0855 |
| KjBC | VBF-preselected gate jets | 0.6810 | 0.6944 | 0.0498 | 0.0516 |

- **1st vs 2nd jet** (Herwig AUC, roles A / B): Kj 0.7527 / 0.7529, KjBC 0.7667 / 0.7665. No slot bias.
- **Spread of the five fold models** on their own Herwig gate jets: Kj 0.7514–0.7544 (rms 0.0012), KjBC 0.7641–0.7694 (rms 0.0017).

Paired AUC differences per fold model, on the same jets (mean ± standard error of the five; every fold model gains):

| Difference | Measures | Gate jets, Herwig | Gate jets, MG5+Pythia | Gate == 0, Herwig | Gate == 0, MG5+Pythia | VBF presel., Herwig | VBF presel., MG5+Pythia |
|---|---|---|---|---|---|---|---|
| NNj − Kj | the constituents | +0.0096 ± 0.0003 | +0.0075 ± 0.0005 | +0.0070 ± 0.0002 | +0.0070 ± 0.0001 | +0.0201 ± 0.0003 | +0.0172 ± 0.0004 |
| KjBC − Kj | the beam connection | +0.0138 ± 0.0004 | +0.0141 ± 0.0004 | +0.0324 ± 0.0003 | +0.0295 ± 0.0001 | +0.0180 ± 0.0009 | +0.0172 ± 0.0004 |
| NNjBC − KjBC | the constituents, given the beam connection | +0.0087 ± 0.0002 | +0.0067 ± 0.0001 | +0.0068 ± 0.0001 | +0.0062 ± 0.0002 | +0.0170 ± 0.0005 | +0.0136 ± 0.0004 |
| NNjBC − NNj | the beam connection, given the constituents (above) | +0.0129 ± 0.0003 | +0.0133 ± 0.0002 | +0.0322 ± 0.0001 | +0.0287 ± 0.0001 | +0.0148 ± 0.0008 | +0.0136 ± 0.0002 |

- **Kj alone has most of NNj's AUC.** On gate jets the constituents add 0.010 (Herwig) and 0.008 (MG5+Pythia) to Kj's 0.753 and 0.756. The beam connection alone adds more: 0.014 on both generators.
- **Both gains are larger for the VBF-preselected jets**: +0.020 / +0.017 for the constituents and +0.018 / +0.017 for the beam connection. The preselection has used up part of the kinematic separation, so the rest counts for more.
- **The two are nearly additive.** On Herwig gate jets, the constituents are worth 0.0087 with the beam connection and 0.0096 without, and the beam connection 0.0129 with the constituents and 0.0138 without.
- **Signal efficiencies at equal background efficiency:** [roc.summary.md](../roc.summary.md), section 5.
