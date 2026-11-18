# nov26 jet networks: training sample, five-fold cross-fitting and scores

Two jet taggers are trained on the nov26 friends to tell a jet from VBF Z production (signal) from a jet in QCD Z+jj production (background):

- **NNjPF**: the DeepSets net NNj from summer26, with each constituent's particle type (CMS PF categories) and charge added to its inputs.
- **ParTj**: a small Particle Transformer. Its attention also uses pairwise features of every pair of constituents.

Both see one jet at a time. The jet's inputs are |η|, m and pT. Each constituent contributes the signed Δη (`jcsDEta`), ΔΦ, pT/pT_jet, its PF type and its charge. Both nets are trained on Herwig only; MG5+Pythia is used purely to test how well they transfer to another generator.

This note explains how the training sample is built, why it is weighted, and how five-fold cross-fitting gives every event a score that no network has seen. The code is in this folder: `common.py` (sample, folds, weights), `NNjPF/model.py`, `ParTj/model.py`, `train_fold.py`, `save_scores.py` and `roc_scores.py`.

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

The assignment depends only on the file content and the seed, so all jobs and both networks get the same folds. It is stored in the h5 files as `subfold`.

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

| | NNjPF | ParTj |
|---|---|---|
| Architecture | DeepSets: φ 9 → 64 → 64, masked sum; ρ (64 + 3) → 128 → 64 → 1 | ParT: embedding 9 → 128 → 64; pair features ln kT, ln z, ln ΔR, ln m² → 4 → 32 → 32 → 4 attention biases; 3 particle-attention blocks; 2 class-attention blocks; head (64 + 3) → 64 → 1 |
| Parameters | 22 209 | 265 421 |
| Optimiser | Adam, learning rate 1e-3 | AdamW, learning rate 5e-4, weight decay 0.01, gradient clip 1.0 |
| Batch size | 512 | 256 |
| Learning-rate schedule | ×0.5 after 5 epochs without improvement | ×0.5 after 3 epochs without improvement |
| Early stopping | 15 epochs without improvement, at most 100 epochs | 8 epochs, at most 40 epochs |

Common to both:
- The PF type enters as a one-hot of the 5 categories, and the charge as a number.
- Each fold fits its own StandardScaler to its training jets. It covers the jet |η|, m and pT, and the continuous constituent features of valid constituents.
- The loss is binary cross-entropy weighted as in section 3.
- The constituent axis is trimmed to the largest constituent count in each batch.
- The networks run as condor GPU jobs: `NNjPF/train_gpu.sub` and `ParTj/train_gpu.sub`, one job per fold.

Each job writes `<net>/best_model_fold{k}.pt`, `scaler_fold{k}.pkl`, `history_fold{k}.npz` and `scores_fold{k}.npz`, plus a loss curve `figs/nov26/<net>_fold{k}_loss.pdf`.

## 8. The scores in the h5 files

After all five jobs of a net finish, `python3 nov26/save_scores.py --net <net>` writes into each of the four `friends/nov26/*.h5`:

| Dataset | Shape | Content |
|---|---|---|
| `NNjPF`, `ParTj` | (N, 2) float32 | out-of-fold score (sigmoid, 0–1) of the jet in each slot; −1 for an empty slot |
| `subfold` | (N,) int8 | the event's sub-fold; network k scored sub-folds 2k and 2k+1 |

Before writing, `save_scores.py` checks that every stored jet was scored exactly once, by the right network.

For example, to get the scores of both jets of the VBF-preselected gate events:

```python
import h5py, numpy as np
with h5py.File('friends/nov26/VBFZ_herwig.h5') as f:
    sel = f['gate'][:] & (f['mjj'][:] > 400) & (np.abs(f['dYjj'][:]) > 2.5)
    s = f['NNjPF'][:][sel]          # (n, 2): both jets of every selected event, all out-of-fold
```

Rerunning `root_to_hdf5_nov26.py` rewrites the h5 files and drops these datasets. Afterwards, rerun `save_scores.py` for each net; no retraining is needed.

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
    gate, lead, s = f['gate'][:], f['leadJetIndex'][:], f['NNjPF'][:]
rows = np.flatnonzero(gate)              # gate events always store two jets
s_lead = s[rows, lead[rows]]             # score of the pT-leading jet
s_sub = s[rows, 1 - lead[rows]]          # score of the subleading jet
```

Both arrays are out-of-fold scores. They can be combined into an event-level score: a product, a mean, or a small second-level model. A second-level model trained on these scores has to reuse the same sub-folds (the `subfold` column) so that it stays honest too.

## 10. Results

### NNjPF (all five folds finished, 2026-11-17)

**Training.** Each fold model trained on about 1.47M jets and was validated on about 211k.

| Fold | GPU | Stopped at epoch | Best validation loss | Wall time |
|---|---|---|---|---|
| 0 | V100S | 60 | 0.6069 | 13 min |
| 1 | V100S | 96 | 0.6075 | 19 min |
| 2 | T4 | 61 | 0.6096 | 13 min |
| 3 | H100 (MIG slice) | 73 | 0.6086 | 15 min |
| 4 | V100S | 63 | 0.6077 | 12 min |

"Stopped at epoch" is where early stopping triggered; the best epoch is 15 earlier. Fold 3 ran on a GPU newer than the LCG_106 PyTorch build supports natively, so its kernels were compiled at start-up and its first epoch took 5 minutes instead of 16 s. Its results agree with the other folds.

**Figures of merit.** Every event is scored by the network that never saw it. Two numbers summarise the ROC curve:
- **AUC**, the area under the ROC curve: the probability that a random VBF jet scores higher than a random QCD jet.
- **ε_S(1 %)**, the signal efficiency at 1 % background efficiency. Cut on the score so that only 1 % of the QCD jets survive; ε_S(1 %) is the fraction of VBF jets that survive the same cut. It measures the high-purity end of the ROC curve, which the AUC barely weighs.

| Jets | AUC Herwig | AUC MG5+Pythia | ε_S(1 %) Herwig | ε_S(1 %) MG5+Pythia |
|---|---|---|---|---|
| Gate jets (roles A + B) | 0.7643 | 0.7637 | 0.0983 | 0.1028 |
| 1st jet of gate events (role A) | 0.7640 | 0.7634 | 0.0980 | 0.1009 |
| 2nd jet of gate events (role B) | 0.7646 | 0.7640 | 0.0990 | 0.1044 |
| Gate == 0 jets (role C) | 0.7199 | 0.7094 | 0.0632 | 0.0596 |
| Gate jets of VBF-preselected events (mjj > 400, \|ΔY\| > 2.5) | 0.6877 | 0.6961 | 0.0497 | 0.0575 |

- **1st vs 2nd jet:** the AUCs agree within 0.0006. For ε_S(1 %) a bootstrap (`nov26/bootstrap_scores.py`, which resamples whole events) gives a statistical uncertainty of about 0.0016 for the 1st or the 2nd jets alone, and 0.0022–0.0024 for their difference. The measured differences are 0.0010 (Herwig, 0.5σ) and 0.0035 (Pythia, 1.5σ). This is the consistency check of section 9: no sign of leakage or of a slot bias.
- **Spread of the five networks** on their own Herwig gate jets: AUC 0.7632, 0.7665, 0.7659, 0.7634, 0.7627 (rms 0.0016); ε_S(1 %) 0.0969, 0.0980, 0.1010, 0.0988, 0.0952 (rms 0.0019).
- **Generator transfer:** the Herwig-trained networks lose nothing on Pythia for gate jets.
- **VBF-preselected events:** the AUC is lower here because the preselection has already removed the easiest kinematic differences.

**Comparison with summer26 NNj.** NNj was trained and tested on slot-0 jets of VBF-preselected events, built from all particles with no pT cut. On a 15 % test split it reached an AUC of 0.685 (Herwig) and 0.699 (MG5+Pythia). NNjPF in the same region gives 0.688 and 0.696.

The two are compatible. The statistical uncertainty is about 0.008 for NNj's test set and about 0.002 for the nov26 out-of-fold scores. The setups differ in three ways (particle threshold, PF inputs, training population), so the near-equality does not isolate the value of the PF information.

### ParTj (all five folds finished, 2026-11-17)

**Training.** ParTj used the same folds as NNjPF, so each fold model trained and validated on the same jets as the NNjPF model of that fold.

| Fold | GPU | Stopped at epoch | Best validation loss | Wall time |
|---|---|---|---|---|
| 0 | H100 (MIG slice) | 29 | 0.6058 | 51 min |
| 1 | H100 (MIG slice) | 33 | 0.6062 | 50 min |
| 2 | H100 NVL | 36 | 0.6083 | 42 min |
| 3 | H100 NVL | 33 | 0.6072 | 39 min |
| 4 | T4 | 35 | 0.6065 | 107 min |

- The best epoch is 8 before the stop.
- An epoch takes about 1 minute on a full H100, 1.5 minutes on a MIG slice and 3 minutes on a T4.
- In every fold, the best validation loss is 0.0011–0.0014 lower than NNjPF's on the same validation jets.

**Figures of merit.**

| Jets | AUC Herwig | AUC MG5+Pythia | ε_S(1 %) Herwig | ε_S(1 %) MG5+Pythia |
|---|---|---|---|---|
| Gate jets (roles A + B) | 0.7656 | 0.7634 | 0.0970 | 0.1015 |
| 1st jet of gate events (role A) | 0.7653 | 0.7631 | 0.0956 | 0.1005 |
| 2nd jet of gate events (role B) | 0.7659 | 0.7637 | 0.0984 | 0.1022 |
| Gate == 0 jets (role C) | 0.7214 | 0.7100 | 0.0650 | 0.0609 |
| Gate jets of VBF-preselected events (mjj > 400, \|ΔY\| > 2.5) | 0.6913 | 0.6969 | 0.0531 | 0.0589 |

- **1st vs 2nd jet:** the AUCs agree within 0.0006, and ε_S(1 %) within 0.0028 (at most 1.3σ), as for NNjPF.
- **Spread of the five networks** on their own Herwig gate jets: AUC 0.7649, 0.7674, 0.7675, 0.7643, 0.7640 (rms 0.0015); ε_S(1 %) 0.0956, 0.0945, 0.1018, 0.0974, 0.0964 (rms 0.0025).

### ParTj vs NNjPF

Both nets scored exactly the same jets, so they can be compared jet by jet. This is a *paired* comparison:

- For each fold model k, take the AUC difference ParTj − NNjPF on model k's inference jets.
- Each AUC on its own fluctuates by about 0.001, from the finite number of jets. The two nets' scores are 98 % correlated, so this fluctuation is nearly the same for both and cancels in the difference.
- The five inference sets are disjoint, so the five differences are independent. The table gives their mean ± its standard error.

| Jets | AUC, Herwig | AUC, MG5+Pythia | ε_S(1 %), Herwig | ε_S(1 %), MG5+Pythia |
|---|---|---|---|---|
| Gate jets | +0.0013 ± 0.0002 | −0.0003 ± 0.0005 | −0.0008 ± 0.0009 | −0.0009 ± 0.0015 |
| Gate == 0 jets | +0.0016 ± 0.0002 | +0.0006 ± 0.0003 | +0.0019 ± 0.0004 | +0.0013 ± 0.0004 |
| Gate jets of VBF-preselected events | +0.0036 ± 0.0003 | +0.0008 ± 0.0006 | +0.0047 ± 0.0016 | +0.0022 ± 0.0014 |

**ROC curves** (`figs/nov26/roc_nets.pdf`): one column each for the gate jets, the gate == 0 jets and the VBF-preselected gate jets.
- Top panels: background rejection 1/ε_B against signal efficiency, for both nets and both generators. A dotted line marks ε_B = 1 %, and the legends give the AUC and ε_S(1 %).
- Lower panels: the ratio of the two nets' rejections at equal signal efficiency.

Background rejection of the VBF-preselected gate jets at fixed signal efficiency:

| Signal efficiency | 5 % | 10 % | 20 % | 30 % | 50 % | 70 % | 90 % |
|---|---|---|---|---|---|---|---|
| NNjPF rejection, Herwig | 97.8 | 41.0 | 15.7 | 8.8 | 4.1 | 2.3 | 1.4 |
| ParTj / NNjPF, Herwig | 1.116 | 1.054 | 1.029 | 1.022 | 1.015 | 1.011 | 1.007 |
| NNjPF rejection, MG5+Pythia | 120.2 | 45.3 | 17.3 | 9.5 | 4.3 | 2.4 | 1.4 |
| ParTj / NNjPF, MG5+Pythia | 1.022 | 1.022 | 1.007 | 1.002 | 0.999 | 1.002 | 1.002 |

- **On Herwig**, the generator both nets were trained on, ParTj has the higher AUC in all five fold models and for all three groups of jets. The gain is small. Its largest relative size is at high purity for the VBF-preselected jets: +12 % rejection at 5 % signal efficiency. Only about 450 background jets survive there, so that point is the least precise.
- **At the 1 % working point** ParTj keeps a gain for the gate == 0 jets and the VBF-preselected jets, but not for the gate jets as a whole (−0.0008 ± 0.0009).
- **On MG5+Pythia** most of the gain disappears. This suggests that much of what ParTj adds through the constituent pairs is specific to how Herwig showers and hadronises, rather than a VBF-vs-QCD difference that both generators share.
- **Cost:** ParTj has 12 times more parameters (265k vs 22k). It took 39–107 minutes per fold, against 12–19 for NNjPF (on different GPUs).

**Scores in the h5.** Both nets' scores are now in the four `friends/nov26/*.h5`: datasets `NNjPF`, `ParTj` and `subfold` (section 8). They were written after all ten jobs had finished. Read back, they are identical to the fold files.
