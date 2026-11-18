# fall26 jet networks: NNjPF and ParTj on the Fall2026 Z samples, with event weights

NNjPF and ParTj are trained here on the Fall2026 Z samples ([fall26.md](../fall26.md)). Everything is done as in nov26, which [nov26/training.md](../nov26/training.md) explains step by step. This note covers only what changes, the event weights, and the results.

The code is the nov26 code with `--campaign fall26`:

```bash
python3 nov26/train_fold.py --campaign fall26 --net NNjPF --fold 0    # one fold model (condor: fall26/<net>/train_gpu.sub)
python3 nov26/save_scores.py --campaign fall26 --net NNjPF            # out-of-fold scores into friends/fall26/*.h5
python3 nov26/roc_scores.py --campaign fall26                         # event-weighted AUC and ε_S(1 %) tables, ROC figure
python3 nov26/bootstrap_scores.py --campaign fall26                   # their statistical uncertainties (bootstrap)

# generator swap (section 6): the same net trained on MG5+Pythia
python3 nov26/train_fold.py --campaign fall26 --net NNjPF --fold 0 --train-gen pythia   # condor: fall26/NNjPF_py/train_gpu.sub
python3 nov26/save_scores.py --campaign fall26 --net NNjPF --train-gen pythia           # h5 column NNjPF_py
python3 nov26/train_gen_matrix.py --campaign fall26                                     # trained-on x tested-on tables, ROC figure
```

Models, scores and logs go to `fall26/<net>/` (`fall26/<net>_py/` when trained on MG5+Pythia), and figures to `figs/fall26/`.

---

## 1. What stays as in nov26

- **One example is one jet.**
  - Role A (the 1st jet of a gate event) and role C (every jet of a gate == 0 event) are trained on.
  - Role B (the 2nd jet of a gate event) is only scored.
  - Training uses Herwig only; MG5+Pythia tests the transfer to another generator.
- **Five-fold cross-fitting.**
  - There are 10 sub-folds, stratified by (sample, gate, nJets), with seed 2026.
  - Fold model k scores sub-folds 2k and 2k+1, validates on (2k+2) mod 10, and trains on the other seven.
  - Every jet is scored once, by a model that never saw its event.
- **Inputs, architectures, hyper-parameters, early stopping:** unchanged (nov26/training.md, section 7).
- **Output:** the datasets `NNjPF` and `ParTj`, shape (N, 2), and `subfold`, shape (N,), in `friends/fall26/<sample>.h5`.

## 2. Why the event weights now matter

Each generated event carries a weight. An unweighted sample is one in which every event counts the same. Herwig generated the Fall2026 VBF Z sample with **variable** weights, between 1 and 150. Phase-space regions that were sampled more densely get smaller weights, and sparsely sampled regions get larger ones. Only the weighted sample describes the physics.

The *effective number of events*, N_eff = (Σw)²/Σw², measures how much statistical power is left after the weights. For the Herwig samples:

| Sample | nov26 inputs, N_eff / N | Fall2026, N_eff / N |
|---|---|---|
| VBFZ Herwig | about 0.96 | 0.36 |
| QCDZjj Herwig | about 0.88 | 0.89 |

In nov26 the weights were close to uniform, so training on unit weights changed little.

In Fall2026, an unweighted training on VBFZ Herwig would learn the generator's sampling density instead of the physical distribution. The unweighted AUC would then measure the separation in that distorted population. The package README states the rule plainly: weights must never be replaced by unit weights.

So in fall26:
- **the loss is weighted** with the event weights (section 3);
- **the AUCs and ROC curves are weighted** with `kWeight`. The unweighted AUCs are printed too, only for comparison with nov26.

## 3. The weighted class balance

nov26 gave every jet of class c in population P (gate jets, or gate == 0 jets) the weight

  w = N_P / (2 N_{P,c}),

so the two classes weigh 50/50 inside each population and each population keeps its size N_P (nov26/training.md, section 3).

fall26 keeps exactly these totals, but shares each class's total among its jets in proportion to their event weights:

  w_i = kWeight_i × N_P / (2 Σ_{j ∈ (P,c)} kWeight_j).

- **The totals per class and per population are the same as in nov26:** 50/50 inside each population, and each population keeps N_P.
- **Inside a class,** a jet from an event with twice the event weight counts twice as much, as the physics requires.

The price is statistical. A few jets with large weights carry a large share of the loss, and the training is noisier. The VBF gate jets (role A) are worth only about a quarter of their count, and the VBF gate == 0 jets a little over half (section 5). The fold logs print N_eff for each class and role.

## 4. Running the jobs

```bash
source init.source
condor_submit -name bigbird24.cern.ch fall26/NNjPF/train_gpu.sub
condor_submit -name bigbird24.cern.ch fall26/ParTj/train_gpu.sub
```

**GPU types.** The LCG_106_cuda PyTorch build has no kernels for compute capability 8.0. On an A100 (or an A30) every job fails at the first GPU operation with "no kernel image is available for execution on the device". The first ParTj submission (cluster 1647771) landed on an A100 node and lost all five jobs this way.

All four submit files (nov26 and fall26) now exclude such GPUs:
```
require_gpus = (Capability =!= 8.0) && !regexp("A100|A30", DeviceName)
```
- V100, T4 and H100 work. An H100 (sm_90) compiles its kernels at start-up: the first epoch of the first nov26 jobs on H100s took 5–9 minutes. The fall26 ParTj jobs on H100s lost only 7–32 s, probably because the compiled kernels are cached in the home directory (`~/.nv/ComputeCache`, last written during the nov26 jobs).
- The meta-comparison `=!=` keeps H100 MIG slices eligible: they do not advertise a capability.

## 5. Results

### NNjPF (all five folds finished, 2026-11-18)

**Training.**

| Fold | GPU | Stopped at epoch | Best validation loss | Wall time |
|---|---|---|---|---|
| 0 | V100S | 58 | 0.6109 | 12 min |
| 1 | V100 | 46 | 0.6103 | 14 min |
| 2 | V100S | 51 | 0.6104 | 12 min |
| 3 | V100S | 45 | 0.6088 | 10 min |
| 4 | V100 | 59 | 0.6081 | 13 min |

The validation loss is event-weighted here, so it cannot be compared with the nov26 values.

**Training statistics with weights.** In fold 0, which is typical:

| Training jets | Jets | N_eff | N_eff / N |
|---|---|---|---|
| VBF, role A (gate) | 239 673 | 54 799 | 0.23 |
| VBF, role C (gate == 0) | 556 024 | 314 525 | 0.57 |
| QCD, role A | 104 027 | 99 941 | 0.96 |
| QCD, role C | 588 736 | 552 713 | 0.94 |

The large VBF Herwig weights sit in the gate events, the region the analysis cares about.

**Figures of merit** (event-weighted; for Pythia, whose weights are constant, the weighted and unweighted values are the same):

| Jets | AUC Herwig | AUC MG5+Pythia | ε_S(1 %) Herwig | ε_S(1 %) MG5+Pythia |
|---|---|---|---|---|
| Gate jets (roles A + B) | 0.7666 | 0.7530 | 0.1006 | 0.0908 |
| 1st jet of gate events (role A) | 0.7660 | 0.7531 | 0.1001 | 0.0900 |
| 2nd jet of gate events (role B) | 0.7671 | 0.7528 | 0.1009 | 0.0914 |
| Gate == 0 jets (role C) | 0.7170 | 0.7019 | 0.0625 | 0.0543 |
| Gate jets of VBF-preselected events | 0.6918 | 0.6791 | 0.0514 | 0.0484 |

- **1st vs 2nd jet:** the AUCs agree within 0.0011 and ε_S(1 %) within 0.0014. The statistical uncertainties of these differences are 0.0009–0.0010 and 0.0024–0.0029 (bootstrap, see *Statistical uncertainties* below). No sign of leakage or of a slot bias.
- **Spread of the five networks** on their own Herwig gate jets:
  - AUC: 0.7665, 0.7681, 0.7659, 0.7646, 0.7678 (rms 0.0013);
  - ε_S(1 %): 0.0979, 0.1057, 0.1013, 0.0949, 0.1078 (rms 0.0048).
  
  The ε_S(1 %) spread is 2.5 times that of nov26: only a quarter of the VBF gate jets' statistical power remains (N_eff 163 000 of 685 000 jets). The spread is what these statistics allow: each model's ε_S(1 %) on its own jets is uncertain by 0.003–0.006 (χ² = 6.7 for 4 degrees of freedom).
- **The weights change the answer.** Unweighted, the Herwig gate-jet AUC would be 0.7312 and ε_S(1 %) 0.0762. The densely sampled regions, which carry small weights, are the harder ones to separate. Only the weighted numbers describe the physics.
- **Herwig vs MG5+Pythia.** The Herwig-trained networks separate Pythia jets less well: the AUC is lower by 0.014 for gate jets, and ε_S(1 %) by 0.010 (about 14σ and 4.5σ with the bootstrap uncertainties). In nov26 there was no such gap. Two explanations are possible: a transfer loss, where the nets learned Herwig-specific features, or a genuine difference in how separable the jets are in the new Pythia samples (whose QCD gate fraction also changed, see fall26.md). Only a network trained on Pythia can tell them apart; section 6 does this. At ε_S(1 %) the loss is entirely a transfer loss, and in AUC about half of it is.
- **Compared with nov26 NNjPF** on Herwig gate jets: AUC 0.7666 (fall26, weighted) against 0.7643, and ε_S(1 %) 0.101 against 0.098. The two are similar, but the samples (generator versions, generation cuts) and the weighting differ.

### ParTj (all five folds finished, 2026-11-18)

**Training.** ParTj used the same folds as NNjPF, so each fold model trained and validated on the same jets as the NNjPF model of that fold. The resubmission (cluster 1647780, section 4) ran on four H100 NVL GPUs and one V100.

| Fold | GPU | Stopped at epoch | Best validation loss | Wall time |
|---|---|---|---|---|
| 0 | H100 NVL | 35 | 0.6092 | 39 min |
| 1 | H100 NVL | 33 | 0.6086 | 42 min |
| 2 | H100 NVL | 31 | 0.6094 | 39 min |
| 3 | V100 | 29 | 0.6074 | 78 min |
| 4 | H100 NVL | 28 | 0.6071 | 30 min |

- The best epoch is 8 before the stop.
- An epoch takes about 1 minute on an H100 NVL and 2.6 minutes on the V100.
- In every fold, the best validation loss is 0.0010–0.0017 lower than NNjPF's on the same validation jets (nov26: 0.0011–0.0014).

**Figures of merit** (event-weighted):

| Jets | AUC Herwig | AUC MG5+Pythia | ε_S(1 %) Herwig | ε_S(1 %) MG5+Pythia |
|---|---|---|---|---|
| Gate jets (roles A + B) | 0.7686 | 0.7544 | 0.1015 | 0.0908 |
| 1st jet of gate events (role A) | 0.7680 | 0.7546 | 0.0989 | 0.0900 |
| 2nd jet of gate events (role B) | 0.7693 | 0.7542 | 0.1046 | 0.0916 |
| Gate == 0 jets (role C) | 0.7185 | 0.7027 | 0.0633 | 0.0549 |
| Gate jets of VBF-preselected events | 0.6952 | 0.6813 | 0.0544 | 0.0519 |

- **1st vs 2nd jet:** the AUCs agree within 0.0013 (uncertainty of the difference: 0.0010). For ε_S(1 %) the 2nd jet is ahead by 0.0056 on Herwig, 1.8σ with the bootstrap uncertainty of 0.0032, and by 0.0016 on Pythia (0.6σ). Both are compatible with statistical fluctuations.
- **Spread of the five networks** on their own Herwig gate jets: AUC 0.7678, 0.7705, 0.7687, 0.7661, 0.7702 (rms 0.0016); ε_S(1 %) 0.0933, 0.1105, 0.0981, 0.1004, 0.1117 (rms 0.0072).
  - The AUC spread is statistical.
  - The ε_S(1 %) values scatter a little more than their statistical uncertainties (0.003–0.007) allow: χ² = 10.2 for 4 degrees of freedom (p = 0.04), against 6.7 (p = 0.15) for NNjPF. At the 1 % working point the five ParTj models seem to differ slightly from one another.
- **Unweighted**, for comparison only: Herwig gate-jet AUC 0.7339 and ε_S(1 %) 0.0766.
- **Compared with nov26 ParTj** on Herwig gate jets: AUC 0.7686 (fall26, weighted) against 0.7656, and ε_S(1 %) 0.102 against 0.097.

### Statistical uncertainties

How much would these numbers change with another sample of the same size? A *bootstrap* estimates this from the sample itself (`python3 nov26/bootstrap_scores.py --campaign fall26`):

1. Make a replica of the sample. Each event gets a random multiplicity drawn from a Poisson distribution of mean 1 (0, 1, 2, … copies), which multiplies its event weight.
2. Recompute every figure of merit on the replica.
3. Repeat 400 times. The standard deviation of the 400 values is the statistical uncertainty.

Three details matter:
- **The multiplicity belongs to the event, not to the jet.** Both jets of a gate event share one event weight, so they must go up and down together. Resampling the jets one by one would treat them as independent and understate the uncertainty.
- **Both nets see the same replicas.** The bootstrap then also gives the uncertainty of their difference, in which most of the sample fluctuation cancels.
- **The trained models stay fixed.** The bootstrap leaves out the randomness of the training, which the spread of the five fold models includes.

Statistical uncertainties (ranges cover the two nets, and in the last row the five models):

| Jets | AUC, Herwig | AUC, MG5+Pythia | ε_S(1 %), Herwig | ε_S(1 %), MG5+Pythia |
|---|---|---|---|---|
| Gate jets | 0.0008 | 0.0006 | 0.0018–0.0019 | 0.0012 |
| 1st jet or 2nd jet alone | 0.0010 | 0.0007–0.0008 | 0.0023–0.0027 | 0.0015–0.0019 |
| Gate == 0 jets | 0.0005 | 0.0004 | 0.0007 | 0.0006 |
| Gate jets of VBF-preselected events | 0.0016–0.0017 | 0.0013 | 0.0021–0.0023 | 0.0020–0.0023 |
| One fold model's own gate jets | 0.0017–0.0021 | 0.0012–0.0014 | 0.0031–0.0066 | 0.0023–0.0032 |

- **The event weights cost precision at the 1 % working point.** For the Herwig gate jets ε_S(1 %) is uncertain by 0.0018–0.0019, against 0.0012–0.0014 in nov26 (same script, unweighted). MG5+Pythia, whose weights are constant, is as precise as nov26.
- **Where the uncertainty comes from.** ε_S(1 %) moves when the signal sample fluctuates. It also moves when the background sample fluctuates, because the cut that keeps 1 % of the QCD jets moves with it. Resampling one sample at a time (NNjPF, gate jets):

  | Samples | Signal sample only | Background sample only |
  |---|---|---|
  | nov26 (both generators), fall26 MG5+Pythia | 0.0003–0.0004 | 0.0012–0.0013 |
  | fall26 Herwig | 0.0013 | 0.0011 |

  With uniform weights the background dominates. The Herwig VBF weights triple the signal part, which becomes as large as the background part.
- **Why triple?** The overall N_eff/N = 0.24 of the VBF gate jets alone would only double it (1/√0.24 ≈ 2). But the large weights sit among the most VBF-like jets:
  - the VBF jets above the 1 % cut weigh 30 % more on average than all VBF gate jets (weighted efficiency 0.101, unweighted 0.077);
  - their weights are also more uneven (N_eff/N = 0.10).

  This is also why the weighted figures of merit are higher than the unweighted ones.

### ParTj vs NNjPF

As in nov26 ([nov26/training.md](../nov26/training.md), *ParTj vs NNjPF*), the comparison is paired:
- For each fold model k, take the difference ParTj − NNjPF on model k's inference jets.
- Both nets scored the same jets, so most of the sample fluctuation cancels in the difference.
- The table gives the mean of the five differences.

The error is the larger of two estimates:
- **The standard error of the five differences**, as in nov26. It includes the training randomness, but five values give only a rough estimate (for one entry it came out below 0.0001).
- **The bootstrap uncertainty of the difference.** It is precise, but leaves out the training randomness.

| Jets | AUC, Herwig | AUC, MG5+Pythia | ε_S(1 %), Herwig | ε_S(1 %), MG5+Pythia |
|---|---|---|---|---|
| Gate jets | +0.0021 ± 0.0003 | +0.0014 ± 0.0004 | +0.0013 ± 0.0021 | +0.0000 ± 0.0016 |
| Gate == 0 jets | +0.0015 ± 0.0001 | +0.0009 ± 0.0003 | +0.0005 ± 0.0005 | +0.0009 ± 0.0004 |
| Gate jets of VBF-preselected events | +0.0033 ± 0.0004 | +0.0023 ± 0.0004 | +0.0041 ± 0.0018 | +0.0010 ± 0.0025 |

**ROC curves** (`figs/fall26/roc_nets.pdf`): event-weighted, with the same layout as in nov26.
- There is one column each for the gate jets, the gate == 0 jets and the VBF-preselected gate jets.
- Top panels: background rejection 1/ε_B against signal efficiency, for both nets and both generators. A dotted line marks ε_B = 1 %, and the legends give the AUC and ε_S(1 %).
- Lower panels: the ratio of the two nets' rejections at equal signal efficiency.

Background rejection of the VBF-preselected gate jets at fixed signal efficiency:

| Signal efficiency | 5 % | 10 % | 20 % | 30 % | 50 % | 70 % | 90 % |
|---|---|---|---|---|---|---|---|
| NNjPF rejection, Herwig | 104.2 | 41.7 | 16.0 | 8.9 | 4.2 | 2.4 | 1.4 |
| ParTj / NNjPF, Herwig | 1.092 | 1.089 | 1.030 | 1.027 | 1.015 | 1.008 | 1.005 |
| NNjPF rejection, MG5+Pythia | 96.0 | 37.7 | 14.4 | 8.2 | 3.9 | 2.2 | 1.4 |
| ParTj / NNjPF, MG5+Pythia | 1.104 | 1.042 | 1.050 | 1.019 | 1.002 | 1.010 | 1.004 |

- **AUC:** ParTj is ahead for all three groups of jets and on both generators, in 29 of the 30 per-model comparisons. On Herwig the gain for the gate jets is +0.0021 ± 0.0003.
- **The AUC gain carries over to MG5+Pythia:** about two thirds of it remains (+0.0014 ± 0.0004 for the gate jets).
  - In nov26 it vanished on Pythia (−0.0003 ± 0.0005), which suggested that ParTj's extra came from Herwig-specific features. That suggestion does not hold for the Fall2026 samples.
  - The two campaigns differ in several ways (generator versions, generation cuts, event-weighted training), so the change cannot be pinned on one of them.
- **At the 1 % working point** the gains are small, at most about 2σ:
  - VBF-preselected jets: +0.0041 ± 0.0018 on Herwig (+8 % relative; all five models agree), and +0.0010 ± 0.0025 on Pythia (the five models scatter from −0.0064 to +0.0086).
  - Gate jets as a whole: no gain (+0.0013 ± 0.0021 and +0.0000 ± 0.0016), as in nov26.
  - Gate == 0 jets: +0.0005 ± 0.0005 and +0.0009 ± 0.0004.
- **Rejection at fixed signal efficiency:** for the VBF-preselected jets ParTj rejects 9–10 % more background at 5 % signal efficiency on both generators, and 1–2 % more at 50 %. Only about 450 (Herwig) and 410 (Pythia) QCD jets survive at 5 % signal efficiency, so that column is the least precise.
- **Generator transfer:** both nets lose the same on MG5+Pythia. For the gate jets the AUC drops by 0.014 (about 14σ), and ε_S(1 %) by 0.010 (NNjPF, 4.5σ) and 0.011 (ParTj, 4.8σ). The gap belongs to the samples or to the Herwig training, not to one architecture. Section 6 splits it between the two.
- **Cost:** ParTj has 12 times more parameters (265k vs 22k). It took 30–78 minutes per fold, against 10–14 for NNjPF.

**Scores in the h5.** Both nets' scores are now in the four `friends/fall26/*.h5`: datasets `NNjPF`, `ParTj` and `subfold`, as in nov26. They were written after all ten jobs had finished. Read back, they are identical to the fold files:
- every stored jet has exactly one score, and empty slots hold −1;
- each event was scored by the model whose inference sub-folds contain it.

---

## 6. NNjPF trained on MG5+Pythia: where the Herwig → Pythia loss comes from

Section 5 found that the Herwig-trained nets separate MG5+Pythia jets less well than Herwig jets (gate jets: AUC −0.014, ε_S(1 %) −0.010). In nov26 they did not. There are two possible reasons:
- a **transfer loss:** the nets learned features that only Herwig has;
- a **sample difference:** the MG5+Pythia jets are harder to separate, whatever the training.

To tell them apart, NNjPF was trained a second time, on MG5+Pythia (2026-11-18), and both versions were tested on both generators.

**Setup.** Only the training generator changes (`--train-gen pythia`; commands at the top of this note):
- **The same sub-folds:** model k of either version scored the same events. The sub-folds are assigned per sample, so the Pythia-trained models are out-of-fold on every Pythia jet, and they never saw a Herwig jet.
- **The same architecture, hyperparameters and class-balanced loss.** The Pythia weights are constant, so the event weights change nothing in this training.
- **Outputs:** `fall26/NNjPF_py/`, and the h5 column `NNjPF_py` in all four `friends/fall26/*.h5`. Read back, it is identical to the fold files, and the `NNjPF` and `ParTj` columns are unchanged.

**Training** (condor cluster 1647840):

| Fold | GPU | Stopped at epoch | Best validation loss | Wall time |
|---|---|---|---|---|
| 0 | V100S | 61 | 0.6177 | 14 min |
| 1 | V100 | 57 | 0.6189 | 15 min |
| 2 | H100 (MIG slice) | 89 | 0.6163 | 18 min |
| 3 | H100 (MIG slice) | 47 | 0.6174 | 13 min |
| 4 | H100 (MIG slice) | 52 | 0.6147 | 13 min |

- **Fold 2's stdout file never reached EOS.** The job ended normally, and its model, loss history and scores did arrive. Its row comes from the loss history and the condor log.
- **With constant weights every training jet counts fully.** In fold 0 there are 260 050 VBF gate jets with N_eff = N, against N_eff = 54 799 for the Herwig training (section 5). The Pythia-trained net has 4.7 times more effective VBF gate jets.

**Figures of merit** (event-weighted, gate jets, ± bootstrap):

| Trained on | AUC, test Herwig | AUC, test MG5+Pythia | ε_S(1 %), test Herwig | ε_S(1 %), test MG5+Pythia |
|---|---|---|---|---|
| Herwig | 0.7666 ± 0.0008 | 0.7530 ± 0.0006 | 0.1006 ± 0.0018 | 0.0908 ± 0.0012 |
| MG5+Pythia | 0.7654 ± 0.0009 | 0.7603 ± 0.0006 | 0.1002 ± 0.0020 | 0.1019 ± 0.0012 |

**Splitting the loss.** The Herwig-trained net's loss on MG5+Pythia (Herwig-trained net: Herwig minus MG5+Pythia) is the sum of two parts:
- **The transfer loss on MG5+Pythia:** the Pythia-trained minus the Herwig-trained net, on the same MG5+Pythia jets. Both versions scored the same jets, so the comparison is paired and precise.
- **The matched gap:** each generator's jets scored by the net trained on that generator, Herwig minus MG5+Pythia. These are independent samples, so their uncertainties add in quadrature.

| Jets | Metric | Loss of the Herwig-trained net | = transfer loss | + matched gap |
|---|---|---|---|---|
| Gate jets | AUC | 0.0136 ± 0.0010 | 0.0073 ± 0.0002 | 0.0063 ± 0.0010 |
| Gate jets | ε_S(1 %) | 0.0098 ± 0.0022 | 0.0111 ± 0.0011 | −0.0012 ± 0.0022 |
| Gate == 0 jets | AUC | 0.0151 ± 0.0007 | 0.0032 ± 0.0001 | 0.0119 ± 0.0007 |
| Gate == 0 jets | ε_S(1 %) | 0.0083 ± 0.0009 | 0.0030 ± 0.0004 | 0.0052 ± 0.0009 |
| VBF-preselected gate jets | AUC | 0.0127 ± 0.0022 | 0.0063 ± 0.0006 | 0.0064 ± 0.0021 |
| VBF-preselected gate jets | ε_S(1 %) | 0.0029 ± 0.0029 | 0.0065 ± 0.0018 | −0.0034 ± 0.0029 |

- **Pooled values.** All values come from the pooled out-of-fold scores, so the parts add up (to within rounding).
- **Uncertainty of the transfer loss:** the larger of the bootstrap σ and the standard error of the five per-model differences. The five fold models agree; for the gate-jet AUC they range from 0.0070 to 0.0079.

The transfer in the other direction: the Pythia-trained net on Herwig jets (Herwig-trained minus Pythia-trained):

| Jets | AUC | ε_S(1 %) |
|---|---|---|
| Gate jets | 0.0011 ± 0.0003 | 0.0004 ± 0.0013 |
| Gate == 0 jets | 0.0040 ± 0.0001 | 0.0042 ± 0.0005 |
| VBF-preselected gate jets | 0.0071 ± 0.0008 | 0.0013 ± 0.0019 |

**Background rejection** 1/ε_B of the gate jets at fixed signal efficiency. "In-generator" means the net trained on the test generator; the ratio rows compare the net trained on the other generator with it.

| Signal efficiency | 5 % | 10 % | 20 % | 30 % | 50 % | 70 % | 90 % |
|---|---|---|---|---|---|---|---|
| Test Herwig: in-generator rejection | 272.9 | 100.9 | 33.8 | 17.1 | 6.8 | 3.3 | 1.6 |
| Pythia-trained / in-generator | 0.979 | 0.997 | 0.991 | 0.987 | 0.986 | 0.994 | 0.997 |
| Test MG5+Pythia: in-generator rejection | 292.5 | 103.0 | 33.2 | 16.1 | 6.4 | 3.2 | 1.6 |
| Herwig-trained / in-generator | 0.837 | 0.841 | 0.846 | 0.895 | 0.948 | 0.965 | 0.985 |

- **Precision:** at 5 % signal efficiency about 1100 (Herwig) and 900 (MG5+Pythia) QCD jets remain.
- **Figure** `figs/fall26/roc_train_gen_NNjPF.pdf`: the colour is the test generator; solid lines show the net trained on the test generator, dashed lines the net trained on the other. The lower panels show the rejection ratio.
- **The full printout** is in `fall26/NNjPF_py/train_gen_matrix.txt`.

**What this shows.**
- **At the working point the whole loss is a transfer loss.**
  - Trained on MG5+Pythia, NNjPF reaches ε_S(1 %) = 0.1019 ± 0.0012 on MG5+Pythia gate jets. That is as high as the Herwig-trained net on Herwig (matched gap −0.0012 ± 0.0022).
  - The Herwig-trained net loses 0.011 of it on MG5+Pythia, 11 % relative. At 5–20 % signal efficiency it rejects 15–16 % less QCD background.
- **In AUC the loss is half transfer, half sample.**
  - Of the 0.0136, 0.0073 ± 0.0002 is transfer, and 0.0063 ± 0.0010 remains with matched training.
  - The remaining gap sits at intermediate signal efficiencies. At 30–50 % the in-generator rejection is about 6 % lower on MG5+Pythia (16.1 against 17.1, 6.4 against 6.8); at 5–10 % it is as high or higher.
- **The transfer is asymmetric.** The Pythia-trained net loses only 0.0011 ± 0.0003 in AUC on Herwig gate jets, and nothing in ε_S(1 %).
- **Gate == 0 jets behave differently.**
  - Their transfer losses are small and similar in both directions (AUC 0.003–0.004).
  - Most of the Herwig-trained net's loss remains with matched training (AUC 0.0119 of 0.0151). Among gate == 0 jets the MG5+Pythia jets are harder to separate.
- **VBF-preselected gate jets:**
  - The transfer costs 0.006–0.007 in AUC in both directions.
  - At ε_S(1 %) only the Herwig-trained net loses (0.0065 ± 0.0018 on MG5+Pythia). The Pythia-trained net is equal within the uncertainty on Herwig (0.0013 ± 0.0019).

**Open points.**
- **The two trainings differ in statistical power, not only in generator.** The Pythia-trained net had 4.7 times more effective VBF gate jets. That alone could make it better and more robust, so part of the asymmetry may come from the statistics rather than from the generators. A Pythia training on a subsample with the Herwig N_eff would separate the two.
- **Why the fall26 Herwig training transfers worse than the nov26 one is not known.** The nov26 one showed no loss. The campaigns differ in:
  - the generator versions;
  - the generation cuts;
  - the event-weighted training, whose heavy Herwig VBF weights concentrate the training on fewer effective jets.

  A Herwig training without event weights would test the last.
