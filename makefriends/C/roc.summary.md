# ROC summary: NNj, NNjPF and ParTj in summer26, nov26 and fall26

The single-jet VBF taggers of the Z channel: VBF Z jets (signal) against QCD Z+jj jets (background).
- **NNj** (summer26): a DeepSets network on the jet's constituents.
- **NNjPF** (nov26, fall26): the same network, with each constituent's particle type and charge added.
- **ParTj** (nov26, fall26): a small Particle Transformer with the same inputs as NNjPF. It also uses pairwise features of every pair of constituents.

**Transformer against DeepSets:** section 4 has the ROC curves and, side by side, the AUC and ε_S(1 %) of both.

Three figures of merit:
- **AUC:** the area under the ROC curve.
- **ε_S(1 %):** the signal efficiency at 1 % background efficiency (higher is better).
- **ε_B(50 %):** the background efficiency at 50 % signal efficiency (lower is better).

All values are event-weighted (kWeight). The uncertainties are statistical, from a bootstrap with 400 replicas: each event gets a Poisson multiplicity, shared by its two jets. The numbers come from the scores stored in the h5 files (`python3 roc_summary.py`, 2026-11-18). NNjPF trained on MG5+Pythia ([fall26/training.md](fall26/training.md), section 6) is listed as an extra row.

---

## 1. VBF-preselected events: shared by all three campaigns

These are both jets of the events that pass the summer26 selection:
- two jets with pT > 30 GeV and |η| < 3, in opposite hemispheres (the nov26 gate);
- mjj > 400 GeV and |ΔY_jj| > 2.5.

The summer26 friends contain only these events, so this is the only category in which all three campaigns can be compared.

| Tested on | Campaign | Net | AUC | ε_S(1 %) | ε_B(50 %) |
|---|---|---|---|---|---|
| Herwig | summer26 | NNj (DeepSets) | 0.6891 ± 0.0036 | 0.0341 ± 0.0045 | 0.2411 ± 0.0056 |
| Herwig | nov26 | NNjPF (DeepSets) | 0.6889 ± 0.0013 | 0.0505 ± 0.0014 | 0.2436 ± 0.0022 |
| Herwig | nov26 | ParTj (transformer) | 0.6926 ± 0.0013 | 0.0540 ± 0.0018 | 0.2399 ± 0.0021 |
| Herwig | fall26 | NNjPF (DeepSets) | 0.6918 ± 0.0016 | 0.0514 ± 0.0023 | 0.2405 ± 0.0026 |
| Herwig | fall26 | ParTj (transformer) | 0.6952 ± 0.0015 | 0.0544 ± 0.0025 | 0.2369 ± 0.0025 |
| Herwig | fall26 | NNjPF (DeepSets), trained on MG5+Pythia | 0.6847 ± 0.0016 | 0.0501 ± 0.0024 | 0.2505 ± 0.0026 |
| MG5+Pythia | summer26 | NNj (DeepSets) | 0.6981 ± 0.0011 | 0.0570 ± 0.0019 | 0.2293 ± 0.0018 |
| MG5+Pythia | nov26 | NNjPF (DeepSets) | 0.6961 ± 0.0011 | 0.0575 ± 0.0015 | 0.2318 ± 0.0019 |
| MG5+Pythia | nov26 | ParTj (transformer) | 0.6969 ± 0.0011 | 0.0589 ± 0.0018 | 0.2321 ± 0.0018 |
| MG5+Pythia | fall26 | NNjPF (DeepSets) | 0.6791 ± 0.0013 | 0.0484 ± 0.0020 | 0.2585 ± 0.0022 |
| MG5+Pythia | fall26 | ParTj (transformer) | 0.6813 ± 0.0013 | 0.0519 ± 0.0021 | 0.2580 ± 0.0022 |
| MG5+Pythia | fall26 | NNjPF (DeepSets), trained on MG5+Pythia | 0.6854 ± 0.0013 | 0.0548 ± 0.0021 | 0.2497 ± 0.0022 |

- **summer26 NNj on Herwig** uses only its 15 % test split (30 788 VBF and 3 423 QCD events), hence the larger uncertainties.
- **Every other entry uses all events.** Each event is scored by a model that never saw it: five-fold cross-fitting for nov26 and fall26, and a net trained on the other generator otherwise.

## 2. All gate jets: nov26 and fall26

These are both jets of gate events, without the mjj and ΔY_jj cuts. summer26 has no entry, because its friends hold only VBF-preselected events.

| Tested on | Campaign | Net | AUC | ε_S(1 %) | ε_B(50 %) |
|---|---|---|---|---|---|
| Herwig | nov26 | NNjPF (DeepSets) | 0.7653 ± 0.0005 | 0.0998 ± 0.0012 | 0.1504 ± 0.0007 |
| Herwig | nov26 | ParTj (transformer) | 0.7666 ± 0.0005 | 0.0987 ± 0.0013 | 0.1483 ± 0.0008 |
| Herwig | fall26 | NNjPF (DeepSets) | 0.7666 ± 0.0008 | 0.1006 ± 0.0018 | 0.1477 ± 0.0011 |
| Herwig | fall26 | ParTj (transformer) | 0.7686 ± 0.0008 | 0.1015 ± 0.0020 | 0.1456 ± 0.0011 |
| Herwig | fall26 | NNjPF (DeepSets), trained on MG5+Pythia | 0.7654 ± 0.0008 | 0.1002 ± 0.0020 | 0.1499 ± 0.0011 |
| MG5+Pythia | nov26 | NNjPF (DeepSets) | 0.7637 ± 0.0005 | 0.1028 ± 0.0011 | 0.1459 ± 0.0007 |
| MG5+Pythia | nov26 | ParTj (transformer) | 0.7634 ± 0.0005 | 0.1015 ± 0.0012 | 0.1484 ± 0.0007 |
| MG5+Pythia | fall26 | NNjPF (DeepSets) | 0.7530 ± 0.0006 | 0.0908 ± 0.0013 | 0.1660 ± 0.0009 |
| MG5+Pythia | fall26 | ParTj (transformer) | 0.7544 ± 0.0006 | 0.0908 ± 0.0012 | 0.1649 ± 0.0008 |
| MG5+Pythia | fall26 | NNjPF (DeepSets), trained on MG5+Pythia | 0.7603 ± 0.0006 | 0.1019 ± 0.0012 | 0.1573 ± 0.0008 |

## 3. NNj against NNjPF and ParTj on the same events

summer26 and nov26 read the same generator files, so their nets can be compared on identical events, matched by the boson variables.
- **Events:** the summer26 evaluation events that also pass the selection in nov26, 99 % of the VBF and 95 % of the QCD events. That is 30 418 VBF and 3 255 QCD events on Herwig, and 205 125 and 29 038 on MG5+Pythia.
- **Bootstrap:** all three nets see the same replicas, so their differences are measured more precisely than the values themselves.
- **Jets:** the jets can differ slightly between the two campaigns, because nov26 clusters only particles with pT > 0.5 GeV and |η| < 4.7.

| Tested on | Metric | NNj (summer26) | NNjPF (nov26) | ParTj (nov26) | NNjPF − NNj | ParTj − NNj |
|---|---|---|---|---|---|---|
| Herwig | AUC | 0.6861 | 0.6866 | 0.6904 | +0.0005 ± 0.0016 | +0.0043 ± 0.0017 |
| Herwig | ε_S(1 %) | 0.0329 | 0.0481 | 0.0473 | +0.0152 ± 0.0046 | +0.0144 ± 0.0054 |
| Herwig | ε_B(50 %) | 0.2433 | 0.2437 | 0.2408 | +0.0005 ± 0.0043 | −0.0025 ± 0.0044 |
| MG5+Pythia | AUC | 0.6956 | 0.6957 | 0.6965 | +0.0001 ± 0.0005 | +0.0009 ± 0.0005 |
| MG5+Pythia | ε_S(1 %) | 0.0560 | 0.0573 | 0.0587 | +0.0013 ± 0.0017 | +0.0027 ± 0.0018 |
| MG5+Pythia | ε_B(50 %) | 0.2324 | 0.2323 | 0.2326 | −0.0001 ± 0.0014 | +0.0002 ± 0.0014 |

## 4. Transformer (ParTj) against DeepSets (NNjPF)

In each campaign the two nets scored exactly the same jets, with the same folds, so their difference is measured precisely.

![ROC curves: transformer against DeepSets](figs/roc_summary.png)

[figs/roc_summary.pdf](figs/roc_summary.pdf), made by `python3 roc_summary.py --figure-only`:
- **Line styles:** the transformer (ParTj) is dashed and DeepSets (NNjPF) solid; blue is nov26 and orange fall26. summer26 NNj, the DeepSets net without particle type, is the grey line in the VBF-preselected panels. The legends give each curve's AUC and ε_S(1 %).
- **Lower panels:** the background rejection of the transformer divided by that of DeepSets, at equal signal efficiency. Above 1 the transformer rejects more background.

| Jets | Tested on | Campaign | AUC DeepSets | AUC transformer | AUC difference | ε_S(1 %) DeepSets | ε_S(1 %) transformer | ε_S(1 %) difference |
|---|---|---|---|---|---|---|---|---|
| VBF-preselected | Herwig | nov26 | 0.6889 | 0.6926 | +0.0037 ± 0.0003 | 0.0505 | 0.0540 | +0.0035 ± 0.0016 |
| VBF-preselected | Herwig | fall26 | 0.6918 | 0.6952 | +0.0034 ± 0.0004 | 0.0514 | 0.0544 | +0.0031 ± 0.0018 |
| VBF-preselected | MG5+Pythia | nov26 | 0.6961 | 0.6969 | +0.0008 ± 0.0006 | 0.0575 | 0.0589 | +0.0014 ± 0.0014 |
| VBF-preselected | MG5+Pythia | fall26 | 0.6791 | 0.6813 | +0.0023 ± 0.0004 | 0.0484 | 0.0519 | +0.0034 ± 0.0025 |
| Gate jets | Herwig | nov26 | 0.7653 | 0.7666 | +0.0013 ± 0.0002 | 0.0998 | 0.0987 | −0.0011 ± 0.0009 |
| Gate jets | Herwig | fall26 | 0.7666 | 0.7686 | +0.0021 ± 0.0003 | 0.1006 | 0.1015 | +0.0009 ± 0.0021 |
| Gate jets | MG5+Pythia | nov26 | 0.7637 | 0.7634 | −0.0003 ± 0.0005 | 0.1028 | 0.1015 | −0.0013 ± 0.0015 |
| Gate jets | MG5+Pythia | fall26 | 0.7530 | 0.7544 | +0.0014 ± 0.0004 | 0.0908 | 0.0908 | 0.0000 ± 0.0016 |

- **Differences** are transformer minus DeepSets, from the pooled out-of-fold values. Their uncertainty is the larger of the bootstrap σ and the standard error of the five per-model differences; the second includes the training randomness. The single values have the uncertainties of sections 1 and 2.
- **Background at 50 % signal efficiency:**
  - On Herwig the transformer passes 1.4–1.5 % less background (ε_B(50 %) lower by 0.0021–0.0036, ± 0.0004–0.0013).
  - On MG5+Pythia the change ranges from −0.0011 to +0.0003. The exception is nov26 gate jets, where the transformer passes 1.7 % more (+0.0025 ± 0.0006).

---

## What the tables show

- **Over most of the ROC curve, NNjPF performs no better than NNj.** On the same events, the AUC and ε_B(50 %) of the two agree within the uncertainties (on MG5+Pythia the AUC difference is +0.0001 ± 0.0005). The particle type, the charge and nine times more training jets make no measurable difference there.
- **At ε_S(1 %) NNjPF is better on Herwig, but not on MG5+Pythia.**
  - On the Herwig test events NNjPF is ahead by 0.015 ± 0.005 (3σ); on MG5+Pythia by 0.001 ± 0.002.
  - The Herwig value rests on about 65 QCD jets above the cut.
  - A gain that shows only on the training generator could come from Herwig-specific features.
- **The transformer (ParTj) is slightly better than DeepSets (NNjPF) on Herwig, and the gain shrinks on MG5+Pythia** (section 4).
  - **On Herwig:** in both campaigns the AUC is higher by 0.001–0.004, and ε_B(50 %) is lower by 0.002–0.004 (1.5 % less background).
  - **At ε_S(1 %):** there is a gain only for the VBF-preselected jets (+0.003 ± 0.002, about 2σ), and none for all gate jets.
  - **On MG5+Pythia:** the nov26 gain disappears; for gate jets ParTj is even slightly worse at ε_B(50 %). fall26 keeps about two thirds of the AUC gain.
  - **Against summer26 NNj** on the same events, ParTj gains 0.004 ± 0.002 in AUC on Herwig and 0.0009 ± 0.0005 on MG5+Pythia.
- **fall26 against nov26:**
  - On Herwig the two campaigns are close: fall26 is higher by 0.001–0.003 in AUC for both nets (NNjPF 0.692 against 0.689 for VBF-preselected jets, 0.767 against 0.765 for gate jets).
  - On MG5+Pythia fall26 is lower for both nets (NNjPF 0.679 against 0.696, and 0.753 against 0.764). This is the transfer loss of [fall26/training.md](fall26/training.md), section 6.
  - Trained on MG5+Pythia, fall26 NNjPF matches nov26 on gate-jet ε_S(1 %) (0.102 against 0.103) and comes closer in AUC (0.760 against 0.764).
- **The preselection keeps the hard cases.** The same nets reach AUC 0.69 on VBF-preselected jets and 0.77 on all gate jets: the mjj and ΔY_jj cuts keep the QCD events that look most like VBF.

## What differs between the rows

| | summer26 NNj | nov26 NNjPF | fall26 NNjPF |
|---|---|---|---|
| Generator samples | Herwig campaign-00023 and MG5+Pythia with dipole recoil | the same files | the Fall2026 production (13.6 TeV; Herwig 7.3.0, MG5 3.5.15 + Pythia 8.317), see [fall26.md](fall26.md) |
| Particles clustered | all final-state particles (pT > 0, \|η\| < 20) | pT > 0.5 GeV, \|η\| < 4.7 | as nov26 |
| Network inputs | per constituent: signed Δη, Δφ, pT/pT_jet; per jet: \|η\|, m, pT | the same DeepSets network, plus the PF type (one-hot) and the charge of each constituent | as nov26 |
| Training jets (Herwig) | the slot-0 jet of VBF-preselected events, 70 % split: about 160 000 (16 000 QCD) | the slot-0 jet of gate events and every jet of gate == 0 events, all events, five folds: 1.47 M per model (680 000 QCD) | as nov26: 1.49 M per model, with event weights in the loss |
| Evaluation | Herwig: 15 % test split; MG5+Pythia: all events | out of fold, all events | as nov26 |

**ParTj** has the same inputs, training jets, folds and weights as NNjPF in its campaign.
- **Architecture:** attention over the constituents, with pairwise features of every constituent pair (ln kT, ln z, ln ΔR, ln m²).
- **Size:** 265 000 parameters, against 22 000 for NNjPF. Details are in [nov26/training.md](nov26/training.md), section 7.

## Notes

- **Weights in summer26 and nov26:** the Herwig weights are narrow (N_eff/N = 0.96–0.98 in these categories). Weighted and unweighted values differ by at most 0.0015.
  - summer26.md and nov26/training.md quote unweighted values: for example NNj AUC 0.6847 (Herwig) and 0.6987 (MG5+Pythia) for the slot-0 jets alone, and NNjPF gate-jet AUC 0.7643 and 0.7637.
- **Weights in fall26:** the Herwig VBF weights are broad (N_eff/N = 0.24 for the VBF-preselected jets). Unweighted values would be far off: for example the VBF-preselected NNjPF AUC would be 0.658 instead of 0.692. Only the weighted values describe the physics.
- **MG5+Pythia weights** are constant in all three campaigns.
- **Uncertainties** in sections 1–3 are statistical only: they leave out the training randomness. For nov26 and fall26, the spread of the five fold models is in [nov26/training.md](nov26/training.md) and [fall26/training.md](fall26/training.md).
- **VBF-preselected events** (VBF / QCD):

  | Generator | summer26 | nov26 | fall26 |
  |---|---|---|---|
  | Herwig | 30 788 / 3 423 (test split) | 203 331 / 21 962 | 168 269 / 23 116 |
  | MG5+Pythia | 207 827 / 30 452 | 205 548 / 29 266 | 202 719 / 19 831 |

  nov26 has 1–4 % fewer events than summer26 on the same inputs, because of its particle cuts.
