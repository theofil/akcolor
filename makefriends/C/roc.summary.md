# ROC summary: NNj and NNjBC (NNj + beam connection)

The single-jet VBF taggers of the Z channel: VBF Z jets (signal) against QCD Z+jj jets (background).
- **NNj** (summer26): a DeepSets network on the jet's constituents, with the jet's |η|, m and pT.
- **NNj** (nov26): the same network, trained with five-fold cross-fitting on the nov26 friends.
- **NNjBC** (nov26): NNj plus the four beam-connection values of the jet (bc_in, bc_out, bc_met, bc_dphi: the activity between the jet and its own beam, defined in [nov26.md](nov26.md)).
- **Kj** and **KjBC** (nov26): NNj and NNjBC without the constituent branch. Kj sees only the jet's |η|, m and pT, so it is the kinematics-only reference.

**What each beam-connection value adds:** section 4. It measures each value on its own (NNj + that value) and inside NNjBC (NNjBC without that value).

**What the constituents and the beam connection add beyond the jet kinematics:** section 5. It compares the four nets Kj, KjBC, NNj and NNjBC at equal background efficiency. Start there for the question "how useful is NNj?".

Four figures of merit:
- **AUC:** the area under the ROC curve.
- **ε_S(1 %):** the signal efficiency at 1 % background efficiency (higher is better).
- **ε_B(50 %):** the background efficiency at 50 % signal efficiency (lower is better).
- **ε_S ratio at equal ε_B** (section 5): ε_S(net) / ε_S(reference) at the same background efficiency. Both nets then keep the same background, so the ratio is directly the gain in S/√B, in S/B and in S/(δB) at that working point.

All values are event-weighted (kWeight). The uncertainties are statistical, from a bootstrap with 400 replicas: each event gets a Poisson multiplicity, shared by its two jets, and all nov26 nets see the same replicas. The numbers come from the scores stored in the h5 files (`python3 roc_summary.py`, 2026-11-18; section 5: `python3 nov26/gain_curves.py`).

---

## 1. VBF-preselected events: shared by both campaigns

These are both jets of the events that pass the summer26 selection:
- two jets with pT > 30 GeV and |η| < 3, in opposite hemispheres (the nov26 gate);
- mjj > 400 GeV and |ΔY_jj| > 2.5.

The summer26 friends contain only these events, so this is the only category in which summer26 NNj can be compared.

| Tested on | Campaign | Net | AUC | ε_S(1 %) | ε_B(50 %) |
|---|---|---|---|---|---|
| Herwig | summer26 | NNj | 0.6891 ± 0.0036 | 0.0341 ± 0.0045 | 0.2411 ± 0.0056 |
| Herwig | nov26 | NNj | 0.6843 ± 0.0013 | 0.0474 ± 0.0022 | 0.2482 ± 0.0020 |
| Herwig | nov26 | NNjBC | 0.6990 ± 0.0013 | 0.0589 ± 0.0020 | 0.2309 ± 0.0021 |
| MG5+Pythia | summer26 | NNj | 0.6981 ± 0.0011 | 0.0570 ± 0.0020 | 0.2293 ± 0.0018 |
| MG5+Pythia | nov26 | NNj | 0.6946 ± 0.0011 | 0.0516 ± 0.0015 | 0.2326 ± 0.0018 |
| MG5+Pythia | nov26 | NNjBC | 0.7081 ± 0.0011 | 0.0590 ± 0.0020 | 0.2159 ± 0.0018 |

- **summer26 NNj on Herwig** uses only its 15 % test split (30 788 VBF and 3 423 QCD events), hence the larger uncertainties.
- **Every other entry uses all events.** Each event is scored by a model that never saw it: five-fold cross-fitting for nov26, and a net trained on the other generator for summer26 on MG5+Pythia.

## 2. All gate jets: nov26

These are both jets of gate events, without the mjj and ΔY_jj cuts. summer26 has no entry, because its friends hold only VBF-preselected events.

| Tested on | Net | AUC | ε_S(1 %) | ε_B(50 %) |
|---|---|---|---|---|
| Herwig | NNj | 0.7633 ± 0.0006 | 0.0963 ± 0.0014 | 0.1510 ± 0.0008 |
| Herwig | NNjBC | 0.7762 ± 0.0006 | 0.1034 ± 0.0015 | 0.1374 ± 0.0008 |
| MG5+Pythia | NNj | 0.7631 ± 0.0005 | 0.1017 ± 0.0012 | 0.1458 ± 0.0007 |
| MG5+Pythia | NNjBC | 0.7764 ± 0.0005 | 0.1098 ± 0.0012 | 0.1330 ± 0.0007 |

Both nets scored exactly the same jets with the same folds, so their difference is known much better than either value: NNjBC − NNj is +0.0129 ± 0.0002 (Herwig) and +0.0133 ± 0.0002 (MG5+Pythia) in AUC, and +0.0071 ± 0.0012 and +0.0082 ± 0.0014 in ε_S(1 %) (section 4).

![ROC curves: NNjBC against NNj](figs/roc_summary.png)

[figs/roc_summary.pdf](figs/roc_summary.pdf): NNj in blue, NNjBC in orange, and summer26 NNj as a dotted grey line in the VBF-preselected panels. The legends give each curve's AUC and ε_S(1 %). The lower panels show the background rejection of NNjBC divided by that of NNj at equal signal efficiency.
- **Gate jets:** NNjBC rejects 10–14 % more background between 10 % and 50 % signal efficiency, on both generators.
- **VBF-preselected jets:** the ratio rises toward high purity, to 1.17 (MG5+Pythia) and 1.32 (Herwig) at 5 % signal efficiency.

## 3. summer26 NNj against nov26 NNj and NNjBC on the same events

summer26 and nov26 read the same generator files, so their nets can be compared on identical events, matched by the boson variables.
- **Events:** the summer26 evaluation events that also pass the selection in nov26, 99 % of the VBF and 95 % of the QCD events. That is 30 418 VBF and 3 255 QCD events on Herwig, and 205 125 and 29 038 on MG5+Pythia.
- **Bootstrap:** all three nets see the same replicas, so their differences are measured more precisely than the values themselves.
- **Jets:** the jets can differ slightly between the two campaigns, because nov26 clusters only particles with pT > 0.5 GeV and |η| < 4.7.

| Tested on | Metric | NNj (summer26) | NNj (nov26) | NNjBC (nov26) | NNj (nov26) − NNj (summer26) | NNjBC − NNj (summer26) | NNjBC − NNj (nov26) |
|---|---|---|---|---|---|---|---|
| Herwig | AUC | 0.6861 | 0.6828 | 0.6971 | −0.0032 ± 0.0014 | +0.0110 ± 0.0018 | +0.0143 ± 0.0014 |
| Herwig | ε_S(1 %) | 0.0329 | 0.0408 | 0.0516 | +0.0079 ± 0.0044 | +0.0187 ± 0.0051 | +0.0108 ± 0.0044 |
| Herwig | ε_B(50 %) | 0.2433 | 0.2491 | 0.2351 | +0.0058 ± 0.0042 | −0.0081 ± 0.0044 | −0.0140 ± 0.0039 |
| MG5+Pythia | AUC | 0.6956 | 0.6941 | 0.7078 | −0.0015 ± 0.0005 | +0.0121 ± 0.0006 | +0.0137 ± 0.0005 |
| MG5+Pythia | ε_S(1 %) | 0.0560 | 0.0514 | 0.0589 | −0.0046 ± 0.0015 | +0.0029 ± 0.0020 | +0.0075 ± 0.0017 |
| MG5+Pythia | ε_B(50 %) | 0.2324 | 0.2331 | 0.2163 | +0.0008 ± 0.0013 | −0.0161 ± 0.0015 | −0.0168 ± 0.0012 |

- **The two NNj trainings are close.** On these events, nov26 NNj is lower by 0.003 (Herwig) and 0.0015 (MG5+Pythia) in AUC. It is trained on all gate and gate == 0 jets, while summer26 NNj is trained on exactly the VBF-preselected population. At ε_S(1 %) the two differ in opposite directions on the two generators, with summer26 resting on only 3 255 Herwig QCD events.
- **NNjBC beats both,** by 0.011–0.014 in AUC on both generators.

## 4. What each beam-connection value adds

Every value is measured twice, on the same jets with the same folds:
- **Add-one:** (NNj + that value) − NNj. This is what the value carries on its own, including what it shares with the other values.
- **Drop-one:** NNjBC − (NNjBC without that value). This is what only that value carries, given the other three.

The values are correlated (bc_in + bc_out is the total pT in the region, and bc_met ≤ that total), so the two numbers bracket a value's contribution. Each difference has the uncertainty max(bootstrap σ, standard error of the five per-fold-model differences); the second includes the training randomness.

![Beam-connection ablation](figs/nov26/bc_ablation.png)

[figs/nov26/bc_ablation.pdf](figs/nov26/bc_ablation.pdf): black is NNjBC − NNj, blue the add-one gains, orange the drop-one losses. Points above zero mean the value helps.

**Gate jets, Δ AUC**

| Value | Add-one, Herwig | Drop-one, Herwig | Add-one, MG5+Pythia | Drop-one, MG5+Pythia |
|---|---|---|---|---|
| all four (NNjBC − NNj) | +0.0129 ± 0.0002 | | +0.0133 ± 0.0002 | |
| bc_in | +0.0001 ± 0.0001 | −0.0002 ± 0.0001 | +0.0005 ± 0.0003 | −0.0005 ± 0.0001 |
| bc_out | +0.0055 ± 0.0003 | +0.0012 ± 0.0002 | 0.0000 ± 0.0002 | +0.0010 ± 0.0001 |
| bc_met | +0.0056 ± 0.0002 | +0.0028 ± 0.0001 | +0.0051 ± 0.0001 | +0.0067 ± 0.0002 |
| bc_dphi | +0.0031 ± 0.0003 | +0.0025 ± 0.0002 | +0.0026 ± 0.0003 | +0.0025 ± 0.0003 |

**Gate jets, Δ ε_S(1 %)**

| Value | Add-one, Herwig | Drop-one, Herwig | Add-one, MG5+Pythia | Drop-one, MG5+Pythia |
|---|---|---|---|---|
| all four (NNjBC − NNj) | +0.0071 ± 0.0012 | | +0.0082 ± 0.0014 | |
| bc_in | +0.0010 ± 0.0008 | −0.0001 ± 0.0007 | +0.0004 ± 0.0010 | −0.0008 ± 0.0006 |
| bc_out | +0.0026 ± 0.0010 | −0.0007 ± 0.0009 | +0.0018 ± 0.0009 | −0.0004 ± 0.0013 |
| bc_met | −0.0013 ± 0.0009 | −0.0016 ± 0.0009 | +0.0013 ± 0.0009 | +0.0004 ± 0.0013 |
| bc_dphi | +0.0064 ± 0.0012 | +0.0041 ± 0.0009 | +0.0040 ± 0.0010 | +0.0029 ± 0.0010 |

**VBF-preselected jets, Δ AUC**

| Value | Add-one, Herwig | Drop-one, Herwig | Add-one, MG5+Pythia | Drop-one, MG5+Pythia |
|---|---|---|---|---|
| all four (NNjBC − NNj) | +0.0147 ± 0.0008 | | +0.0136 ± 0.0004 | |
| bc_in | +0.0015 ± 0.0002 | 0.0000 ± 0.0002 | +0.0014 ± 0.0002 | −0.0004 ± 0.0001 |
| bc_out | +0.0049 ± 0.0006 | +0.0015 ± 0.0003 | +0.0011 ± 0.0004 | 0.0000 ± 0.0003 |
| bc_met | +0.0012 ± 0.0006 | +0.0019 ± 0.0003 | +0.0004 ± 0.0003 | +0.0036 ± 0.0003 |
| bc_dphi | +0.0071 ± 0.0004 | +0.0047 ± 0.0003 | +0.0076 ± 0.0004 | +0.0044 ± 0.0003 |

**VBF-preselected jets, Δ ε_S(1 %)**

| Value | Add-one, Herwig | Drop-one, Herwig | Add-one, MG5+Pythia | Drop-one, MG5+Pythia |
|---|---|---|---|---|
| all four (NNjBC − NNj) | +0.0115 ± 0.0020 | | +0.0074 ± 0.0022 | |
| bc_in | +0.0022 ± 0.0016 | +0.0025 ± 0.0015 | +0.0037 ± 0.0013 | −0.0009 ± 0.0014 |
| bc_out | +0.0029 ± 0.0016 | +0.0033 ± 0.0016 | +0.0020 ± 0.0020 | +0.0015 ± 0.0012 |
| bc_met | −0.0006 ± 0.0017 | +0.0018 ± 0.0014 | −0.0009 ± 0.0014 | +0.0008 ± 0.0011 |
| bc_dphi | +0.0080 ± 0.0020 | +0.0040 ± 0.0014 | +0.0065 ± 0.0015 | +0.0013 ± 0.0016 |

ε_B(50 %) follows the AUC; its values are in the `roc_summary.py` printout.

## 5. What the constituents and the beam connection add beyond the jet kinematics

Two problems make the AUC a poor measure of how useful a tagger is:
- **It depends on the population, not only on the net.** NNj has AUC 0.76 on gate jets and 0.68 on VBF-preselected jets (sections 1 and 2), and most of that difference is the jet kinematics (below).
- **It averages over every working point.** Most of the area comes from loose cuts that an analysis would not use.

This section measures each source of information against the kinematics-only net, on the same jets:

| Net | Jet \|η\|, m, pT | Constituents | Beam connection |
|---|---|---|---|
| Kj | ✓ | | |
| KjBC | ✓ | | ✓ |
| NNj | ✓ | ✓ | |
| NNjBC | ✓ | ✓ | ✓ |

All four share the folds, the training jets and the recipe ([nov26/training.md](nov26/training.md)); Kj and KjBC are NNj and NNjBC without the constituent branch.

**The figure of merit is the signal efficiency ratio at equal background efficiency**, R(ε_B) = ε_S(net) / ε_S(reference).
- At equal ε_B both nets keep the same background, so R is the ratio of S/√B, of S/B and of S/(δB) at that working point (for S ≪ B; S/B ≈ 3 % after the VBF preselection).
- It therefore holds whether the analysis is limited by statistics or by a background systematic. At 300 fb⁻¹ the preselected background is about 2.3 M events (summer26.md, Z_only_table), so any background uncertainty above 1/√B ≈ 0.07 % already dominates its statistical error.
- R² is the luminosity equivalent at that working point.

![Signal efficiency ratios at equal background efficiency](figs/nov26/gain_curves.png)

[figs/nov26/gain_curves.pdf](figs/nov26/gain_curves.pdf):
- **Upper panels:** ε_S/√ε_B against ε_B, the factor by which a cut on the score multiplies S/√B of the category (the "performance" of the summer26 table; 1 = no cut).
- **Lower panels:** the ratios to Kj, with ±1σ bootstrap bands. Herwig is solid, MG5+Pythia dashed.
- **Categories:** gate jets; VBF-preselected jets; and VBF-preselected events with the two jet scores combined as logit(s₀) + logit(s₁), the same combination for every net.
- The curves stop where 100 raw background entries remain: at ε_B = 0.03 % (gate jets), 0.23 % / 0.17 % (VBF-preselected jets) and 0.46 % / 0.34 % (VBF-preselected events), Herwig / MG5+Pythia.

**VBF-preselected jets, R at equal ε_B**

| Ratio | Measures | Herwig, ε_B = 10 % | 1 % | 0.3 % | MG5+Pythia, ε_B = 10 % | 1 % | 0.3 % |
|---|---|---|---|---|---|---|---|
| KjBC / Kj | beam connection | 1.114 ± 0.013 | 1.252 ± 0.064 | 1.409 ± 0.104 | 1.094 ± 0.009 | 1.141 ± 0.033 | 1.246 ± 0.081 |
| NNj / Kj | constituents | 1.103 ± 0.012 | 1.170 ± 0.062 | 1.237 ± 0.083 | 1.088 ± 0.007 | 1.142 ± 0.031 | 1.148 ± 0.079 |
| NNjBC / Kj | both | 1.206 ± 0.017 | 1.453 ± 0.060 | 1.583 ± 0.148 | 1.154 ± 0.010 | 1.306 ± 0.045 | 1.457 ± 0.106 |
| NNjBC / NNj | beam connection on top of the constituents | 1.094 ± 0.009 | 1.242 ± 0.056 | 1.279 ± 0.117 | 1.061 ± 0.007 | 1.143 ± 0.040 | 1.269 ± 0.085 |
| NNjBC / KjBC | constituents on top of the beam connection | 1.082 ± 0.010 | 1.160 ± 0.035 | 1.123 ± 0.135 | 1.055 ± 0.006 | 1.144 ± 0.042 | 1.169 ± 0.066 |

**VBF-preselected events, both jets, R at equal ε_B** (0.3 % is below the 100-entry floor)

| Ratio | Herwig, ε_B = 10 % | 1 % | MG5+Pythia, ε_B = 10 % | 1 % |
|---|---|---|---|---|
| KjBC / Kj | 1.128 ± 0.012 | 1.359 ± 0.074 | 1.093 ± 0.010 | 1.200 ± 0.048 |
| NNj / Kj | 1.134 ± 0.014 | 1.257 ± 0.076 | 1.116 ± 0.014 | 1.337 ± 0.058 |
| NNjBC / Kj | 1.255 ± 0.015 | 1.702 ± 0.086 | 1.189 ± 0.012 | 1.506 ± 0.060 |
| NNjBC / NNj | 1.106 ± 0.009 | 1.354 ± 0.083 | 1.065 ± 0.008 | 1.126 ± 0.056 |
| NNjBC / KjBC | 1.113 ± 0.011 | 1.252 ± 0.048 | 1.088 ± 0.007 | 1.254 ± 0.074 |

**Gate jets, R at equal ε_B**

| Ratio | Herwig, ε_B = 10 % | 1 % | 0.3 % | MG5+Pythia, ε_B = 10 % | 1 % | 0.3 % |
|---|---|---|---|---|---|---|
| KjBC / Kj | 1.071 ± 0.003 | 1.109 ± 0.015 | 1.168 ± 0.033 | 1.058 ± 0.004 | 1.123 ± 0.012 | 1.149 ± 0.026 |
| NNj / Kj | 1.045 ± 0.003 | 1.128 ± 0.013 | 1.198 ± 0.045 | 1.027 ± 0.002 | 1.094 ± 0.010 | 1.148 ± 0.021 |
| NNjBC / Kj | 1.112 ± 0.005 | 1.210 ± 0.016 | 1.304 ± 0.040 | 1.087 ± 0.004 | 1.182 ± 0.013 | 1.271 ± 0.028 |
| NNjBC / NNj | 1.064 ± 0.005 | 1.073 ± 0.013 | 1.088 ± 0.027 | 1.058 ± 0.003 | 1.081 ± 0.014 | 1.106 ± 0.023 |
| NNjBC / KjBC | 1.038 ± 0.003 | 1.091 ± 0.015 | 1.116 ± 0.029 | 1.027 ± 0.002 | 1.053 ± 0.011 | 1.106 ± 0.021 |

Each uncertainty is the larger of the paired bootstrap σ (200 replicas, shared by all nets) and the standard error of the five per-fold-model ratios; the second includes the training randomness.

**Maximum of ε_S/√ε_B, with the ε_B where it sits** (the best single cut for a statistics-limited analysis)

| Category | Generator | Kj | KjBC | NNj | NNjBC |
|---|---|---|---|---|---|
| VBF-preselected jets | Herwig | 1.041 (69 %) | 1.061 (56 %) | 1.068 (57 %) | 1.091 (48 %) |
| | MG5+Pythia | 1.056 (54 %) | 1.084 (44 %) | 1.086 (45 %) | 1.112 (40 %) |
| VBF-preselected events, both jets | Herwig | 1.137 (42 %) | 1.216 (24 %) | 1.221 (25 %) | 1.310 (11 %) |
| | MG5+Pythia | 1.216 (18 %) | 1.309 (14 %) | 1.321 (12 %) | 1.407 (11 %) |
| Gate jets | Herwig | 1.249 (19 %) | 1.308 (18 %) | 1.291 (17 %) | 1.349 (14 %) |
| | MG5+Pythia | 1.285 (15 %) | 1.343 (12 %) | 1.311 (15 %) | 1.373 (11 %) |

**Where the AUC differences between populations come from** (AUC, Herwig / MG5+Pythia)

| Input | Gate jets | Gate, \|ΔY_jj\| > 2.5 | Gate, mjj > 400 GeV | VBF-preselected jets |
|---|---|---|---|---|
| pT alone | 0.7211 / 0.7261 | 0.7348 / 0.7355 | 0.5873 / 0.6260 | 0.6367 / 0.6613 |
| \|η\| alone | 0.6275 / 0.6270 | 0.5574 / 0.5574 | 0.5362 / 0.5137 | 0.5119 / 0.5048 |
| m alone | 0.5918 / 0.6081 | 0.6029 / 0.6172 | 0.4962 / 0.5435 | 0.5268 / 0.5635 |
| Kj | 0.7538 / 0.7556 | 0.7525 / 0.7508 | 0.6397 / 0.6569 | 0.6643 / 0.6773 |
| KjBC | 0.7675 / 0.7697 | 0.7621 / 0.7613 | 0.6635 / 0.6776 | 0.6821 / 0.6944 |
| NNj | 0.7633 / 0.7631 | 0.7629 / 0.7596 | 0.6624 / 0.6773 | 0.6843 / 0.6946 |
| NNjBC | 0.7762 / 0.7764 | 0.7714 / 0.7684 | 0.6823 / 0.6938 | 0.6990 / 0.7081 |

For pT, |η| and m alone, larger values count as VBF-like in every column.

**The same nets at fixed jet kinematics** (AUC within 15 × 20 bins of \|η\| and log pT, averaged with signal weights; the bins with ≥ 50 raw jets of each class hold at least 94.8 % of the signal)

| Input | Gate jets | Gate, \|ΔY_jj\| > 2.5 | Gate, mjj > 400 GeV | VBF-preselected jets |
|---|---|---|---|---|
| m alone (lighter = VBF-like) | 0.5815 / 0.5616 | 0.5776 / 0.5567 | 0.5856 / 0.5570 | 0.5819 / 0.5550 |
| Kj | 0.5887 / 0.5696 | 0.5862 / 0.5660 | 0.5889 / 0.5624 | 0.5867 / 0.5617 |
| KjBC | 0.6201 / 0.6105 | 0.6216 / 0.6113 | 0.6230 / 0.6075 | 0.6226 / 0.6095 |
| NNj | 0.6231 / 0.6062 | 0.6235 / 0.6059 | 0.6291 / 0.6082 | 0.6279 / 0.6071 |
| NNjBC | 0.6466 / 0.6347 | 0.6487 / 0.6350 | 0.6527 / 0.6366 | 0.6522 / 0.6365 |

|η| and pT alone give 0.50–0.53 here, which is what remains of them inside a bin.

**What section 5 shows**
- **A single-jet tagger's AUC is mostly the jet's kinematics.**
  - On gate jets, the jet pT alone reaches AUC 0.72, and Kj (|η|, m, pT) 0.754. All the constituents add 0.009 on top (NNj 0.763), and NNjBC 0.022 (0.776).
  - The mjj > 400 GeV cut removes most of the pT separation (pT alone 0.72 → 0.59 on Herwig, 0.73 → 0.63 on MG5+Pythia) and costs every net 0.08–0.11 in AUC. The |ΔY_jj| > 2.5 cut alone costs at most 0.008.
  - That is the difference between the gate-jet and the VBF-preselected AUCs of sections 1 and 2: a property of the population, not of the tagger.
- **At fixed jet kinematics the nets no longer depend on the population.**
  - Within (|η|, pT) bins NNj has AUC 0.623–0.629 (Herwig) and 0.606–0.608 (MG5+Pythia) in all four populations, and NNjBC 0.647–0.653 and 0.635–0.637. The raw AUCs of the same nets differ by 0.08–0.10.
  - What Kj keeps at fixed (|η|, pT) is the jet mass: at fixed pT, lighter jets are VBF-like (m alone 0.578–0.586 on Herwig, Kj 0.586–0.589).
- **The beam connection is worth as much as the constituents, and the two add almost independently.**
  - At fixed kinematics: KjBC 0.620 against NNj 0.623 on Herwig, and 0.611 against 0.606 on MG5+Pythia (gate jets).
  - At ε_B = 1 % on VBF-preselected jets: KjBC / Kj = 1.25 and NNj / Kj = 1.17 on Herwig, both 1.14 on MG5+Pythia.
  - Together, NNjBC / Kj = 1.45 (Herwig) and 1.31 (MG5+Pythia), close to the products 1.46 and 1.30.
- **The gain grows toward high purity.**
  - NNjBC / Kj on VBF-preselected jets: 1.21 at ε_B = 10 %, 1.45 at 1 % and 1.58 at 0.3 % (Herwig); 1.15, 1.31 and 1.46 (MG5+Pythia).
  - With both jets of the event: 1.255 at 10 % and 1.702 at 1 % (Herwig); 1.189 and 1.506 (MG5+Pythia).
  - As a luminosity equivalent, R²: NNjBC over Kj at ε_B = 1 % on VBF-preselected jets is worth 2.1 × (Herwig) and 1.7 × (MG5+Pythia) the luminosity.
- **The beam connection on top of NNj (NNjBC / NNj), Herwig / MG5+Pythia:**
  - VBF-preselected jets: 1.09 / 1.06 at ε_B = 10 %, and 1.24 / 1.14 at 1 %.
  - VBF-preselected events, both jets: 1.11 / 1.07 at 10 %, and 1.35 / 1.13 at 1 %.
  - The ΔAUC of +0.013–0.015 (section 4) hides most of this. The area is dominated by loose cuts, where the gain is smallest.
- **A statistics-limited cut gains little from any single-jet net.**
  - On VBF-preselected jets the best ε_S/√ε_B sits at loose cuts (ε_B = 40–69 %): 1.04 (Kj) to 1.09 (NNjBC) on Herwig, 1.06 to 1.11 on MG5+Pythia.
  - With both jets of the event it is 1.14 to 1.31 (Herwig) and 1.22 to 1.41 (MG5+Pythia), at ε_B = 11–42 %.
  - The nets pay off where the background systematic limits the analysis: there the ratios at high purity apply directly.
- **On MG5+Pythia** (the nets are trained on Herwig):
  - At fixed kinematics (gate jets), what the constituents add over Kj transfers (AUC +0.034 on Herwig, +0.036 on MG5+Pythia), and so does the beam connection (+0.031 and +0.041). The loss on MG5+Pythia is in the jet-mass part: Kj 0.589 → 0.570.
  - At high purity the beam connection gains less on MG5+Pythia: KjBC / Kj at ε_B = 1 % on VBF-preselected jets is 1.25 on Herwig and 1.14 on MG5+Pythia.

---

## What the tables show

- **At equal background efficiency the beam connection is a sizeable gain, which the AUC understates** (section 5).
  - On VBF-preselected jets NNjBC keeps 1.09 (Herwig) and 1.06 (MG5+Pythia) times NNj's signal at ε_B = 10 %, and 1.24 and 1.14 times at ε_B = 1 %. Using both jets of the event, the factor at 1 % is 1.35 and 1.13. At equal ε_B this is the gain in significance.
  - Beyond the jet kinematics, the beam connection is worth as much as all the constituents, and the two add almost independently: over the kinematics-only Kj, NNjBC keeps 1.45 (Herwig) and 1.31 (MG5+Pythia) times the signal at ε_B = 1 %.
  - The AUC gain of +0.013–0.015 is dominated by loose cuts, where the gain is smallest, and the AUC itself depends on the population (see the last bullet).
- **The beam connection improves NNj on both generators.**
  - On the same jets, NNjBC − NNj is +0.013 in AUC for gate jets and +0.014–0.015 for VBF-preselected jets.
  - ε_S(1 %) rises by 0.007–0.012. ε_B(50 %) falls by 0.013–0.017: 9 % less background for gate jets, 7 % less for VBF-preselected jets.
  - Unlike the particle-type inputs of the discontinued NNjPF, the gain is the same on MG5+Pythia as on Herwig, the generator the nets were trained on.
- **bc_dphi is the most robust value.**
  - It is the only value that raises the AUC in every category, on both generators, both alone and inside NNjBC: AUC +0.003 (gate jets) and +0.007 (VBF-preselected) on its own, and still +0.0025 and +0.0045 when the other three are present.
  - It carries most of the gain at high purity: ε_S(1 %) +0.004 to +0.008 on its own, and removing it costs 0.003–0.004 for gate jets.
  - This matches the distributions in [nov26.md](nov26.md): for VBF the beam-side activity lines up with the jet in φ, while for QCD its direction is nearly flat.
- **bc_met carries the largest unique AUC gain for gate jets, but none at high purity.**
  - Removing it from NNjBC costs 0.0028 (Herwig) and 0.0067 (MG5+Pythia) in AUC, the largest drop-one loss.
  - At ε_S(1 %) it adds nothing measurable, alone or inside NNjBC (all within 1.5σ of zero).
  - In VBF-preselected jets it adds little on its own (≤ 0.0012), but removing it from NNjBC costs 0.002–0.004. It works through its combination with the other values: bc_met against bc_in + bc_out tells how collimated the activity in the region is.
- **bc_out helps on its own only on Herwig.**
  - Alone it adds +0.0055 (gate jets) and +0.0049 (VBF-preselected) in AUC on Herwig, but nothing (gate jets) or +0.0011 (VBF-preselected) on MG5+Pythia.
  - Inside NNjBC it adds only about 0.001. Its standalone gain is mostly shared with bc_met, and the part that does not transfer is consistent with the larger beam-side activity of QCD Z+jj in MG5+Pythia ([nov26.md](nov26.md), bc distributions).
- **bc_in adds nothing.**
  - Its effect is ≤ 0.0005 in AUC for gate jets, either way. Alone it adds +0.0015 for VBF-preselected jets, but nothing once the other values are present.
  - NNjBC without bc_in (`NNjBC_noin`) is as good as NNjBC, and on MG5+Pythia marginally better (by 0.0005 in AUC). It is the leaner choice.
- **The contributions are not additive.**
  - The four drop-one losses add up to 0.006 of the 0.013 total for Herwig gate jets, and the four add-one gains to 0.014. About half of the gain is shared between values, mostly between bc_out and bc_met, which both measure the activity in the region.
- **The AUC depends on the population, not only on the net.** The same nets reach AUC 0.68–0.71 on VBF-preselected jets and 0.76–0.78 on all gate jets.
  - The difference is jet kinematics: the mjj > 400 GeV cut removes most of the separation in jet pT (pT alone: AUC 0.72 → 0.59 on Herwig), while the |ΔY_jj| cut alone barely matters.
  - At fixed jet (|η|, pT) every net has the same AUC in all four populations of section 5. Compare nets on the same jets, and measure usefulness where the tagger would be used.

## What differs between the rows

| | summer26 NNj | nov26 NNj | nov26 NNjBC |
|---|---|---|---|
| Generator samples | Herwig campaign-00023 and MG5+Pythia with dipole recoil | the same files | the same files |
| Particles clustered | all final-state particles (pT > 0, \|η\| < 20) | pT > 0.5 GeV, \|η\| < 4.7 | as nov26 NNj |
| Network inputs | per constituent: signed Δη, Δφ, pT/pT_jet; per jet: \|η\|, m, pT | the same | the same, plus log(1 + bc_in), log(1 + bc_out), log(1 + bc_met) and cos, sin of bc_dphi |
| Training jets (Herwig) | the slot-0 jet of VBF-preselected events, 70 % split: about 160 000 (16 000 QCD) | the slot-0 jet of gate events and every jet of gate == 0 events, all events, five folds: 1.47 M per model (680 000 QCD) | as nov26 NNj |
| Evaluation | Herwig: 15 % test split; MG5+Pythia: all events | out of fold, all events | as nov26 NNj |

The add-one and drop-one nets of section 4 are nov26 NNj with the stated bc values; they share the training jets, folds and recipe ([nov26/training.md](nov26/training.md), section 7).

## Notes

- **Weights:** the nov26 and summer26 Herwig weights are narrow (N_eff/N = 0.96–0.98 in these categories), and the MG5+Pythia weights are constant. Weighted and unweighted values differ by at most 0.0015. [nov26/training.md](nov26/training.md) quotes unweighted values, for example NNj and NNjBC gate-jet AUC 0.7624 and 0.7754 on Herwig.
- **Section 5** uses 200 bootstrap replicas (`nov26/gain_curves.py`, 2026-11-18) and quotes the larger of the bootstrap σ and the standard error of the five per-fold-model values, like section 4. Its AUC and ε_S(1 %) values reproduce sections 1 and 2 exactly.
- **Uncertainties** in sections 1–3 are statistical only and leave out the training randomness. The five fold models of NNj and NNjBC spread by an rms of 0.0016–0.0017 in AUC on Herwig gate jets. Section 4 includes that spread.
- **The discontinued PF inputs, for scale.** On the events of section 3, NNjPF (NNj + particle type and charge, nov26, deleted 2026-11-18) scored 0.6866 (Herwig) and 0.6957 (MG5+Pythia) in AUC, 0.004 and 0.0016 above nov26 NNj. On gate jets it was 0.002 and 0.0006 above (unweighted 0.7643 and 0.7637, against 0.7624 and 0.7631). The beam connection gains several times more, and the gain does not fade on MG5+Pythia.
- **VBF-preselected events** (VBF / QCD):

  | Generator | summer26 | nov26 |
  |---|---|---|
  | Herwig | 30 788 / 3 423 (test split) | 203 331 / 21 962 |
  | MG5+Pythia | 207 827 / 30 452 | 205 548 / 29 266 |

  nov26 has 1–4 % fewer events than summer26 on the same inputs, because of its particle cuts.
