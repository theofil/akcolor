# ROC summary: NNj and NNjBC (NNj + beam connection)

The single-jet VBF taggers of the Z channel: VBF Z jets (signal) against QCD Z+jj jets (background).
- **NNj** (summer26): a DeepSets network on the jet's constituents, with the jet's |η|, m and pT.
- **NNj** (nov26): the same network, trained with five-fold cross-fitting on the nov26 friends.
- **NNjBC** (nov26): NNj plus the four beam-connection values of the jet (bc_in, bc_out, bc_met, bc_dphi: the activity between the jet and its own beam, defined in [nov26.md](nov26.md)).

**What each beam-connection value adds:** section 4. It measures each value on its own (NNj + that value) and inside NNjBC (NNjBC without that value).

Three figures of merit:
- **AUC:** the area under the ROC curve.
- **ε_S(1 %):** the signal efficiency at 1 % background efficiency (higher is better).
- **ε_B(50 %):** the background efficiency at 50 % signal efficiency (lower is better).

All values are event-weighted (kWeight). The uncertainties are statistical, from a bootstrap with 400 replicas: each event gets a Poisson multiplicity, shared by its two jets, and all nov26 nets see the same replicas. The numbers come from the scores stored in the h5 files (`python3 roc_summary.py`, 2026-11-18).

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

---

## What the tables show

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
- **The preselection keeps the hard cases.** The same nets reach AUC 0.68–0.71 on VBF-preselected jets and 0.76–0.78 on all gate jets: the mjj and ΔY_jj cuts keep the QCD events that look most like VBF.

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
- **Uncertainties** in sections 1–3 are statistical only and leave out the training randomness. The five fold models of NNj and NNjBC spread by an rms of 0.0016–0.0017 in AUC on Herwig gate jets. Section 4 includes that spread.
- **The discontinued PF inputs, for scale.** On the events of section 3, NNjPF (NNj + particle type and charge, nov26, deleted 2026-11-18) scored 0.6866 (Herwig) and 0.6957 (MG5+Pythia) in AUC, 0.004 and 0.0016 above nov26 NNj. On gate jets it was 0.002 and 0.0006 above (unweighted 0.7643 and 0.7637, against 0.7624 and 0.7631). The beam connection gains several times more, and the gain does not fade on MG5+Pythia.
- **VBF-preselected events** (VBF / QCD):

  | Generator | summer26 | nov26 |
  |---|---|---|
  | Herwig | 30 788 / 3 423 (test split) | 203 331 / 21 962 |
  | MG5+Pythia | 207 827 / 30 452 | 205 548 / 29 266 |

  nov26 has 1–4 % fewer events than summer26 on the same inputs, because of its particle cuts.
