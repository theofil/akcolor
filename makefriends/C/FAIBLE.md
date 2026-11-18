Beyond DeepSets for the akcolor VBF-vs-QCD taggers: answer and experiment plan
Status: literature research finished (started 2026-11-15). Nothing in the akcolor code has been changed.

Context
Question: the NN taggers use DeepSets, and an earlier Transformer gave the same ROC curves. Before concluding that DeepSets is as far as it goes, is there a more modern NN technology, or anything else in the model or the training, worth trying?

I took "the NN model" to be the akcolor one (not the run3/giorgos dimuon analysis):

Code: /afs/cern.ch/work/t/theofil/akcolor/makefriends/C/summer26/, nets NNkin, NNj, NNjB, NNjj, NNjjB, NNjjBj, their _H/_W variants, and the no-constituent MLPs Kj, KjB, KjjBj.
Notes: makefriends/C/summer26.md, studies/README.md, studies/O_jet.md.
What the present setup looks like (read from the code and notes):

VBF V/H+jj against QCD V/H+jj at particle level; trained on Herwig, checked on MG5+Pythia. Preselection mjj > 400 GeV, |dY_jj| > 2.5.
Herwig QCD after preselection: 22.6k to 23.3k events per channel (VBF: 147k to 205k). With the 70/15/15 split that is about 16k background events to train on and 3.4k to test on.
Model: phi 3-64-64 on (d_eta, d_phi, pT/pT_jet) per constituent, masked sum per jet, rho 128-64-1 with jet and boson scalars. Adam 1e-3, one seed, one split.
Inputs: constituents of the two or three leading R = 0.4 jets (pT > 30 GeV, |eta| < 3). Nothing from the rest of the event.
Your own results: the largest gain came from adding the third jet (Z channel S/sqrtB 76.8 to 108.9). Constituents add little on top of the jet four-vectors (KjjBj 100.4 against NNjjBj 108.9, on 40 and 34 raw background events). Constituents lower the generator resilience.
The earlier Transformer (NNtrans, git commit 8244e23) was a plain nn.TransformerEncoder: 2 layers, width 64, 4 heads, mean pooling, no pairwise features, one jet, Adam 1e-3.
Answer
On the same inputs and the same statistics, no newer architecture is likely to beat DeepSets by an amount you can measure today. The closest public benchmark to your per-jet problem is quark/gluon tagging: there, seven years of architectures moved the AUC from 0.9005 (DeepSets) to 0.916 (best foundation model), with 1.6M training jets, where you have about 16k background events. Your test split resolves the tight working point to about 9%.
The earlier Transformer result is what the literature predicts for a plain encoder. It does not show a mistake, and it does not show that there is no headroom. What makes the Particle Transformer better than DeepSets-class models is pairwise features plus very large training sets or pre-training.
The remaining information is more likely in what the network is given than in the network: the particles outside the tag-jet cones, particle type, and inputs written in a frame that respects the symmetries of the event.
Order of work: make differences measurable; improve the inputs and the statistics with DeepSets unchanged; then run one or two modern architectures as a cross-check. A statement "DeepSets is sufficient" is defensible once a pre-trained or symmetry-aware model agrees with it within errors, worded as "for these inputs and this sample size".
Evidence (primary sources read for this plan)
Point	Numbers	Source
Quark/gluon benchmark, R = 0.4, 1.6M training jets	AUC and gluon rejection at 50% quark efficiency: DeepSets (PFN-Ex) 0.9005, 34.7; ParticleNet 0.9116, 39.8; ParT 0.9121, 41.3; ParT fine-tuned 0.9151, 42.4; OmniLearn 0.9159, 43.7; OmniLearned-m 0.9162, 43.2	arXiv:2202.03772 Tab. 6; arXiv:2510.24066 Tab. 3
Particle type is worth as much as an architecture change	DeepSets without particle type 0.8911; with charge 0.8924; with realistic particle type 0.9005; with full particle type 0.9052. Latent size: "good performance achieved by l = 16"	arXiv:1810.05165 Tab. 3, Sec. 4.2
Pairwise features make ParT work	JetClass accuracy 0.861 with them, 0.849 without; rejections 20 to 30% lower without	arXiv:2202.03772 Tab. 1
Transformers need data	JetClass with 2M, 10M, 100M jets: ParT 0.836, 0.850, 0.861. Top tagging with 1.2M jets: ParT from scratch has the same AUC as ParticleNet (0.9858)	arXiv:2202.03772 Tab. 3, 5
Small samples favour symmetry-aware networks	Top tagging with about 6k, 12k, 60k jets: LorentzNet AUC 0.9793, 0.9812, 0.9839; ParticleNet 0.9687, 0.9734, 0.9807	arXiv:2201.08187 Tab. 3
L-GATr data efficiency	L-GATr on 10M JetClass jets (accuracy 0.859) is close to ParT on 100M (0.861)	arXiv:2411.00446 Tab. 5
LLoCa	Accuracy 0.864 against ParT 0.861 and L-GATr 0.866, at about a fifth of L-GATr's training time. For small datasets "data augmentation outperforms equivariance" (amplitude regression)	arXiv:2505.20280 Tab. 1, Fig. 4
Foundation models with little data	OmniLearn fine-tuned on 10% of the ATLAS top dataset matches the best earlier model on 100%. Two-jet events: "a DeepSets to mix the two jets". Pre-training uses R = 0.8 and 1.0 jets; R = 0.4 only downstream	arXiv:2404.16091; arXiv:2510.24066
VBF against gluon fusion: the information is in the whole event	AUC: kinematics BDT 0.820; plus jet-CNN of the two leading jets 0.870; BDT with third jet 0.905; plain self-attention on 100 particles 0.900; full-event CNN 0.940. About 200k training events. Jet images add nothing on top of the event image. Strong dependence on the Pythia recoil setting	arXiv:2209.05518 Tab. 6, App. B, C
VBF H to invisible with full-event calorimeter images	Limit on the invisible branching ratio improves "by a factor of three" over the shape analysis	arXiv:2008.05434 (abstract)
CMS practice for VBF	Network inputs include the quark-gluon likelihood of both tag jets and soft track-jet activity outside the tag jets	arXiv:2009.04363
Performance against resilience	Complex models reach higher raw performance "at the cost of resiliency"; IRC-safe DeepSets (EFN) and angularities are the most robust	arXiv:2509.19431
Generator gap exceeds architecture gap	Quark/gluon AUC, Pythia-trained: 0.900 on Pythia, 0.804 on Herwig. Multiplicity is the leading feature and the main generator difference	arXiv:2212.10493 Fig. 5
What quark/gluon networks learn	Leading latent directions: multiplicity and particle type, radial profile, fragmentation; seven observables reach AUC 0.874 to 0.884 against 0.902 for the network	arXiv:2507.21214
Distance to the optimum	Best tagger within a factor 2 in rejection of the estimated optimum for top tagging; for W, Z, H to gg the gap is "nearly closed in some cases"	arXiv:2511.15779; arXiv:2607.29508
Lund-plane inputs	Adding Lund tokens to ParT: 12% more rejection for H to bb, 7% for top, none for H to cc or H to 4q. A kt cut trades performance for resilience	arXiv:2605.26821 Tab. 1; arXiv:2112.09140
Colour-singlet optimal observable	Likelihood ratio for one soft emission depends on the emission angles to both colour partners: (theta_ak^2 + theta_bk^2) / theta_ab^2	arXiv:2006.10480 Eq. 13
Not found: a published constituent-level tagger for VBF aimed at colour flow, and any 2026 architecture that clearly beats L-GATr or ParT on JetClass.

What to try, in order
Step 0. Make differences measurable (first; small code change)
Five-fold cross-fitting in place of the single split, so every background event gets a score from a model that never saw it. Raw background at the tight point goes from 34 to about 230 events; the error on eff_S/sqrt(eff_B) from about 9% to about 3%; the AUC error from about 0.0024 to about 0.0009.
Three seeds per configuration; average their scores.
Choose the working point on validation folds, or quote eff_S at fixed eff_B of 10% and 1%. Today sb_estimate.py takes the maximum of S/sqrtB on the test ROC itself, which biases the quoted maximum upwards when few background events remain.
Paired bootstrap on common events for the difference between two nets.
Learning curve: train on 12.5, 25, 50, 100% of the background. Still rising at 100% means limited by statistics; flat means limited by what the inputs contain.
Check one ROC with sample_weight=kWeight: training and ROC are unweighted, while the Herwig VBF weights have a relative spread of 20% (Z) and 44% (H) after selection.
Step 1. Keep DeepSets, change how the inputs are presented (dataset.py, model.py)
Remove absolute phi: measure all azimuths from the boson and feed cos and sin of the differences. Today raw phi of each jet and of the boson is standardised and the network has to learn rotation invariance from 16k events.
Reflect eta so that a fixed tag jet is at positive eta; order tag-jet slots by pT (leadJetIndex) instead of alternating them, or process both with a shared block.
Orient constituent coordinates by the event: +d_eta away from the other tag jet, +d_phi towards it. Multi-jet nets now pass raw d_eta to phi while the jet's signed eta only reaches rho.
Give phi the scale and the context: log z, log pT_i, dR_i, plus log pT_jet, |eta_jet| and a jet-role flag appended to every constituent.
Use log pT and log m for jet and boson scalars; add the third-jet centrality.
Train W, Z, H together with a channel flag: three times the background.
Recipe: latent 128, sum and mean pooling, AdamW with weight decay.
Step 2. More background from existing Monte Carlo
Rerun makefriends with a loose dijet selection into a separate directory (the notes quote 15 to 26% QCD efficiency for the old loose selection against 2.3% now). Pre-train on the loose sample, fine-tune on the preselected one.
Step 3. Particle type and charge of the constituents
The PDG code is already read in makefriends.cpp (objects_arr[4]). Store a detector-realistic class (charged hadron, neutral hadron, photon, electron, muon) and the charge per constituent.
Step 4. The rest of the event
Store particles outside the tag jets (for example pT > 1 GeV, |eta| < 5, up to about
and all jets down to about 10 GeV. The input tree holds every final-state particle, so no new generation is needed.
Event-level DeepSets over them, with coordinates relative to the tag jets.
Judge on performance and on Herwig-to-Pythia resilience together: this region is the most sensitive to underlying-event modelling.
Step 5. Modern architectures as cross-checks (after Steps 0 to 2)
Fine-tuned foundation model as per-jet encoder: OmniLearned-s (3M weights, needs only d_eta, d_phi, log pT, log E), followed by the present rho.
Lorentz-equivariant event-level transformer: LLoCa-Transformer or L-GATr-slim (pip install lloca, pip install lgatr), with beam reference.
ParT with pairwise features, from the JetClass kinematics-only weights.
Step 6. More QCD Monte Carlo
Generator-level mjj and rapidity cuts, looser than the analysis cuts, to reach about 200k preselected QCD events per channel. Needed for any claim at eff_B near 1%.
For resilience, separately from performance
IRC-safe DeepSets (EFN: sum of z_i times a function of the angles) next to the present one; the difference measures how much of the gain is non-perturbative.
Lund-plane declusterings with a kt cut (the CambridgeAachen/ code already reclusters).
Soft and collinear augmentation during training; input perturbations as in the CMS UParT tagger.
What approving this plan starts
Steps 0 and 1 only, for NNjjBj in the Z channel, as a new directory summer26/NNjjBj2/ copied from NNjjBj. Existing nets, ntuples and makefriends.cpp stay untouched. I smoke-test on CPU with a few thousand events, then show the condor_submit lines; jobs are submitted only when you say so.

Files
New summer26/NNjjBj2/: dataset.py (event frame, oriented constituents, log features, optional joint channels), model.py (context into phi, pooling, latent size), train.py (folds, seeds, out-of-fold scores), inference.py (weighted ROC, fixed-background working points, Pythia transfer), plus the usual train_gpu.*.
New summer26/compare_nets.py: paired bootstrap on out-of-fold scores of two nets.
Later steps only: makefriends.cpp, root_to_hdf5.py, run.summer26.sh.
Verification
With the new input options switched off and one fold, NNjjBj2 reproduces the NNjjBj Herwig AUC (0.870, Z) within the statistical error.
Input checks: rotating every phi by a random angle, reflecting eta, and swapping the tag-jet slots leave the network inputs unchanged to float precision.
Result table for old and new: AUC, eff_S at eff_B of 10% and 1%, on Herwig (out-of-fold) and on Pythia transfer, with paired-bootstrap errors, plus zeta.
Caveats
The expected sizes of gains in Steps 1 to 4 are my judgement from the cited results; none has been measured on your samples.
The argument that DeepSets suits colour flow (independent soft emissions make the log-likelihood ratio a sum over emissions) is my extension of the single-emission result in arXiv:2006.10480, which does not treat multiple emissions.
Side remark on studies/O_jet.md: the sum of z_i^2 r_i^2 is infrared safe but not collinear safe (a collinear splitting changes the sum of z^2), so it should not be called IRC safe.
