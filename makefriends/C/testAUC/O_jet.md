# O_jet — Learning a Classical Jet Observable from DeepSets

## Motivation

The VBF H→invisible analysis uses two jet taggers — DeepSets (`NNpull`) and a Transformer
(`NNtrans`) — that both reach AUC ≈ 0.813 despite very different architectures.  A feature
ablation study showed that removing the 6 jet-level scalars (which include the classical pull
vector components `jetPVM`, `jetSPVA`) costs only ~0.9% AUC, and removing the constituent
pull weight `jcsW` costs essentially nothing.  This raised the question: **what is the model
actually computing?**

To answer it, we trained a minimal `DeepSets (no scale)` variant using only 3 constituent
features per jet and then ran an interpretability pipeline to extract the learned observable.

---

## The Minimal Model

**Architecture:** standard DeepSets — `phi` (per-constituent MLP) → masked sum pool → `rho`
(jet-level MLP) → logit.

**Input features (3 per constituent, up to 80 per jet):**

| Symbol | Variable | Definition |
|--------|----------|------------|
| $\Delta\eta_i$ | `jcsDEta` | pseudorapidity of constituent $i$ relative to the jet axis |
| $\Delta\phi_i$ | `jcsDPhi` | azimuthal angle of constituent $i$ relative to the jet axis |
| $z_i$ | `jcsPt/jetPt` | transverse momentum fraction: $z_i = p_{T,i} / p_{T}^{\rm jet}$ |

where $r_i = \sqrt{\Delta\eta_i^2 + \Delta\phi_i^2}$ is the angular distance of constituent $i$
from the jet axis.

**Result:** AUC = 0.7810, compared to 0.8132 for the full-feature model.

---

## Interpretability Pipeline

The DeepSets forward pass factorises as:

$$\text{logit} = \rho\!\left(\sum_i \phi(x_i)\right)$$

The effective contribution of constituent $i$ to the logit is $w \cdot \phi(x_i)$, where

$$w = \mathbb{E}\!\left[\frac{\partial\,\rho(h)}{\partial h}\right]$$

is the mean Jacobian of $\rho$ with respect to its input $h$, averaged over the test set
(27 077 jets).  This projects the 64-dimensional phi output onto a single discriminating
direction, giving a **scalar per-constituent kernel**

$$f(\Delta\eta_i,\,\Delta\phi_i,\,z_i) \;\approx\; w \cdot \phi(x_i).$$

We then fit $f$ with an 18-term polynomial basis in $(\Delta\eta, \Delta\phi, z, r)$
using Ridge regression (α = 10⁻⁴) on 200 000 randomly sampled valid constituents from the
test set.

### Figure 1 — Effective kernel heatmaps and radial profiles

[phi_kernel_heatmap.pdf](phi_kernel_heatmap.pdf)

2×2 panel: $w \cdot \phi(\Delta\eta, \Delta\phi)$ at fixed $z \in \{0.05, 0.10, 0.20, 0.50\}$.
Second page: radial profiles $w \cdot \phi$ vs $r$ (mean ± std over angles) for each $z$.

### Figure 2 — Gradient attribution maps

[phi_gradient_maps.pdf](phi_gradient_maps.pdf)

For 5 000 test jets: $|\partial\,\text{logit}/\partial\Delta\eta| + |\partial\,\text{logit}/\partial\Delta\phi|$
as a hexbin in $(\Delta\eta, \Delta\phi)$ space (signal vs background), and
$|\partial\,\text{logit}/\partial z|$ vs $z$.

---

## Polynomial Fit Results

The fit quality is very high: **R² = 0.978** (train and validation), confirming that the
linearisation $w \cdot \phi(x)$ captures nearly all of the model's per-constituent logic.

### Top polynomial terms

| Rank | Term | Coefficient |
|------|------|-------------|
| 1 | $z^2 r^2$ | −133.15 |
| 2 | $z\,\Delta\eta$ | −70.48 |
| 3 | $z^2 r$ | +38.26 |
| 4 | $z\,r$ | −12.54 |
| 5 | $z\,r^2$ | −4.60 |
| 6 | $z\,\Delta\phi$ | −2.74 |
| 7 | $z$ | +1.63 |
| 8 | $z^2$ | +0.90 |
| 9 | $r$ | −0.35 |
| 10 | $r^2$ | +0.24 |

The two dominant terms together already define a physically transparent observable:

$$O_{\rm jet} = \sum_i \left(-133\, z_i^2\, r_i^2 \;-\; 70\, z_i\,\Delta\eta_i \;+\; \cdots \right)$$

- **$z_i^2 r_i^2$** — a pT²-weighted energy–energy correlator, sensitive to the radial
  spread of the hardest constituents.  Unlike the classical pull vector (which weights by
  $z_i r_i$), the quadratic pT weight suppresses soft constituents more aggressively.
- **$z_i\,\Delta\eta_i$** — a pT-weighted η asymmetry.  This term is odd in $\Delta\eta$
  and therefore sensitive to whether the jet is "pulled" toward positive or negative η,
  encoding colour-flow directionality along the beam axis.

All remaining significant terms also carry at least one power of $z$, confirming that the
model relies on **momentum-weighted spatial moments** rather than pure geometry.

---

## AUC Comparison

| Observable | AUC |
|------------|-----|
| DeepSets (no scale) — full network | **0.7810** |
| $O_{\rm jet}$ — polynomial (this work) | **0.7621** |
| $\|\vec{t}\,\|$ — pull vector magnitude (`jetPVM`) | 0.5714 |
| $\|\theta_s\|$ — signed pull angle (`jetSPVA`) | 0.5470 |

The polynomial observable recovers **97% of the gap** between random (0.5) and the full
network, while the classical pull-vector components are far weaker.

### Figure 3 — ROC comparison

[observable_roc.pdf](observable_roc.pdf)

---

## Proposed Classical Observable

The simplest analytically-defined observable suggested by the fit is:

$$\boxed{O_{\rm jet} = \sum_{i=1}^{N_c} z_i^2\, r_i^2}$$

where $z_i = p_{T,i}/p_{T}^{\rm jet}$ and $r_i = \sqrt{\Delta\eta_i^2 + \Delta\phi_i^2}$.

This is a **pT²-weighted jet width squared** — an IRC-safe, infrared and collinear safe
observable that can be computed from calorimeter or particle-flow constituents, measured
at detector level, and unfolded to particle level with standard methods (e.g. OmniFold or IBU).

Including the $z\,\Delta\eta$ term would require specifying an η sign convention relative to
the colour partner (as the original pull vector does), which is natural in VBF but makes the
observable less generic.

---

## Files

| File | Description |
|------|-------------|
| [interpret_utils.py](interpret_utils.py) | Shared helpers (data loading, w computation, poly basis) |
| [phi_visualize.py](phi_visualize.py) | Kernel heatmaps |
| [phi_gradient.py](phi_gradient.py) | Gradient attribution |
| [phi_pysr.py](phi_pysr.py) | Polynomial fit + PySR attempt |
| [observable_roc.py](observable_roc.py) | ROC comparison |
| [poly_model.npz](poly_model.npz) | Saved Ridge coefficients |
| [pysr_summary.txt](pysr_summary.txt) | Full coefficient table and fit statistics |
| [phi_kernel_heatmap.pdf](phi_kernel_heatmap.pdf) | Figure 1 |
| [phi_gradient_maps.pdf](phi_gradient_maps.pdf) | Figure 2 |
| [observable_roc.pdf](observable_roc.pdf) | Figure 3 |
