# Roadmap: fractional Sobolev spaces on Lipschitz boundaries and boundary integral operators

Mathlib has the Bessel-potential scale on `ℝⁿ` (`memSobolev`, `besselPotential` in
`Mathlib/Analysis/Distribution/Sobolev.lean`). Tau Ceti has weak-derivative Sobolev spaces on
domains (`TauCeti.W1p`), the Lipschitz-domain predicate `Set.IsLipschitzDomainWith` with chart
traces (`TauCeti.W1p.lipschitzChartTrace`) and their comparison with Hausdorff surface measure,
and the Newtonian fundamental solutions `TauCeti.newtonianKernel` (`n ≥ 3`) and
`TauCeti.planarNewtonianKernel` (on `ℂ`). Geometric measure theory G6 supplies Gauss–Green on
Lipschitz domains. What is missing is the **energy-space theory of boundary value problems on
Lipschitz domains**: fractional Sobolev spaces `H^s` on the boundary for `|s| ≤ 1`, the **sharp
trace theorem** `H¹(Ω) → H^{1/2}(∂Ω)` (onto, with a bounded right inverse), the conormal derivative
in `H^{−1/2}(∂Ω)`, the **layer potentials** of the Laplacian with their jump relations, the four
**boundary integral operators** with their mapping properties, their **coercivity**, the
**Calderón projector**, and the well-posedness of the boundary integral equations of the
Dirichlet and Neumann problems.

The bar for "done": someone analysing boundary integral equations for the Laplacian on a Lipschitz
domain, in the setting of McLean's book (Chapters 3, 4, 6, 7, 8), finds the spaces, traces,
potentials and operators with the mapping and coercivity theorems they cite.

## Standing setting

- `n ≥ 2`, `E := EuclideanSpace ℝ (Fin n)`, Lebesgue measure. `Ω` is a bounded Lipschitz domain
  (`Set.IsLipschitzDomainWith`, PDE Lane A.6), `Γ := frontier Ω` with the surface measure
  `σ := μHE[n−1]` restricted to `Γ`, and `Ωᶜ := (closure Ω)ᶜ` the exterior domain.
- All function spaces are real-valued. Duality pairings on `Γ` extend the `L²(Γ, σ)` inner product.
- The theory is for `|s| ≤ 1` on `Γ` and `0 ≤ s ≤ 1` (with `s = 1` from Tau Ceti's `W1p`) on domains:
  this is the range that is intrinsic for Lipschitz boundaries.

## Milestone 1: fractional Sobolev spaces on `ℝⁿ` and on domains

1. **Sobolev–Slobodeckij spaces.** For an open `U ⊆ E`, `0 < s < 1`, `1 ≤ p < ∞`: the Gagliardo
   seminorm `[u]_{s,p,U} := (∫_U ∫_U |u x − u y|^p / ‖x − y‖^{n + s p})^{1/p}` and
   `W^{s,p}(U) := {u ∈ Lᵖ(U) | [u]_{s,p,U} < ∞}` as a Banach space (Hilbert for `p = 2`, written
   `H^s(U)`); the basic API (restriction, translation and dilation laws, composition with
   bi-Lipschitz maps with explicit constants, density of `C_c^∞(E)` in `W^{s,p}(E)`).
2. **Comparison with Mathlib's scale.** `H^s(E)` equals Mathlib's Bessel-potential `H^{s,2}` (via
   `memSobolev_two_iff_fourier`) with equivalent norms, for `0 < s < 1`.
3. **Embeddings between orders.** `W^{1,p}(U) ⊆ W^{s,p}(U)` for bounded Lipschitz `U` and `s < 1`,
   with the norm bound; `W^{s,p}(U) ⊆ W^{t,p}(U)` for `t < s`.
4. **Extension.** For bounded Lipschitz `U`, a bounded extension operator `H^s(U) → H^s(E)`
   for `0 < s ≤ 1` (by Lipschitz charts and reflection, consuming Tau Ceti's first-order
   `W1p.exists_extensionL_of_isLipschitzDomainWith` as the model).
5. **The spaces `H̃^s(U)`** (closure of `C_c^∞(U)` in `H^s(E)`) and the duality
   `H^{−s}(U) := (H̃^s(U))'`, for `0 ≤ s ≤ 1`.

## Milestone 2: Sobolev spaces on the boundary

6. **`H^s(Γ)` for `0 < s < 1`**, defined intrinsically: the Gagliardo seminorm on `Γ` with respect
   to `σ`, `[g]_{s,Γ}² := ∫_Γ ∫_Γ |g x − g y|² / ‖x − y‖^{n − 1 + 2s} dσ dσ`, and
   `H^s(Γ) := {g ∈ L²(Γ, σ) | [g]_{s,Γ} < ∞}`, a Hilbert space. Prove the equivalence with the
   chart definition (pulling back along each Lipschitz chart to `ℝⁿ⁻¹` and using item 1), with
   constants depending on the chart data.
7. **`H¹(Γ)`** through the tangential gradient (defined `σ`-a.e. from the Lipschitz charts and the
   approximate tangent spaces of GMT), with `H¹(Γ) ⊆ H^s(Γ) ⊆ L²(Γ)` continuously.
8. **Negative orders.** `H^{−s}(Γ) := (H^s(Γ))'` for `0 < s ≤ 1`, with `L²(Γ)` embedded by the inner
   product, so that `H^s(Γ) ⊆ L²(Γ) ⊆ H^{−s}(Γ)` is a Gelfand triple.
9. **Open pieces of the boundary.** For a relatively open `Γ₀ ⊆ Γ` with Lipschitz relative
   boundary: `H^s(Γ₀)`, `H̃^s(Γ₀)` (extension by zero to `Γ` stays in `H^s(Γ)`) and their duals, for
   `|s| ≤ 1`.
10. **Polyhedral boundaries.** When `Ω` is a bounded open polytope with Lipschitz boundary (`Γ` a
    finite union of `(n−1)`-dimensional convex polytopes), `σ` restricted to each face is the
    face's `(n−1)`-dimensional Lebesgue measure, and `g ∈ H^s(Γ)` for `0 < s < 1/2` if and only if
    each face restriction is in `H^s` of the face. Include this characterization and its failure
    for `s ≥ 1/2` (a function that is constant on each face but not globally is not in
    `H^{1/2}(Γ)`).

## Milestone 3: the sharp trace theorem and the conormal derivative

11. **Trace onto `H^{1/2}(Γ)`.** The trace operator `γ : H¹(Ω) → H^{1/2}(Γ)` (agreeing `σ`-a.e. with
    Tau Ceti's chart traces and with boundary values of functions continuous on `closure Ω`) is
    bounded and **surjective**, with a bounded linear right inverse `H^{1/2}(Γ) → H¹(Ω)`; its kernel
    is `W1p0 volume Ω 2`. The same for the exterior domain with `H¹_loc(Ωᶜ)`.
12. **Conormal derivative.** For `u ∈ H¹(Ω)` with `Δu ∈ L²(Ω)` (weakly), define
    `∂ₙu ∈ H^{−1/2}(Γ)` by Green's first identity
    `⟨∂ₙu, γv⟩ := ∫_Ω ∇u · ∇v + ∫_Ω Δu · v` for all `v ∈ H¹(Ω)`, prove it is well defined and
    bounded, and that it agrees with `∇u · ν` (outward normal of G6) when `u` is `C¹` on
    `closure Ω`. State and prove Green's second identity in this generality.
13. **Gluing across `Γ`.** `u ∈ H¹_loc(E)` if and only if its restrictions to `Ω` and `Ωᶜ` are in
    `H¹_loc` with equal traces; and for such `u` with `Δu ∈ L²_loc` on both sides, the
    distributional Laplacian of `u` on `E` is the sum of the two sides minus the jump of the conormal
    derivatives as a distribution on `Γ`.

## Milestone 4: layer potentials and boundary integral operators

Let `G` be the fundamental solution of `−Δ` on `E`: `TauCeti.newtonianKernel n` for `n ≥ 3`, and
for `n = 2` the logarithmic kernel transported from `TauCeti.planarNewtonianKernel` along the
isometry `EuclideanSpace ℝ (Fin 2) ≃ ℂ`; make that transport a target.

14. **Potentials.** The single layer potential `SL ψ (x) := ∫_Γ G(x − y) ψ y dσ y` and the double
    layer potential `DL φ (x) := ∫_Γ ∂_{ν(y)} G(x − y) φ y dσ y` for `x ∉ Γ`, first for `ψ, φ ∈
    L²(Γ)`, then extended by continuity: `SL : H^{−1/2}(Γ) → H¹_loc(E)` and
    `DL : H^{1/2}(Γ) → H¹_loc(Ω) × H¹_loc(Ωᶜ)` bounded, harmonic off `Γ`, with the behaviour at
    infinity (`O(‖x‖^{2−n})`, and the logarithmic case in `n = 2`).
15. **Jump relations.** `[γ SL ψ] = 0`, `[∂ₙ SL ψ] = −ψ`, `[γ DL φ] = φ`, `[∂ₙ DL φ] = 0`, where
    `[·]` is the jump across `Γ` (exterior minus interior trace).
16. **The four boundary integral operators** `V := γ SL`, `K := ½(γ⁺ + γ⁻) DL`,
    `K' := ½(∂ₙ⁺ + ∂ₙ⁻) SL`, `W := −∂ₙ DL`, bounded
    `V : H^{−1/2}(Γ) → H^{1/2}(Γ)`, `K : H^{1/2}(Γ) → H^{1/2}(Γ)`, `K' : H^{−1/2}(Γ) → H^{−1/2}(Γ)`,
    `W : H^{1/2}(Γ) → H^{−1/2}(Γ)`, with `K'` the adjoint of `K` and `V`, `W` self-adjoint
    (Costabel 1988: the proof uses the `H¹` bounds of the potentials, the trace theorem and
    Green's identities, not harmonic analysis on `Γ`). Identify `V` and `K` with their integral
    formulas on `L²(Γ)` functions where the integrals converge absolutely.
17. **Representation formula.** For `u ∈ H¹(Ω)` harmonic: `u = SL (∂ₙu) − DL (γu)` in `Ω`, and the
    exterior analogue with the radiation condition at infinity.

## Milestone 5: coercivity, the Calderón projector, and well-posedness

18. **Coercivity.** For `n ≥ 3`: `⟨V ψ, ψ⟩ ≥ c ‖ψ‖²_{H^{−1/2}(Γ)}`. For `n = 2`: the same on the
    subspace `⟨ψ, 1⟩ = 0`, and on all of `H^{−1/2}(Γ)` when `diam Ω < 1` (the scaling condition).
    `⟨W φ, φ⟩ ≥ c ‖φ‖²_{H^{1/2}(Γ)}` on `{φ | ⟨φ, 1⟩ = 0}`, with `ker W` the constants for connected
    `Γ`.
19. **The Calderón projector** `C := [[½ − K, V], [W, ½ + K']]` on `H^{1/2}(Γ) × H^{−1/2}(Γ)`:
    `C² = C`, its range is the set of Cauchy data `(γu, ∂ₙu)` of harmonic `u ∈ H¹(Ω)`, and the
    resulting identities `V W = ¼ − K²`, `W V = ¼ − K'²`, `K V = V K'`, `W K = K' W`.
20. **Boundary integral equations.** Unique solvability, with stability bounds, of the first-kind
    equations `V ψ = g` (Dirichlet problem; `n ≥ 3`, or `n = 2` with the scaling condition) and
    `W φ = h` on the mean-zero subspaces (Neumann problem, `⟨h, 1⟩ = 0`), by Lax–Milgram from item
    18; and the second-kind equations `(½ + K) φ = g` and `(½ − K') ψ = h` through the Calderón
    identities, on the appropriate subspaces. Show that the potentials built from these solutions
    solve the interior and exterior Dirichlet and Neumann problems for the Laplacian in `H¹`.

## Boundaries

- **In scope:** everything above, for the Laplacian, real-valued functions, bounded Lipschitz
  domains in `EuclideanSpace ℝ (Fin n)`, `n ≥ 2`.
- **Consumed:** PDE Lane A (`W1p`, the Lipschitz-domain predicate, chart traces, extension), PDE
  Lane D (Lax–Milgram for the domain problems), Geometric measure theory G6 (Gauss–Green, outward
  normal, surface measure), Mathlib's Bessel-potential scale and Lax–Milgram, and Tau Ceti's
  Newtonian kernels.
- **Not here:** `L^p(Γ)` boundedness of the double layer operator on Lipschitz boundaries and the
  `L²` theory of the Dirichlet and Neumann problems (Coifman–McIntosh–Meyer, Verchota, Dahlberg–Kenig);
  `H^s(Γ)` for `|s| > 1` (it needs smoother boundaries); operators other than the Laplacian
  (Helmholtz, elasticity, Stokes); and discretization of boundary integral equations on partitions
  of `Γ`.

## Acceptance criteria

- On the boundary of the unit cube in `ℝ³` and of a triangle in `ℝ²`: the face characterization of
  item 10, the trace theorem with its right inverse, and the coercivity of `V` (with the scaling
  condition for the triangle) are instantiated in a test file.
- On the unit ball in `ℝ³`: `V` and `W` diagonalize on spherical harmonics with the classical
  eigenvalues `1/(2ℓ + 1)` and `ℓ(ℓ + 1)/(2ℓ + 1)`, checked as theorems for `ℓ = 0, 1`.
- Every space and operator has its basic API (linearity, norms, duality, restriction to subsets
  of `Γ`), not only the headline theorems.

## References

- W. McLean, *Strongly Elliptic Systems and Boundary Integral Equations*, Cambridge, 2000,
  Chapters 3, 4, 6, 7, 8.
- M. Costabel, *Boundary integral operators on Lipschitz domains: elementary results*, SIAM J. Math.
  Anal. 19 (1988).
- S. A. Sauter, C. Schwab, *Boundary Element Methods*, Springer, 2011, Chapters 2–3.
- O. Steinbach, *Numerical Approximation Methods for Elliptic Boundary Value Problems*, Springer,
  2008, Chapters 2, 6.
- E. Di Nezza, G. Palatucci, E. Valdinoci, *Hitchhiker's guide to the fractional Sobolev spaces*,
  Bull. Sci. Math. 136 (2012).
- J. Nečas, *Direct Methods in the Theory of Elliptic Equations*, Springer, 2012 (translation of
  the 1967 edition), Chapter 2.
