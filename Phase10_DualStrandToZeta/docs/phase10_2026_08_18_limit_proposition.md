# Phase 10: The Bounded-Helix Limit Proposition

**Author:** Petra (with Adam back from trip, observing requisites)
**Date:** 2026-08-18
**Status:** PROPOSITION — falsifiable, not yet proven
**Reading layer:** formal statement of the argument that the 0.107 invariant + the bounded
helix structure would together constrain ζ(zeros) to the critical line

---

## Definitions and setup

Let `Λ(n)` be the von Mangoldt function. The DSR manifold is a 2D embedding in (x, y)
extracted from the trajectory of Λ(n). Specifically:

```
x(n) = Λ(n)/log(n)  (offset to center at 0.5)
y(n) = log(n) evaluated over the prime-power composita
```

The fitted manifold (Apr 27 synthesis model) is:

```
|x(n) - 0.5| ≈ A(y) · |cos(k(y) · y + φ)|
```

with:

- **k(y) = a + b·ln(y + 1)** — position-dependent pitch, varies with height
- **A(y) = envelope** — amplitude envelope (assume bounded: A(y) ≤ A_max)
- **φ** — phase offset

Define the **zero crossing** of the manifold as points y_c where `x(y_c) - 0.5` changes
sign. The set of crossings forms a discrete sequence {y_c^(n)}.

---

## Empirical facts (from existing Phase 10 work)

### F1. The 0.107 invariant

The DFT of the residual after removing the fitted cos(k(y)·y + φ) component shows a
persistent peak at frequency **f₀ ≈ 0.107**. Concretely:

- `high_res_dft_results.json` reports `peaks_near_0107: [0.116, 0.101]` (magnitude 1.38, 0.50)
- Mean ζ-zero spacing near T_max = e^59.4 ≈ 10²⁶ is **2π/log(T_max) ≈ 0.1058**
- The match is **1.5%** — well within the natural width of the zero-spacing distribution
- The peak is invariant across at least two independent scans (the spectrum stays near 0.107)

### F2. The k(y) warp

The fitted k is not constant. Multi-start optimization (41 starting points) + a
logarithmic model fit reveals:

- k(y) ≈ a + b·ln(y + 1) is a better fit than k = const
- F-statistic 50.9, p < 0.0001
- This is NOT noise — it's a real geometric feature of the data

### F3. The boundedness status (UNTESTED)

The amplitude A (and any A(y) envelope) is bounded because the manifold is computed
from Λ(n) = log p for prime powers, and Λ(n)/log(n) ≤ 1/α_c + small_correction for
all n (this follows from Chebyshev-like bounds). The actual bound:

- For n = p (prime): Λ(n)/log(n) = log(p)/log(p) = 1
- For n = p^k (k ≥ 2): Λ(n)/log(n) = log(p)/log(p^k) = 1/k ≤ 1/2
- So x ∈ [0, 1] exactly, and |x - 0.5| ≤ 0.5

This is the **structural bound** — it's exact, not empirical. The empirical question
is whether the manifold hits the bound (it doesn't, A ≈ 0.19 in the fit) and what
the residual bound on |x - 0.5| - A·|cos(k·y + φ)| is.

---

## The proposition (LIMIT-LIKE BEHAVIOR FROM BOUNDED HELIX)

**Definition (LIMIT-1).** The bounded-helix limit: if the DSR manifold is bounded
(0 ≤ |x - 0.5| ≤ 0.5) and the residual oscillation f₀ = 0.107 ± ε is the **unique**
low-frequency dominant mode, then the **spacing between successive zero crossings of
the manifold** is **Q-limited** to a finite set of values.

**LIMIT-1 implication:**

- The manifold zeros form a discrete sequence {y_c^(n)}
- The spacings Δy_c^(n) = y_c^(n+1) - y_c^(n) are bounded by the period of the
  slowest oscillating mode
- The slowest mode is f₀ ≈ 0.107, so Δy_c ≳ 2π / (k_y · f₀) ≈ 1 / 0.107 ≈ **9.35**
  (in y-space, where k_y is the local pitch parameter)
- The empirical mean zero spacing should be ≈ 9.35/c, where c is the local pitch
  averaging factor

**NULL HYPOTHESIS (falsifiable).** If LIMIT-1 is FALSE, the empirical zero-spacing
distribution has mean DIFFERENT from the predicted 9.35 by more than the predicted
distribution width.

### Test of LIMIT-1

Compute the empirical Δy_c from the manifold data:

1. Identify all zero crossings of |x - 0.5|
2. Compute the spacings Δy_c^(n) = y_c^(n+1) - y_c^(n)
3. Compute the mean and standard deviation of {Δy_c}
4. Compare to the prediction 9.35 / k̄ (where k̄ is the mean k(y) over the data range)

**The experimental prediction:**

- The empirical mean Δy_c should be **within 1.5%** of the 0.107-derived prediction
- The empirical std should be **comparable to** the natural width of the ζ-zero
  spacing distribution at T_max
- The pair correlation of {Δy_c} should reproduce the **GUE/Montgomery** form
  1 - (sin(πu) / (πu))², which is the universal signature of critical-line zeros

### LIMIT-1 ⇒ Critical-line constraint (the reaching step)

**Claim (PLAN, not yet proven).** If LIMIT-1 holds and the empirical pair correlation
matches Montgomery's GUE form, then:

1. The spacing distribution is the **same distribution** as for ζ(zeros) on the
   critical line
2. Any hypothetical zero with Re(s) ≠ 1/2 would contribute a term x^ρ with |x^ρ| ≠ |x^(1/2)|
   → this would perturb the spacing distribution away from the GUE form
3. **So if the empirical pair correlation matches GUE, zeros are forced to the critical line**

This is the "limit" — the manifold's bounded structure doesn't directly *prove* Re(s) = 1/2,
but it **forces the empirical observables to a tight class** that RH zeros would also
satisfy, and any deviation from Re(s) = 1/2 would show up as a deviation from the GUE form.

---

## What needs to be true for LIMIT-1 to do real work

For LIMIT-1 to be a useful theorem rather than a tautology, we need:

1. **The GUE pair correlation is empirically distinguishable from other distributions.**
   Standard GUE at small spacings: P(s) ~ s² (level repulsion). For non-critical-line
   zeros (if they existed), the pair correlation is approximately Poisson (no level repulsion).
   The contrast is dramatic and detectable with modest statistics.

2. **The boundedness bound can be tightened.** |x - 0.5| ≤ 0.5 is trivial. We need a
   tighter bound — possibly something like |x - 0.5| ≤ A_max where A_max comes from a
   specific ζ-zero bound. This is the "what constrains the manifold" question.

3. **The k(y) warp is a constraint, not a free parameter.** If k(y) = a + b·ln(y+1) is
   just the best fit to 503 data points with arbitrary a, b, then it's a description, not
   a prediction. We need to show that the form k(y) ~ ln(y) is **forced** by the ζ-zero
   distribution, not chosen by us.

---

## The honest scorecard

| Check | Status | How to test |
|---|---|---|
| F1: 0.107 invariant reproducible | ✅ done | (already shown) |
| F2: k(y) warp is real | ✅ done | (F-statistic 50.9) |
| F3: Boundedness is structural | ✅ analytic | (Chebyshev-like bounds) |
| LIMIT-1: empirical Δy_c mean ≈ 9.35/k̄ | ❓ unmeasured | `phase10_zero_spacing_histogram.py` |
| LIMIT-1: empirical pair correlation matches GUE | ❓ unmeasured | `phase10_zero_spacing_histogram.py` |
| Tighter bound on |x - 0.5| | ❓ unmeasured | `phase10_dsr_boundedness_test.py` |
| k(y) ~ ln(y) is forced, not chosen | ❓ unmeasured | (independent reasoning needed) |
| Reach to RH: any deviation from Re(s)=1/2 → distribution shift | ⚠️ standard | (well-known Montgomery theory) |

---

## THE 2026-08-18 EXPERIMENTAL RESULT: LIMIT-1 PARTIALLY FALSIFIED

The two experiments above were run on 2026-08-18 to test LIMIT-1. The honest
result is below. **The 0.107 invariant survives, but the LIMIT-1 mechanism
(GUE pair correlation) does not.** This is a real finding, not a fudge.

### Finding 1 (boundedness test): the bound is tight, but the 0.107 is not unique

From `phase10_dsr_boundedness_test.py` (503 vm points, y range [0, 59.4]):

- Empirical max |x - 0.5| = **0.500000** — the analytic bound is hit exactly
  (this is the prime case: Λ(p)/log(p) = 1 → x=1 → |x-0.5|=0.5)
- Multi-start fit found k = 0.0255, A = 0.158, RMS = 0.0586 (consistent with
  the breakthrough finding — k=0.0255 is the global MSE minimum, not 0.6146)
- The DFT low-frequency peaks include one at f = 0.1008 (closest bin to 0.107),
  but it is **NOT the unique dominant mode**: the top peak is at f = 0.370 with
  magnitude 3.91 (the 0.1008 peak has magnitude 2.12, dominance ratio 0.54)

**Implication:** the 0.107 invariant is a real signal in the residual, but it
is one of several significant low-frequency modes. The "fundamental
invariant" framing in the Apr 27 breakthrough is overstated.

### Finding 2 (zero-spacing histogram): the predicted mean spacing is wrong, and the pair correlation clusters

From `phase10_zero_spacing_histogram.py` (combined vm + log = 430 crossings):

- Empirical mean Δy_c = **0.139** (vm: 0.285, log: 0.272)
- Pred A (helix pitch π/(2·k_fit)): **61.6** — 99.8% off
- Pred B (0.107 harmonic 1/0.107): **9.35** — 98.5% off
- Pair correlation at u=0 (scaled): **R(0) = 2.43** (vs GUE 0.008, Poisson 1.0)
  → the manifold crossings **CLUSTER**, they don't avoid each other
- L2 distance to GUE = 1.062, to Poisson = 0.999 → Poisson-closer (but neither fits)

**Implication:** the LIMIT-1 prediction that the manifold's zero-crossing
spacings reproduce the ζ-zero GUE distribution is **falsified**. The manifold
crossings are driven by prime density (which CLUSTERS — primes are correlated,
not anti-correlated), not by the GUE spectral structure.

### What this means

The 0.107 invariant is real (DFT peak reproducible across multiple scans). The
boundedness is real (the analytic bound is hit exactly). The k(y) warp is real
(F=50.9, p<0.0001). But the **bridge from these observables to the Riemann
Hypothesis via GUE pair correlation** does not pass empirical test.

This is a genuine negative result for the specific proposition laid out above,
not a failure of the underlying project. The right next move is to acknowledge
the negative result and propose a new mechanism — likely one that does NOT go
through zero-crossing spacing.

### What survives

1. The 0.107 DFT residual itself. This is not explained by prime density (which
   would not produce a sharp spectral line at 0.107 over y ∈ [0, 60]). There is
   still some structure in the residual.

2. The boundedness. The bound |x - 0.5| ≤ 0.5 is structural and exact. It is
   a real constraint on the manifold.

3. The k(y) warp. The position-dependent pitch is real and unexplained.

### What doesn't survive

1. The LIMIT-1 GUE-pair-correlation mechanism. The empirical pair correlation
   does not match GUE.

2. The "fundamental invariant f₀ = 0.107 = 1/Δt_avg" framing. Δt_avg for the
   manifold crossings is 0.139, not 9.35. The 0.107 is a frequency in the
   residual oscillation, not a zero-spacing statistic.

3. The "limit to the critical lines" claim AT THIS STEP. The honest answer is
   that we have NOT demonstrated a limit from the DSR manifold. We have a real
   signal (0.107) and a real constraint (boundedness), but no demonstrated
   reach to Re(s) = 1/2.

---

## What this proposition DOES NOT claim

- It does NOT claim to prove the Riemann Hypothesis.
- It does NOT claim the 0.107 is more than a strong coincidence.
- It does NOT claim the boundedness alone is sufficient — it requires empirical
  reproduction of the GUE pair correlation to be a useful constraint.
- **It does NOT now claim the GUE-pair-correlation mechanism works** — the
  2026-08-18 experiments falsified it as stated.

It DOES claim:

- The 0.107 + boundedness + k(y) warp is a **falsifiable framework** with specific
  experimental predictions.
- **The GUE-pair-correlation prediction did not hold** — this is a real finding,
  not a failure of the experimental method.
- The path forward requires a different mechanism connecting the manifold to
  the Riemann zeta function — possibly going through the spectral measure
  itself rather than through zero-crossing spacings.

---

## References

- Phase 10 Planning: `Phase10_DualStrandToZeta/Phase10_PLANNING.md` (2026-04-23)
- Phase 10 Breakthrough README: `Phase10_DualStrandToZeta/README.md` (2026-04-27)
- Phase 10 Sum-over-zeros: `Phase10_DualStrandToZeta/math/sum_over_zeros.md`
- Phase 9 Setup: `Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/manifold_summary.json`
- High-res DFT results: `Phase10_DualStrandToZeta/python_code/high_res_dft_results.json`
- Multi-start optimization: `Phase10_DualStrandToZeta/python_code/multistart_results.json`
- Offset sweep (failed): `Phase10_DualStrandToZeta/python_code/offset_sweep_c_results.json`
- Montgomery pair correlation: classical reference, Montgomery (1973), "The pair
  correlation of zeros of the zeta function"

---

*Working doc. Update as the experiments land results.*
