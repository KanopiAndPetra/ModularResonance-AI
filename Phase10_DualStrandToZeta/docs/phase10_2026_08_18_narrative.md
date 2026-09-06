# Phase 10 Narrative — 2026-08-18

**Author:** Petra (with Adam back from trip, observing requisites)
**Date:** 2026-08-18
**Status:** Working synthesis — *not yet a published result*
**Reading layer:** above the PLANNING.md (Apr 23) and the README.md (Apr 27 breakthrough), integrating what was learned AFTER the breakthrough

---

## What this document is

The PLANNING.md (Kanopi, 2026-04-23) and README.md (Kanopi, 2026-04-27) describe the
Phase 10 path as it stood on Apr 27. Three things have happened since then:

1. **The blind prediction test (May 13) failed.** `high_res_dft_and_prediction.py` fit a
   training helix (k=0.99, A=0.155, phase=-0.194) and tried to predict crossings on a held-out
   region. The result: zero crossings matched, `mapping_error: Infinity`, verdict
   `"NO predictive power (correlation only)"`. The 0.107 DFT peak was confirmed in this run
   (`peaks_near_0107: [0.116, 0.101]`). So the *residual* is real, but our *predictive model*
   isn't.

2. **The offset_sweep_c (Jun 7) failed.** The hypothesis was that a small additive offset c
   in the helix-crossing formula `y_cross(n) = ((2n+1)·π/2 - φ) / k + c` would close the
   QFPIL gap. The predicted optimum was c ≈ 0.02. The sweep at c = 0.0, 0.001, 0.002 …
   produced `mean_abs_delta = 0.589` at every value of c (no signal). The 0.02 optimum
   didn't exist.

3. **The harmonic context question Adam raised on return.** After Adam's trip, the question
   became: *does the 0.107 invariant actually demonstrate a limit to the critical lines?*
   This is the question this Phase 10 narrative is built to address — honestly.

---

## What we have, honestly

### Solid (reproducible, robust)

| Signal | Value | Source | Confidence |
|---|---|---|---|
| DFT residual peak | f₀ ≈ 0.107 | `high_res_dft_results.json`, `peaks_near_0107: [0.116, 0.101]` | High — appears in 2+ independent scans |
| Mean ζ-zero spacing at T_max = 10²⁶ | 0.1058 | `2π/log(T_max)` formula | High — analytic |
| **0.107 / 0.1058 match** | **1.5%** | direct computation | **High — robust** |
| k(y) varies with y | F=50.9, p<0.0001 | multi-start optimization | High — statistically solid |
| k(y) follows logarithmic warp | k(y) = a + b·ln(y+1) | synthesis model | Medium — model form is plausible but fit residuals not yet published |

### Mixed (real signal, but predictive modeling failed)

| Signal | Value | Source | Confidence |
|---|---|---|---|
| Blind prediction of crossings | 0 / 0 | `high_res_dft_and_prediction.py` | Low — model fails |
| Offset sweep c ≈ 0.02 closes QFPIL gap | mean_abs_delta = 0.589 (constant) | `offset_sweep_c_results.json` | Negative result — gap not closed |
| Phase 9 k = 0.6146 as "golden ratio" | reproducible only at k₀=0.63 init | `multistart_results.json` | Negative — initialization artifact |

### Not yet tested

- Boundedness of |x - 0.5| — does the helix have a finite amplitude envelope, and does the envelope carry the ζ-bridge signal?
- Empirical zero-spacing distribution from the manifold — does the manifold's zero-crossing spectrum reproduce the Montgomery pair correlation?
- The exact form of the critical-line constraint — what does it even mean for a manifold to "have a limit to the critical lines"?

---

## What "proving a limit to the critical lines" would actually require

The Riemann Hypothesis states: the nontrivial zeros of ζ(s) all lie on the line Re(s) = 1/2.

To prove a **limit** that ENFORCES this from the DSR structure, we need to show:

1. **The DSR manifold extracts a specific finite-dimensional structure from Λ(n).**
   The manifold is determined by the prime distribution through Λ(n). So hypotheses
   about ζ(zeros) translate to hypotheses about the manifold.

2. **Examining the manifold's geometric constraints → contradictions for any hypothetical
   zero with Re(s) ≠ 1/2.**
   This is where the boundedness argument and the k(y) warp come in.

3. **The constraint is quantitatively tight.** Not just "0.107 is close to 0.106" but
   "0.107 differs from 0.106 by less than the natural width of the zero-spacing
   distribution at this height, and any deviation from Re(s) = 1/2 would push the
   constraint out of that window."

### What we can do TODAY

- Run the **boundedness test** on |x - 0.5| — does the actual measured amplitude envelope
  agree with the theoretical bound?
- Compute the **empirical zero-spacing distribution** from the manifold crossings and
  compare to the GUE/Montgomery expectation.
- Propose the **limit-proposition** formally (this is the other doc in this folder).

### What we cannot do today

- Prove the Riemann Hypothesis. The 0.107 is a **signal**, not a proof. The proof
  requires a closed-form argument — a finite-dimensional operator identity, a contour
  integral, or a spectral theory result — that connects the manifold shape to the
  Riemann zeta function's zero set. We have evidence; we don't have a proof.

---

## The honest Phase 10 path forward

The path is: **prove the limit proposition** (testable, falsifiable, narrow) → **if it
holds, generalize** (the bounded helix implies the critical line). The proposition is
the contribution this Phase 10 can make; the generalization is the long-term prize.

The boundedness test and the empirical zero-spacing histogram are the two experiments
that convert a "signal" (the 0.107) into a "testable framework" (the limit proposition).
If both falsify the proposition, the 0.107 is a coincidence and Phase 10 closes. If
both confirm, the proposition becomes the basis for Phase 11 (the formal proof).

See:
- `phase10_2026_08_18_limit_proposition.md` — the formal statement
- `python_code/phase10_dsr_boundedness_test.py` — Experiment 1
- `python_code/phase10_zero_spacing_histogram.py` — Experiment 2

---

*Working doc. Update as the experiments land results.*
