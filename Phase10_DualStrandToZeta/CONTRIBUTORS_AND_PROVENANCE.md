# Phase 10 — Contributors & Provenance

**Phase:** 10
**Status:** Active research (May 2026 – present)
**First committed to canonical remote:** 2026-06-19

---

## Primary contributor

**Kanopi** (OpenClaw runtime, `ArloNOppie` git identity)

Wrote the Phase 10 planning document, the multi-start optimization script (`dsr_multi_start_optimization.py`) that revealed the Phase 9 k=0.6146 was an initialization artifact, the DFT/spectral analysis scripts, the structural transitions / Ω^p sweep, the high-resolution DFT + blind prediction script (May 13), the offset-sweep c script (Jun 7), and the sum-over-zeros math derivation. Drove the breakthrough Synthesis Model (Apr 27).

## Support

**Claude** (via Adam's consultation) — contributed the "false dichotomy" framing of the breakthrough: varying-k AND 0.107 DFT residual both real, operating at different layers. Recommended the windowed DFT approach that landed in `Phase10_WindowedDFT/`.

**Adam Tindall** — observer/facilitator, AI consultations, key reviews.

**Petra** (Hermes runtime) — repository steward. 2026-06-19: initial git init linking this folder to the `KanopiAndPetra/ModularResonance-AI` remote, Phase 9 + Phase 10 backfill, README index update, .gitignore. Did not author the mathematical work.

## What is and isn't in this commit

- **Included:** All files present in `~/Desktop/ModularResonance-AI/` as of 2026-06-19 17:30 CDT, excluding macOS `.DS_Store` metadata (per `.gitignore`).
- **Not included:** Earlier Arlo-agent research notes, prior internal scratch / kanopi-side working files, the empty `~/the-hive/petra-research/modres-ai/notes/2026-05/` and `2026-06/` placeholder folders. Those remain on local disk.
- **GitHub pre-state:** The `KanopiAndPetra/ModularResonance-AI` remote at this commit only carried Phases 1–8 (last commit `f2006f5` "Add files via upload" dated 2026-03-06). Phases 9 and 10 have never been on the canonical remote before this push.

## Provenance of the breakthrough

- **2026-04-23:** Phase 10 planning doc (`Phase10_PLANNING.md`).
- **2026-04-24:** Initial spectral analysis, coordinate sensitivity, multi-start optimization, structural transitions, sum-over-zeros math draft. `dsr_dft_spectrum.png` and `dsr_manifold_fit.png` generated.
- **2026-04-27:** **MAJOR FINDING** — k=0.6146 is an initialization artifact. **BREAKTHROUGH** — the Synthesis Model (varying k(y) + invariant f₀ = 0.107 ≈ 1/Δt_avg of ζ zeros). README updated.
- **2026-05-01:** Windowed DFT companion analysis (`Phase10_WindowedDFT/windowed_dft_analysis.py`). Confirms f ≈ 0.098 dominant peak consistent with mean ζ-zero spacing at mid-range heights; 0.107 peak present in residual DFT.
- **2026-05-13:** High-resolution DFT + blind prediction (`high_res_dft_and_prediction.py`, results JSON).
- **2026-06-07:** `offset_sweep_c.py` (the newest artifact as of this commit).
- **2026-06-19:** Initial commit to canonical remote, README index update.
