#!/usr/bin/env python3
"""
Phase 10 Re-run v3 — Additive Offset c Sweep (direct residual-crossing match)
=============================================================================

Falsification test of the calibration-bias hypothesis from the 04:50
data-provenance archaeology (2026-06-06), RE-FRAMED for the actual data structure.

Key insight (from v2 archaeology): the trajectory data has:
- 228/248 (92%) "zero residuals" — pairs where midpoint = 0.5 exactly
- 20/248 (8%) "non-zero residuals" — pairs where midpoint ≠ 0.5
  (these cluster at specific Y values near predicted helix crossings)

These non-zero residuals are the actual "signal" — the points where the helix
model *should* predict crossings but the data shows asymmetric strand positions.
The "QFPIL gap" is the systematic offset between these actual-strand-asymmetry
Y values and the predicted helix crossings.

Direct test: sweep an additive offset c in the crossing formula
    y_cross(n) = ((2n+1)*pi/2 - phi) / k + c
and for each c, compute the mean distance between predicted crossings and the
20 non-zero-residual Y values. The optimum c is the calibration bias.

This is a DIRECT test (not via ζ-zero mapping). The mapping hypothesis is
loaded separately as a downstream check.

Author: Kanopi ✋
Date:   2026-06-07 04:50 CDT (Sunday, focused-window hour)
"""

import json
import numpy as np

WORKSPACE = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase10_DualStrandToZeta"
TRAJ_FILE = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/traj_vmPhase9.csv"


# ============================================================
# DATA LOAD
# ============================================================

def load_residual_data():
    """Load residuals (midpoint - 0.5) for each paired strand point."""
    import csv

    with open(TRAJ_FILE) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = []
        for row in reader:
            if len(row) >= 4 and row[2] and row[3]:
                try:
                    rows.append((int(row[0]), int(row[1]), float(row[2]), float(row[3])))
                except Exception:
                    pass

    frames = {}
    for traj_id, frame, x, y in rows:
        frames.setdefault(frame, []).append((x, y))

    all_residuals, all_ys = [], []
    for frame in sorted(frames.keys()):
        pts = sorted(frames[frame], key=lambda p: p[1])
        i = 0
        while i < len(pts) - 1:
            y0, y1 = pts[i][1], pts[i+1][1]
            if abs(y1 - y0) < 0.1:
                x0, x1 = pts[i][0], pts[i+1][0]
                if x0 > x1:
                    x0, x1 = x1, x0
                mid_x = (x0 + x1) / 2.0
                all_residuals.append(mid_x - 0.5)
                all_ys.append(y0)
                i += 2
            else:
                i += 1

    return np.array(all_ys), np.array(all_residuals)


def find_non_zero_crossing_ys(ys, residuals):
    """Find unique Y values where residual is non-zero (the 'signal' points).

    The non-zero residuals cluster at specific Y values. Each unique Y
    represents a candidate 'crossing-like' event in the data.
    """
    nz_mask = residuals != 0
    nz_ys = ys[nz_mask]
    unique_nz_ys = np.unique(nz_ys)
    return unique_nz_ys


def find_predicted_crossings(k, phi, c, y_min, y_max, n_max=200):
    """All crossings of cos(k*y + phi) = 0 with additive offset c."""
    crossings = []
    for n in range(n_max):
        y_cross = (np.pi / 2 - phi + n * np.pi) / k + c
        if y_cross > y_max:
            break
        if y_cross >= y_min:
            crossings.append(y_cross)
    return np.array(crossings)


def mean_nearest_distance(predicted, actual_ys):
    """For each actual Y, find nearest predicted crossing. Mean abs distance."""
    if len(predicted) == 0 or len(actual_ys) == 0:
        return float("inf"), []
    deltas = []
    for y_act in actual_ys:
        nearest = predicted[np.argmin(np.abs(predicted - y_act))]
        deltas.append(y_act - nearest)
    deltas = np.array(deltas)
    return float(np.mean(np.abs(deltas))), deltas.tolist()


# ============================================================
# MAIN
# ============================================================

def main():
    print("="*60)
    print("Phase 10 Re-run v3 — Additive Offset c Sweep")
    print("(direct residual-crossing match)")
    print("="*60)

    print("\nLoading data...")
    ys, residuals = load_residual_data()
    n_zero = int(np.sum(residuals == 0))
    n_nonzero = int(np.sum(residuals != 0))
    print(f"  Total points: {len(residuals)}")
    print(f"  Zero residuals (midpoint=0.5 exactly): {n_zero} ({100*n_zero/len(residuals):.1f}%)")
    print(f"  Non-zero residuals: {n_nonzero} ({100*n_nonzero/len(residuals):.1f}%)")

    # Identify the "signal" points — unique Y values with non-zero residual
    signal_ys = find_non_zero_crossing_ys(ys, residuals)
    print(f"  Unique Y values with non-zero residual: {len(signal_ys)}")
    print(f"  Signal Y values: {sorted(signal_ys.tolist())}")

    # Load ζ-zeros and fitted helix from spectral_analysis
    with open(f"{WORKSPACE}/python_code/spectral_analysis_results.json") as f:
        spec = json.load(f)
    k_use = spec["k_fitted"]
    phi_use = spec["phase"]
    A_use = spec["A"]
    zeta_gammas = spec["zeta_gammas"]
    print(f"\n  Loaded helix: k={k_use:.6f}, A={A_use:.6f}, phi={phi_use:.6f}")
    print(f"  Loaded {len(zeta_gammas)} ζ-zeros")

    # Baseline (c=0) — predicted crossings vs signal Ys
    print(f"\n--- Baseline (c=0) ---")
    pred_cross_base = find_predicted_crossings(k_use, phi_use, 0.0, ys.min(), ys.max())
    base_mean, base_deltas = mean_nearest_distance(pred_cross_base, signal_ys)
    print(f"  Predicted crossings in Y range: {len(pred_cross_base)}")
    print(f"  Mean |signal_y - nearest predicted|: {base_mean:.6f}")
    print(f"  First 5 deltas (signal - predicted): {base_deltas[:5]}")
    print(f"  All deltas (signal - predicted): {[f'{d:+.4f}' for d in base_deltas]}")

    # Sweep c
    print(f"\n--- Sweep c in [0.0, 0.05] (51 points) ---")
    c_values = np.linspace(0.0, 0.05, 51)
    sweep = []
    for c in c_values:
        pred_cross = find_predicted_crossings(k_use, phi_use, float(c), ys.min(), ys.max())
        mean_d, deltas = mean_nearest_distance(pred_cross, signal_ys)
        sweep.append({
            "c": float(c),
            "n_predicted_crossings": len(pred_cross),
            "mean_abs_delta": mean_d,
            "deltas": deltas,
        })

    # Optimum
    finite_sweep = [s for s in sweep if s["mean_abs_delta"] != float("inf")]
    if not finite_sweep:
        print("  All c values produced inf — aborting.")
        return

    best = min(finite_sweep, key=lambda r: r["mean_abs_delta"])
    best_c = best["c"]
    best_mean = best["mean_abs_delta"]

    # c closest to 0.02 (the predicted value)
    c_002_idx = int(np.argmin(np.abs(c_values - 0.02)))
    c_002 = sweep[c_002_idx]

    # Also sweep a WIDER range to see if the optimum is outside [0, 0.05]
    print(f"\n--- Extended sweep c in [-0.5, 1.0] (151 points) — to find global optimum ---")
    c_ext = np.linspace(-0.5, 1.0, 151)
    ext_sweep = []
    for c in c_ext:
        pred_cross = find_predicted_crossings(k_use, phi_use, float(c), ys.min(), ys.max())
        mean_d, deltas = mean_nearest_distance(pred_cross, signal_ys)
        ext_sweep.append({
            "c": float(c),
            "mean_abs_delta": mean_d,
        })

    ext_best = min([s for s in ext_sweep if s["mean_abs_delta"] != float("inf")],
                    key=lambda r: r["mean_abs_delta"])
    ext_best_c = ext_best["c"]
    ext_best_mean = ext_best["mean_abs_delta"]

    # Reduction vs baseline
    if base_mean > 0 and best_mean != float("inf"):
        reduction_pct = 100.0 * (base_mean - best_mean) / base_mean
    else:
        reduction_pct = 0.0

    # Verdict
    if ext_best_mean < 0.1:
        verdict = "CONFIRMED"
        verdict_note = (
            f"Optimum c={ext_best_c:.4f} (extended sweep) gives mean |delta| = {ext_best_mean:.4f}. "
            f"QFPIL gap closes cleanly with this offset."
        )
    elif ext_best_mean < base_mean * 0.5:
        verdict = "PARTIAL"
        verdict_note = (
            f"Optimum c={ext_best_c:.4f} reduces mean |delta| from {base_mean:.4f} to {ext_best_mean:.4f} "
            f"({reduction_pct:.1f}% reduction). Offset is real but the residual gap is non-trivial."
        )
    elif reduction_pct > 5:
        verdict = "WEAK"
        verdict_note = (
            f"Optimum c={ext_best_c:.4f} reduces mean |delta| by only {reduction_pct:.1f}%. "
            f"Offset is not the dominant source of the QFPIL gap."
        )
    else:
        verdict = "FALSIFIED"
        verdict_note = (
            f"Optimum c={ext_best_c:.4f} reduces mean |delta| by only {reduction_pct:.1f}%. "
            f"Additive offset does not close the QFPIL gap."
        )

    # Print summary
    print(f"\n{'='*60}")
    print(f"RESULT")
    print(f"{'='*60}")
    print(f"  Baseline (c=0) mean |delta|:           {base_mean:.6f}")
    print(f"  c=0.02 (predicted) mean |delta|:       {c_002['mean_abs_delta']:.6f}")
    print(f"  Optimum c in [0,0.05]:                 {best_c:.4f}, mean |delta|={best_mean:.6f}")
    print(f"  Optimum c in [-0.5, 1.0] (extended):  {ext_best_c:.4f}, mean |delta|={ext_best_mean:.6f}")
    print(f"  Reduction vs baseline:                 {reduction_pct:.1f}%")
    print(f"  Verdict:                               {verdict}")
    print(f"  Note:                                  {verdict_note}")

    # Save results
    output = {
        "task": "QFPIL Gap — Additive Offset c Sweep (v3, direct match)",
        "date": "2026-06-07 04:50 CDT",
        "author": "Kanopi ✋",
        "version_note": (
            "v1: deduped residuals → 0% signal (collapsed transitions). "
            "v2: raw residuals, helix fit → 0% variance explained, 0 crossings in train region. "
            "v3 (this): direct test — match signal Y values (non-zero residual locations) "
            "to predicted helix crossings, sweep c."
        ),
        "hypothesis": (
            "A small constant additive offset c in the helix-crossing formula "
            "y_cross(n) = ((2n+1)*pi/2 - phi) / k + c closes the QFPIL gap "
            "(systematic offset between actual strand-asymmetry Y values and "
            "predicted helix crossings)."
        ),
        "predicted_optimum": "c ≈ 0.02",
        "data_profile": {
            "n_total": int(len(residuals)),
            "n_zero_residuals": n_zero,
            "n_nonzero_residuals": n_nonzero,
            "n_unique_signal_ys": int(len(signal_ys)),
            "signal_ys": [float(y) for y in sorted(signal_ys.tolist())],
            "residual_std": float(residuals.std()),
        },
        "helix_params_used": {
            "k": float(k_use),
            "A": float(A_use),
            "phi": float(phi_use),
            "source": "spectral_analysis_results.json",
        },
        "baseline": {
            "c": 0.0,
            "n_predicted_crossings": int(len(pred_cross_base)),
            "mean_abs_delta": float(base_mean),
            "all_deltas": [float(d) for d in base_deltas],
        },
        "sweep_narrow": sweep,
        "c_at_predicted_optimum": {
            "c": 0.02,
            "mean_abs_delta": float(c_002["mean_abs_delta"]),
        },
        "optimum_narrow": {
            "c": float(best_c),
            "mean_abs_delta": float(best_mean),
        },
        "sweep_extended": ext_sweep,
        "optimum_extended": {
            "c": float(ext_best_c),
            "mean_abs_delta": float(ext_best_mean),
        },
        "reduction_pct_vs_baseline": float(reduction_pct),
        "verdict": verdict,
        "verdict_note": verdict_note,
        "interpretation": {
            "CONFIRMED": (
                "The QFPIL gap is a clean calibration-bias artifact of the form "
                "y_cross + c. The varying-k model + DFT residual 0.107 framework is preserved; "
                "4-line README addendum (with c value) is sufficient documentation."
            ),
            "PARTIAL": (
                "An offset helps but is not the whole story. The 0.02 gap is partly a "
                "calibration artifact and partly a real second-order effect. "
                "Recommend deeper investigation before handoff."
            ),
            "WEAK": (
                "An offset is not the dominant source. The gap is structural "
                "(k local, A local, or non-constant offset). "
                "Next move: revisit k(y) varying model and check whether the gap tracks ln(y+1)."
            ),
            "FALSIFIED": (
                "The gap is not a simple additive offset. The bias is structural. "
                "Next move: revisit the k(y) varying model and check whether the gap tracks ln(y+1)."
            ),
        }[verdict],
    }

    out_path = f"{WORKSPACE}/python_code/offset_sweep_c_results.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")

    return output


if __name__ == "__main__":
    main()
