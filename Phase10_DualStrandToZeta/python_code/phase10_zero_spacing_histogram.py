#!/usr/bin/env python3
"""
phase10_zero_spacing_histogram.py
==================================

Test of LIMIT-1 from the Phase 10 limit proposition (2026-08-18).

Background:
    The Limit Proposition hypothesizes that the DSR manifold's zero crossings
    form a discrete sequence whose spacings reproduce the ζ-zero spacing
    distribution on the critical line (modulo an overall scale factor).

The empirical test:
    1. Compute the zero crossings of |x - 0.5| from the manifold data
    2. Compute the spacings Δy_c between successive crossings
    3. Compute the **scaled spacings** u_n = (Δy_c^(n) - mean) / std
    4. Compute the **pair correlation** of {u_n}
    5. Compare to the GUE/Montgomery form: 1 - (sin(πu) / (πu))²

If the empirical pair correlation matches GUE, the manifold's zero crossings
behave like ζ-zeros on the critical line. If it doesn't, the proposition
falsifies.

Inputs:
    ~/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/traj_vmPhase9.csv
    boundedness_results.json (from the boundedness test)

Outputs:
    zero_spacing_results.json
    zero_spacing_histogram.png
    zero_spacing_pair_correlation.png

Author: Petra (2026-08-18)
"""

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.fft import fft, fftfreq

BASE = Path("/Users/oppie1.kanopi/Desktop/ModularResonance-AI")
VM_PATH = BASE / "Phase9_DualStrandResonance" / "DSR_Phase9" / "trajectory_data" / "traj_vmPhase9.csv"
LOG_PATH = BASE / "Phase9_DualStrandResonance" / "DSR_Phase9" / "trajectory_data" / "traj_logPhase9.csv"
BOUNDEDNESS_PATH = BASE / "Phase10_DualStrandToZeta" / "python_code" / "boundedness_results.json"
OUTPUT_DIR = BASE / "Phase10_DualStrandToZeta" / "python_code"


def load_vm_data(path):
    """Load (x, y) pairs from a trajectory CSV. Skips empty rows."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            xs = row.get("x", "")
            ys = row.get("y", "")
            if xs and ys:
                try:
                    rows.append((float(xs), float(ys)))
                except ValueError:
                    continue
    return np.array(rows)


def find_zero_crossings(xs, ys, x_threshold_low=0.45, x_threshold_high=0.55):
    """
    Find zero crossings of |x - 0.5|, i.e., points where x transitions
    from below 0.5 to above 0.5 (or vice versa).

    A simple approach: sort by y, then look for points where x crosses 0.5.
    Linear interpolation gives the y value of the crossing.

    The thresholds (0.45, 0.55) allow rounding: a point with x = 0.499 is
    "essentially on the critical line" but we only count crossings that
    actually move across 0.5.
    """
    order = np.argsort(ys)
    y_s = ys[order]
    x_s = xs[order]

    crossings = []
    for i in range(len(x_s) - 1):
        x1, x2 = x_s[i], x_s[i + 1]
        y1, y2 = y_s[i], y_s[i + 1]
        # Crossing of x = 0.5
        if (x1 - 0.5) * (x2 - 0.5) < 0:
            # Linear interpolation
            t = (0.5 - x1) / (x2 - x1)
            yc = y1 + t * (y2 - y1)
            crossings.append(yc)

    # Deduplicate crossings that are too close (numerical noise)
    if len(crossings) > 1:
        dedup = [crossings[0]]
        for c in crossings[1:]:
            if c - dedup[-1] > 0.01:  # minimum Δy of 0.01
                dedup.append(c)
        crossings = dedup

    return np.array(crossings)


def spacing_distribution(crossings):
    """Compute Δy_c and the scaled (normalized) spacings."""
    if len(crossings) < 2:
        return None
    deltas = np.diff(crossings)
    mean_d = np.mean(deltas)
    std_d = np.std(deltas)
    if std_d == 0:
        return None
    scaled = (deltas - mean_d) / std_d
    return {
        "n_crossings": int(len(crossings)),
        "n_spacings": int(len(deltas)),
        "raw_mean": float(mean_d),
        "raw_std": float(std_d),
        "raw_min": float(deltas.min()),
        "raw_max": float(deltas.max()),
        "raw_median": float(np.median(deltas)),
        "scaled_mean": float(np.mean(scaled)),
        "scaled_std": float(np.std(scaled)),
        "scaled_deltas": scaled.tolist(),
        "raw_deltas": deltas.tolist(),
    }


def pair_correlation(deltas, n_bins=50, max_u=5.0):
    """
    Compute the empirical pair correlation function.

    Pair correlation R(u) = (1/N) Σ_{i≠j} δ(u - (s_i - s_j))

    where s_i are the scaled spacings. We compute this on a uniform grid
    and compare to the GUE form: R_GUE(u) = 1 - (sin(πu) / (πu))².
    """
    if len(deltas) < 4:
        return None

    s = deltas
    n = len(s)

    # All pair differences (avoid double-counting and self-pairs)
    diffs = []
    for i in range(n):
        for j in range(i + 1, n):
            diffs.append(s[j] - s[i])
    diffs = np.array(diffs)

    # Bin into the empirical pair correlation
    bins = np.linspace(-max_u, max_u, 2 * n_bins + 1)
    hist, edges = np.histogram(diffs, bins=bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2

    return {
        "n_pairs": len(diffs),
        "bin_centers": centers.tolist(),
        "hist_values": hist.tolist(),
        "max_u": max_u,
        "n_bins": n_bins,
    }


def gue_pair_correlation(u):
    """GUE / Montgomery pair correlation: 1 - (sin(πu)/(πu))² ."""
    u = np.asarray(u, dtype=float)
    # Avoid division by zero at u=0
    result = np.ones_like(u)
    mask = u != 0
    result[mask] = 1 - (np.sin(np.pi * u[mask]) / (np.pi * u[mask])) ** 2
    return result


def poisson_pair_correlation(u):
    """Poisson pair correlation: 1 (no level repulsion)."""
    return np.ones_like(u, dtype=float)


def compare_pair_correlations(empirical_centers, empirical_values):
    """Compare empirical pair correlation to GUE and Poisson."""
    u = np.array(empirical_centers)
    gue = gue_pair_correlation(u)
    poisson = poisson_pair_correlation(u)

    # Use only the central region (-2 < u < 2) for the comparison (where
    # the GUE vs Poisson distinction is strongest)
    mask = (u > -2) & (u < 2) & (u != 0)
    emp = np.array(empirical_values)[mask]
    if len(emp) == 0:
        return None, None, None

    # GUE has a strong DIP at u=0 (level repulsion). Poisson does not.
    # Compute the central-bin value (closest to u=0)
    central_idx = np.argmin(np.abs(u))
    emp_central = empirical_values[central_idx]
    gue_central = gue[central_idx]
    poisson_central = poisson[central_idx]

    # Compute L2 distance from GUE and from Poisson in the central region
    gue_l2 = float(np.sqrt(np.mean((emp - gue[mask]) ** 2)))
    poisson_l2 = float(np.sqrt(np.mean((emp - poisson[mask]) ** 2)))

    # Central-bin test: is the empirical value closer to GUE or Poisson?
    gue_distance = abs(emp_central - gue_central)
    poisson_distance = abs(emp_central - poisson_central)

    return {
        "central_bin_u": float(u[central_idx]),
        "central_bin_empirical": float(emp_central),
        "central_bin_GUE": float(gue_central),
        "central_bin_Poisson": float(poisson_central),
        "distance_to_GUE": float(gue_distance),
        "distance_to_Poisson": float(poisson_distance),
        "L2_distance_to_GUE": gue_l2,
        "L2_distance_to_Poisson": poisson_l2,
        "preferred_model": "GUE" if gue_distance < poisson_distance else "Poisson",
    }


def plot_histogram(deltas, output_path):
    """Histogram of raw spacings."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        # Raw spacing histogram
        axes[0].hist(deltas, bins=30, color="steelblue", alpha=0.7, edgecolor="k")
        axes[0].axvline(np.mean(deltas), color="red", linestyle="--",
                        label=f"mean = {np.mean(deltas):.3f}")
        axes[0].set_xlabel("Δy_c (raw spacing)")
        axes[0].set_ylabel("count")
        axes[0].set_title(f"Zero-crossing spacings, n={len(deltas)}")
        axes[0].legend()

        # Scaled spacings
        scaled = (deltas - np.mean(deltas)) / np.std(deltas)
        axes[1].hist(scaled, bins=30, color="seagreen", alpha=0.7, edgecolor="k")
        axes[1].set_xlabel("Δy_c (scaled)")
        axes[1].set_ylabel("count")
        axes[1].set_title("Scaled spacings (mean 0, std 1)")

        fig.suptitle("Phase 10: DSR manifold zero-crossing spacing distribution")
        fig.tight_layout()
        fig.savefig(output_path, dpi=120)
        plt.close(fig)
        return True
    except Exception as e:
        return f"histogram plot failed: {e}"


def plot_pair_correlation(empirical_centers, empirical_values, output_path):
    """Pair correlation plot: empirical vs GUE vs Poisson."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        u = np.array(empirical_centers)
        emp = np.array(empirical_values)
        gue = gue_pair_correlation(u)
        poisson = poisson_pair_correlation(u)

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(u, emp, "o-", color="steelblue", alpha=0.7, label="empirical", markersize=4)
        ax.plot(u, gue, "r--", lw=2, label="GUE / Montgomery: 1 - (sin(πu)/(πu))²")
        ax.plot(u, poisson, "g--", lw=2, label="Poisson: 1 (no level repulsion)")
        ax.axhline(1, color="k", lw=0.5, alpha=0.5)
        ax.set_xlabel("scaled spacing difference u")
        ax.set_ylabel("pair correlation R(u)")
        ax.set_title("Phase 10: Pair correlation of DSR zero-crossing spacings")
        ax.set_xlim(-5, 5)
        ax.set_ylim(-0.5, 2.5)
        ax.legend(loc="upper right")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(output_path, dpi=120)
        plt.close(fig)
        return True
    except Exception as e:
        return f"pair correlation plot failed: {e}"


def main():
    print("=" * 70)
    print("Phase 10: Zero-Crossing Spacing Distribution")
    print("=" * 70)

    # Load manifold data
    print("\nLoading data...")
    vm = load_vm_data(VM_PATH)
    log = load_log_data = load_vm_data(LOG_PATH)
    print(f"  vm: {len(vm)} points, y range [{vm[:, 1].min():.3f}, {vm[:, 1].max():.3f}]")
    print(f"  log: {len(log)} points, y range [{log[:, 1].min():.3f}, {log[:, 1].max():.3f}]")

    # Load boundedness results to use the fitted parameters
    with open(BOUNDEDNESS_PATH) as f:
        boundedness = json.load(f)
    k_fit = boundedness["manifold_fit"]["k"]
    A_fit = boundedness["manifold_fit"]["A"]
    phi_fit = boundedness["manifold_fit"]["phi"]
    print(f"\nUsing boundedness fit: k={k_fit:.6f}, A={A_fit:.6f}, φ={phi_fit:.6f}")

    # Test 1: VM trajectory zero crossings
    print("\n--- Test 1: VM trajectory crossings ---")
    vm_crossings = find_zero_crossings(vm[:, 0], vm[:, 1])
    print(f"  crossings: {len(vm_crossings)}")
    if len(vm_crossings) > 1:
        vm_delta = np.diff(vm_crossings)
        print(f"  spacings: mean={np.mean(vm_delta):.4f}, std={np.std(vm_delta):.4f}")
        print(f"  min={np.min(vm_delta):.4f}, max={np.max(vm_delta):.4f}")
        print(f"  median={np.median(vm_delta):.4f}")

    # Test 2: LOG trajectory zero crossings
    print("\n--- Test 2: LOG trajectory crossings ---")
    log_crossings = find_zero_crossings(log[:, 0], log[:, 1])
    print(f"  crossings: {len(log_crossings)}")
    if len(log_crossings) > 1:
        log_delta = np.diff(log_crossings)
        print(f"  spacings: mean={np.mean(log_delta):.4f}, std={np.std(log_delta):.4f}")
        print(f"  min={np.min(log_delta):.4f}, max={np.max(log_delta):.4f}")
        print(f"  median={np.median(log_delta):.4f}")

    # Combined crossings (union of both trajectories)
    print("\n--- Test 3: Combined VM + LOG crossings ---")
    combined = np.sort(np.concatenate([vm_crossings, log_crossings]))
    if len(combined) > 1:
        combined_delta = np.diff(combined)
        print(f"  combined crossings: {len(combined)}")
        print(f"  spacings: mean={np.mean(combined_delta):.4f}, std={np.std(combined_delta):.4f}")
        print(f"  min={np.min(combined_delta):.4f}, max={np.max(combined_delta):.4f}")

    # LIMIT-1 prediction
    # The slowest DFT mode in the boundedness test was f₀ ≈ 0.107
    # Predicted mean spacing: 1 / (k_fit * f₀) — but k_fit is in cycles/unit y
    # and f₀ is in cycles/unit y, so the product f₀ * k_fit is...wait, k_fit
    # is the *helix pitch*, and the zero crossing occurs at half-period of
    # the cos. The spacing between crossings is Δy = π / (k_fit * local_factor)
    # multiplied by some factor from the A envelope.
    # A simpler prediction: Δy_c ≈ π / k_fit / 2 if k_fit is the fundamental
    # frequency, OR Δy_c ≈ 1 / f₀ if the 0.107 is the dominant harmonic.
    #
    # The two predictions are not the same:
    #   Pred A (from helix fundamental):  π / (2 * k_fit) = π / (2 * 0.0255) ≈ 61.6
    #   Pred B (from 0.107 harmonic):   1 / f₀ ≈ 1 / 0.107 ≈ 9.35
    #
    # Pred A is in unit y, where 1 unit = log(n).
    # Pred B is in cycles per unit y, so the spacing in unit y is 1 / 0.107 ≈ 9.35.
    # The 9.35 prediction is much smaller than 61.6 — different physics.
    #
    # Empirical: the actual spacings will tell us which is right.

    print("\n--- LIMIT-1 predictions ---")
    helix_pitch_pred = math.pi / (2 * k_fit)
    harmonic_pred = 1 / 0.107
    print(f"  Pred A (from helix fundamental k_fit={k_fit:.4f}): Δy_c ≈ {helix_pitch_pred:.4f}")
    print(f"  Pred B (from 0.107 DFT residual): Δy_c ≈ {harmonic_pred:.4f}")
    print(f"  Empirical (vm): mean={np.mean(vm_delta):.4f}, std={np.std(vm_delta):.4f}")
    print(f"  Empirical (log): mean={np.mean(log_delta):.4f}, std={np.std(log_delta):.4f}")
    print(f"  Empirical (combined): mean={np.mean(combined_delta):.4f}, std={np.std(combined_delta):.4f}")

    # Pair correlation
    print("\n--- Pair correlation (combined dataset) ---")
    pc = pair_correlation(combined_delta, n_bins=50, max_u=5.0)
    if pc is None:
        print("  insufficient data for pair correlation")
        return
    print(f"  n_pairs: {pc['n_pairs']}")
    print(f"  bin range: [-{pc['max_u']}, {pc['max_u']}], n_bins: {pc['n_bins']}")

    cmp_result = compare_pair_correlations(pc["bin_centers"], pc["hist_values"])
    if cmp_result:
        print(f"  central bin u = {cmp_result['central_bin_u']:.4f}")
        print(f"  central empirical = {cmp_result['central_bin_empirical']:.4f}")
        print(f"  central GUE = {cmp_result['central_bin_GUE']:.4f}")
        print(f"  central Poisson = {cmp_result['central_bin_Poisson']:.4f}")
        print(f"  L2 distance to GUE = {cmp_result['L2_distance_to_GUE']:.4f}")
        print(f"  L2 distance to Poisson = {cmp_result['L2_distance_to_Poisson']:.4f}")
        print(f"  preferred model: {cmp_result['preferred_model']}")

    # Plots
    out_dir = OUTPUT_DIR / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_histogram(combined_delta, out_dir / "phase10_zero_spacing_histogram.png")
    plot_pair_correlation(pc["bin_centers"], pc["hist_values"],
                          out_dir / "phase10_zero_spacing_pair_correlation.png")
    print(f"\n  plots: {out_dir}/phase10_zero_spacing_*.png")

    # Save results
    results = {
        "task": "Phase 10 Zero-Crossing Spacing Distribution",
        "date": "2026-08-18",
        "author": "Petra",
        "version": "v1.0",
        "vm_crossings": {
            "n": int(len(vm_crossings)),
            "mean_spacing": float(np.mean(vm_delta)) if len(vm_crossings) > 1 else None,
            "std_spacing": float(np.std(vm_delta)) if len(vm_crossings) > 1 else None,
        },
        "log_crossings": {
            "n": int(len(log_crossings)),
            "mean_spacing": float(np.mean(log_delta)) if len(log_crossings) > 1 else None,
            "std_spacing": float(np.std(log_delta)) if len(log_crossings) > 1 else None,
        },
        "combined_crossings": {
            "n": int(len(combined)),
            "mean_spacing": float(np.mean(combined_delta)),
            "std_spacing": float(np.std(combined_delta)),
            "median_spacing": float(np.median(combined_delta)),
            "min_spacing": float(np.min(combined_delta)),
            "max_spacing": float(np.max(combined_delta)),
        },
        "limit_proposition_predictions": {
            "pred_A_helix_pitch": float(helix_pitch_pred),
            "pred_B_0107_harmonic": float(harmonic_pred),
            "p_fit": float(k_fit),
            "A_fit": float(A_fit),
            "phi_fit": float(phi_fit),
        },
        "pair_correlation_comparison": cmp_result,
        "limit_proposition_status": {
            "n_crossings_adequate": len(combined) >= 10,
            "favored_model": cmp_result["preferred_model"] if cmp_result else "indeterminate",
        },
    }

    out_path = OUTPUT_DIR / "zero_spacing_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {out_path}")

    # Verdict
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    pred_A_dist = abs(np.mean(combined_delta) - helix_pitch_pred) / helix_pitch_pred
    pred_B_dist = abs(np.mean(combined_delta) - harmonic_pred) / harmonic_pred
    print(f"  Mean empirical Δy_c = {np.mean(combined_delta):.4f}")
    print(f"  Pred A (helix pitch {helix_pitch_pred:.4f}): {pred_A_dist*100:.1f}% off")
    print(f"  Pred B (0.107 harmonic {harmonic_pred:.4f}): {pred_B_dist*100:.1f}% off")
    if cmp_result:
        print(f"  Pair correlation: {cmp_result['preferred_model']} preferred")
        if cmp_result["preferred_model"] == "GUE":
            print("  → consistent with critical-line ζ-zero distribution")
        else:
            print("  → NOT consistent with critical-line GUE; LIMIT-1 partially falsified")


if __name__ == "__main__":
    main()
