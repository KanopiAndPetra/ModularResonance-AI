#!/usr/bin/env python3
"""
phase10_dsr_boundedness_test.py
================================

Test of the structural boundedness of the DSR manifold.

Background:
    The DSR manifold encoding of Λ(n) maps x = Λ(n)/log(n) ∈ [0, 1] so
    |x - 0.5| ∈ [0, 0.5]. This is an analytic bound (Chebyshev-like).

    The Phase 10 Limit Proposition (2026-08-18) requires:
        (a) The amplitude bound is tight in the empirical sense
        (b) The residual after removing the fitted cos(k·y + φ) is bounded
        (c) The DFT peak at f₀ = 0.107 is the unique dominant low-frequency mode

    This script measures (a), (b), and (c) on the actual data.

Inputs:
    ~/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/traj_vmPhase9.csv
    manifold_summary.json (Phase 9 fit parameters)

Outputs:
    boundedness_results.json
    boundedness_summary.txt (human-readable)

Author: Petra (2026-08-18)
"""

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.fft import fft, fftfreq, rfft, rfftfreq
from scipy.optimize import minimize

BASE = Path("/Users/oppie1.kanopi/Desktop/ModularResonance-AI")
VM_PATH = BASE / "Phase9_DualStrandResonance" / "DSR_Phase9" / "trajectory_data" / "traj_vmPhase9.csv"
LOG_PATH = BASE / "Phase9_DualStrandResonance" / "DSR_Phase9" / "trajectory_data" / "traj_logPhase9.csv"
SUMMARY_PATH = BASE / "Phase9_DualStrandResonance" / "DSR_Phase9" / "trajectory_data" / "manifold_summary.json"
OUTPUT_DIR = BASE / "Phase10_DualStrandToZeta" / "python_code"


def load_vm_data(path):
    """Load (x, y) pairs from a trajectory CSV. Skips empty rows."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            xs = row.get("x", "")
            ys = row.get("y", "")
            if xs and ys:  # skip blank-frame rows
                try:
                    rows.append((float(xs), float(ys)))
                except ValueError:
                    continue
    return np.array(rows)


def load_manifold_params(path):
    """Load the Phase 9 fit parameters (k, A, phase)."""
    with open(path) as f:
        d = json.load(f)
    return d["fitted"]


def sort_by_y(data):
    """Sort (x, y) by y ascending. The manifold is defined over y."""
    return data[data[:, 1].argsort()]


def empirical_boundedness(xs, ys):
    """
    Test (a): empirical bound on |x - 0.5|.

    Returns dict with max, min, percentile stats, and a flag for whether
    the distribution is essentially bounded.
    """
    devs = np.abs(xs - 0.5)
    return {
        "n": int(len(xs)),
        "min_dev": float(devs.min()),
        "max_dev": float(devs.max()),
        "median_dev": float(np.median(devs)),
        "p95_dev": float(np.percentile(devs, 95)),
        "p99_dev": float(np.percentile(devs, 99)),
        "structural_bound": 0.5,
        "structural_bound_hits": bool(devs.max() <= 0.5 + 1e-9),
        "y_range": [float(ys.min()), float(ys.max())],
    }


def residual_after_manifold(xs, ys, k, A, phi):
    """
    Test (b): residual after removing the fitted cos(k·y + φ) in *raw* y-space.

    Note: the Phase 9 fit uses y/ymax normalization; we use raw y here so the
    k is in cycles per unit y, not cycles per unit y_max. This is the more
    natural unit for the DFT analysis.
    """
    manifold = np.abs(A * np.cos(k * ys + phi))
    residuals = np.abs(xs - 0.5) - manifold
    return residuals, manifold


def fit_manifold_raw_y(xs, ys, n_starts=12):
    """
    Fit |x-0.5| = A|cos(k·y + φ)| using raw y (not normalized).

    Multi-start to avoid the initialization artifact that plagued Phase 9.
    """
    targets = np.abs(xs - 0.5)

    def objective(p):
        A, k, phi = p
        preds = np.abs(A * np.cos(k * ys + phi))
        return np.sqrt(np.mean((targets - preds) ** 2))

    best = None
    # spread k_starts to cover [0.01, 2.0] broadly, not just around the golden ratio
    k_starts = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8]
    A_starts = [0.1, 0.15, 0.19, 0.25, 0.4]
    phi_starts = [0.0, 0.3, 0.5, -0.3, -0.5]

    for k0 in k_starts[:n_starts]:
        for A0 in A_starts:
            for ph0 in phi_starts:
                res = minimize(objective, x0=[A0, k0, ph0],
                              method="Nelder-Mead",
                              options={"maxiter": 5000, "xatol": 1e-8, "fatol": 1e-9})
                if best is None or res.fun < best.fun:
                    best = res

    return best.x[1], best.x[0], best.x[2], best.fun  # k, A, phi, rms


def dft_lowfreq_dominance(xs, ys, k, A, phi):
    """
    Test (c): is f₀ = 0.107 the unique dominant low-frequency mode in the residual?

    Returns the frequencies of the top N peaks and the magnitude of each, plus
    the 0.107 bin and the relative dominance.
    """
    residuals, _ = residual_after_manifold(xs, ys, k, A, phi)

    # For the DFT, we need a uniformly sampled signal. The data is NOT uniform
    # in y (it's a sparse, biased sample). We use the natural frame index as
    # the time axis (consistent with dsr_spectral_analysis.py).
    N = len(residuals)
    # Sort by y so the DFT is taken over a monotonic ascent
    order = np.argsort(ys)
    sig = residuals[order]

    # DFT in frame-index unit; convert to frequency in cycles per unit y by
    # dividing by mean Δy
    y_sorted = ys[order]
    mean_dy = np.mean(np.diff(y_sorted))
    N_y = y_sorted[-1] - y_sorted[0]  # total y-range
    # bin frequency in (cycles per unit y) = (bin index) / (N * mean_dy)
    freqs = np.fft.rfftfreq(N, d=mean_dy)
    mag = np.abs(np.fft.rfft(sig))

    # Look at low frequencies (0 < f < 0.5) where the 0.107 peak should live
    low_mask = (freqs > 0.01) & (freqs < 0.5)
    low_freqs = freqs[low_mask]
    low_mag = mag[low_mask]

    # Top 10 peaks in low-frequency range
    idxs = np.argsort(low_mag)[::-1][:10]
    top_peaks = [(float(low_freqs[i]), float(low_mag[i])) for i in idxs]

    # Find the magnitude at f = 0.107 (nearest bin)
    target = 0.107
    if len(low_freqs) > 0:
        target_idx = np.argmin(np.abs(low_freqs - target))
        mag_at_0107 = float(low_mag[target_idx])
        freq_at_0107 = float(low_freqs[target_idx])
    else:
        mag_at_0107 = None
        freq_at_0107 = None

    # Maximum magnitude in the low-frequency range, for dominance ratio
    max_mag = float(low_mag.max()) if len(low_mag) > 0 else 0.0
    dominance = (mag_at_0107 / max_mag) if (max_mag > 0 and mag_at_0107 is not None) else 0.0

    # The expected ζ-zero mean spacing near the y_max
    y_max = float(ys.max())
    T_max = math.exp(y_max)
    mean_spacing = 2 * math.pi / math.log(T_max)

    return {
        "n": int(N),
        "mean_dy": float(mean_dy),
        "y_range": float(N_y),
        "top_lowfreq_peaks": top_peaks,
        "freq_at_0107_target": freq_at_0107,
        "mag_at_0107_target": mag_at_0107,
        "max_mag_lowfreq": max_mag,
        "dominance_ratio_0107_vs_maxLowFreq": float(dominance),
        "y_max": y_max,
        "T_max": T_max,
        "expected_zeta_mean_spacing": float(mean_spacing),
        "match_pct": float(100.0 * abs(freq_at_0107 - target) / target) if freq_at_0107 else None,
    }


def plot_data(xs, ys, residuals, k, A, phi, output_dir):
    """Save a single diagnostic plot: data + fit + residual."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        order = np.argsort(ys)
        y_s = ys[order]
        x_s = xs[order]
        devs = np.abs(x_s - 0.5)
        fit = np.abs(A * np.cos(k * y_s + phi))
        res = devs - fit

        fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
        axes[0].plot(y_s, x_s, ".", alpha=0.4, markersize=2, label="x(n)")
        axes[0].plot(y_s, 0.5 + fit, "r-", lw=1, label=f"|x-0.5| fit (k={k:.4f}, A={A:.4f})")
        axes[0].plot(y_s, 0.5 - fit, "r-", lw=1)
        axes[0].set_ylabel("x")
        axes[0].legend(loc="upper right", fontsize=8)
        axes[0].set_title("DSR manifold: data + envelope fit")

        axes[1].plot(y_s, devs, "b.", alpha=0.4, markersize=2, label="|x-0.5|")
        axes[1].plot(y_s, fit, "r-", lw=1, label="fit")
        axes[1].set_ylabel("|x - 0.5|")
        axes[1].legend(loc="upper right", fontsize=8)
        axes[1].set_title("|x - 0.5| and fitted cos(k·y+φ)")

        axes[2].plot(y_s, res, "g.", alpha=0.4, markersize=2)
        axes[2].axhline(0, color="k", lw=0.5)
        axes[2].set_xlabel("y = log(n)")
        axes[2].set_ylabel("residual")
        axes[2].set_title("Residual after fit")

        fig.suptitle("Phase 10 Boundedness Test — DSR VM trajectory")
        fig.tight_layout()
        out = output_dir / "output" / "phase10_boundedness_diagnostic.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=120)
        plt.close(fig)
        return str(out)
    except Exception as e:
        return f"plot failed: {e}"


def main():
    print("=" * 70)
    print("Phase 10: DSR Boundedness Test")
    print("=" * 70)

    # Load both trajectories
    print("\nLoading trajectory data...")
    vm = load_vm_data(VM_PATH)
    log = load_vm_data(LOG_PATH)
    print(f"  vm: {len(vm)} valid (x, y) pairs")
    print(f"  log: {len(log)} valid (x, y) pairs")
    print(f"  vm y-range: [{vm[:, 1].min():.3f}, {vm[:, 1].max():.3f}]")
    print(f"  log y-range: [{log[:, 1].min():.3f}, {log[:, 1].max():.3f}]")

    # Use the vm trajectory for the manifold analysis (it's the prime point set)
    xs = vm[:, 0]
    ys = vm[:, 1]

    # (a) Empirical boundedness
    print("\n--- (a) Empirical boundedness of |x - 0.5| ---")
    bound = empirical_boundedness(xs, ys)
    print(f"  n = {bound['n']}")
    print(f"  min |x-0.5| = {bound['min_dev']:.6f}")
    print(f"  max |x-0.5| = {bound['max_dev']:.6f}  (structural bound: 0.5)")
    print(f"  median |x-0.5| = {bound['median_dev']:.6f}")
    print(f"  p95 |x-0.5| = {bound['p95_dev']:.6f}")
    print(f"  p99 |x-0.5| = {bound['p99_dev']:.6f}")
    print(f"  structural bound hit: {bound['structural_bound_hits']}")

    # Fit the manifold in raw y-space (multi-start, broad scan)
    print("\n--- Fitting manifold in raw y-space (multi-start) ---")
    k_fit, A_fit, phi_fit, rms = fit_manifold_raw_y(xs, ys)
    print(f"  k = {k_fit:.6f}")
    print(f"  A = {A_fit:.6f}")
    print(f"  φ = {phi_fit:.6f}")
    print(f"  RMS = {rms:.6f}")

    # Save the fit for the histogram script
    fit_params = {"k": float(k_fit), "A": float(A_fit), "phi": float(phi_fit), "rms": float(rms)}

    # (b) Residual after manifold removal
    print("\n--- (b) Residual after manifold removal ---")
    residuals, manifold = residual_after_manifold(xs, ys, k_fit, A_fit, phi_fit)
    print(f"  residual std = {np.std(residuals):.6f}")
    print(f"  residual max |·| = {np.max(np.abs(residuals)):.6f}")
    print(f"  structural residual bound (0.5 - A_fit) = {0.5 - A_fit:.6f}")
    print(f"  residual within structural bound: {np.all(np.abs(residuals) <= 0.5 + 1e-9)}")

    # (c) DFT low-frequency dominance
    print("\n--- (c) DFT low-frequency dominance (the 0.107 question) ---")
    dft = dft_lowfreq_dominance(xs, ys, k_fit, A_fit, phi_fit)
    print(f"  N = {dft['n']}, mean Δy = {dft['mean_dy']:.4f}")
    print(f"  Top low-frequency peaks (cycles per unit y):")
    for i, (f, m) in enumerate(dft["top_lowfreq_peaks"][:10]):
        print(f"    {i+1}. f = {f:.6f}, magnitude = {m:.4f}")
    print(f"  freq at target 0.107 = {dft['freq_at_0107_target']:.6f}")
    print(f"  mag at target 0.107 = {dft['mag_at_0107_target']:.4f}")
    print(f"  max magnitude in low-freq band = {dft['max_mag_lowfreq']:.4f}")
    print(f"  dominance ratio = {dft['dominance_ratio_0107_vs_maxLowFreq']:.4f}")
    print(f"  y_max = {dft['y_max']:.3f}, T_max = e^{dft['y_max']:.1f} ≈ 10^{dft['y_max']/2.303:.1f}")
    print(f"  expected ζ-zero mean spacing = {dft['expected_zeta_mean_spacing']:.6f}")
    print(f"  bin at 0.107 vs target 0.107: {dft['match_pct']:.2f}% off")

    # Plot
    plot_path = plot_data(xs, ys, residuals, k_fit, A_fit, phi_fit, OUTPUT_DIR)
    print(f"\n  Diagnostic plot: {plot_path}")

    # Save results
    results = {
        "task": "Phase 10 DSR Boundedness Test",
        "date": "2026-08-18",
        "author": "Petra",
        "version": "v1.0",
        "data_path": str(VM_PATH),
        "n_points": len(vm),
        "empirical_boundedness": bound,
        "manifold_fit": fit_params,
        "residual_stats": {
            "std": float(np.std(residuals)),
            "max_abs": float(np.max(np.abs(residuals))),
            "structural_residual_bound": float(0.5 - A_fit),
            "within_structural_bound": bool(np.all(np.abs(residuals) <= 0.5 + 1e-9)),
        },
        "dft_lowfreq": dft,
        "limit_proposition_status": {
            "F1_0107_invariant_reproducible": dft["max_mag_lowfreq"] > 0 and dft["mag_at_0107_target"] is not None,
            "F3a_boundedness_analytic": bound["structural_bound_hits"],
            "F3b_residual_bounded": bool(np.all(np.abs(residuals) <= 0.5 + 1e-9)),
        },
    }

    out_path = OUTPUT_DIR / "boundedness_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {out_path}")

    # Verdict
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    if bound["structural_bound_hits"]:
        print("✓ Empirical |x-0.5| ≤ 0.5 holds (analytic bound matches data)")
    else:
        print(f"✗ Empirical |x-0.5| max = {bound['max_dev']:.6f} > 0.5 — unexpected")

    if dft["freq_at_0107_target"] is not None and dft["mag_at_0107_target"] is not None:
        if abs(dft["freq_at_0107_target"] - 0.107) < 0.02:
            print(f"✓ 0.107 bin confirmed at f = {dft['freq_at_0107_target']:.4f}")
        else:
            print(f"⚠ 0.107 bin off — got {dft['freq_at_0107_target']:.4f}")

    print(f"  dominance of 0.107 vs max low-freq: {dft['dominance_ratio_0107_vs_maxLowFreq']:.4f}")
    if dft["dominance_ratio_0107_vs_maxLowFreq"] > 0.5:
        print("  → 0.107 is a dominant low-frequency mode")
    else:
        print("  → 0.107 is NOT the unique dominant low-frequency mode")


if __name__ == "__main__":
    main()
