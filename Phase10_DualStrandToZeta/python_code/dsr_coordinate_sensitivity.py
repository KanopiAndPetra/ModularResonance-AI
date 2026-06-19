#!/usr/bin/env python3
"""
dsr_coordinate_sensitivity.py
==============================
Test whether the fitted k = 0.618413 ≈ 1/φ is invariant across different
coordinate parametrizations of the DSR manifold.

If k changes substantially → coordinate artifact
If k·log(N_max) remains constant → k may be real (geometric invariant)

References:
- ~/Desktop/ModularResonance-AI/Phase10_DualStrandToZeta/Phase10_PLANNING.md
- ~/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/
"""

import csv
import math
import numpy as np
from scipy.optimize import minimize

DATA_PATH = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/traj_vmPhase9.csv"
OUTPUT_PATH = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase10_DualStrandToZeta/python_code/coordinate_sensitivity_results.csv"


def load_vm_data():
    """Load von Mangoldt trajectory data."""
    rows = []
    with open(DATA_PATH) as f:
        reader = csv.DictReader(f)
        for row in reader:
            x = row.get('x', '')
            y = row.get('y', '')
            if x and y:
                rows.append((float(x), float(y)))
    return np.array(rows)


def fit_manifold(xs, ys, model_label):
    """
    Fit |x - 0.5| = A * |cos(k * y + phase)| to data.
    Returns (k, A, phase, RMS) or (None,...) on failure.
    """
    if len(xs) < 10:
        return None, None, None, None

    targets = np.abs(xs - 0.5)
    ys_norm = ys / ys.max()  # normalize y to [0, 1]

    def objective(params):
        A, k, phase = params
        preds = np.abs(A * np.cos(k * ys_norm + phase))
        return np.sum((targets - preds) ** 2)

    best = None
    for k0 in [0.5, 0.618, 0.7, 0.8]:
        for A0 in [0.15, 0.19, 0.25]:
            for phase0 in [0.0, 0.3, 0.6]:
                res = minimize(
                    objective,
                    x0=[A0, k0, phase0],
                    method='Nelder-Mead',
                    options={'maxiter': 2000}
                )
                if best is None or res.fun < best.fun:
                    best = res

    if best is None:
        return None, None, None, None

    A, k, phase = best.x
    preds = np.abs(A * np.cos(k * ys_norm + phase))
    rms = np.sqrt(np.mean((targets - preds) ** 2))
    return k, A, phase, rms


def run_sensitivity():
    """
    10 coordinate parametrizations for the von Mangoldt manifold.
    Each gives a different (k, A, phase) fit.
    """
    data = load_vm_data()
    xs_raw = data[:, 0]
    ys_raw = data[:, 1]

    N_max = ys_raw.max()
    log_N_max = math.log(N_max)

    print(f"Loaded {len(data)} points, N_max={N_max:.2f}, log(N_max)={log_N_max:.4f}")
    print(f"Golden ratio inverse: {1/((1+math.sqrt(5))/2):.6f}")
    print()

    parametrizations = [
        # Label, x_func, y_func
        ("baseline: |x-0.5| vs y/ymax", lambda x, y: x, lambda x, y: y/ y.max()),
        ("x/n vs log(n)/log(N)", lambda x, y: x/y if y > 0 else 0, lambda x, y: math.log(y)/math.log(N_max) if y > 1 else 0),
        ("log(x)/log(N_max) vs log(y)/log(N_max)", lambda x, y: math.log(max(abs(x),1e-10))/log_N_max if x != 0 else 0, lambda x, y: math.log(max(y,1))/log_N_max),
        ("x/sqrt(y) vs log(y)", lambda x, y: x/math.sqrt(y) if y > 0 else 0, lambda x, y: math.log(max(y,1))),
        ("x*log(y)/y vs log(y)", lambda x, y: x*math.log(max(y,2))/y if y > 0 else 0, lambda x, y: math.log(max(y,1))),
        ("sqrt(|x|) vs y^0.5", lambda x, y: math.sqrt(abs(x)), lambda x, y: math.sqrt(y)),
        ("x/y^0.3 vs log(y)", lambda x, y: x/(y**0.3) if y > 0 else 0, lambda x, y: math.log(max(y,1))),
        ("log(1+|x|) vs log(y)", lambda x, y: math.log(1+abs(x)), lambda x, y: math.log(max(y,1))),
        ("sign(x)*|x|^0.7 vs y^0.5", lambda x, y: math.copysign(abs(x)**0.7, x), lambda x, y: math.sqrt(y)),
        ("x*exp(y/Nmax) vs log(y)", lambda x, y: x*math.exp(y/N_max) if y > 0 else 0, lambda x, y: math.log(max(y,1))),
    ]

    results = []
    for label, xf, yf in parametrizations:
        xs = np.array([xf(px, py) for px, py in data])
        ys = np.array([yf(px, py) for px, py in data])

        # Filter valid
        valid = np.isfinite(xs) & np.isfinite(ys) & (np.abs(xs) < 1e6) & (ys < 1e6)
        xs_v = xs[valid]
        ys_v = ys[valid]

        k, A, phase, rms = fit_manifold(xs_v, ys_v, label)

        # Invariant: k * log(N_max) (for log(y) parametrizations)
        k_log_invariant = k * log_N_max if k is not None else None

        results.append({
            'parametrization': label,
            'k': k,
            'A': A,
            'phase': phase,
            'RMS': rms,
            'k_log_invariant': k_log_invariant,
            'n_points': len(xs_v)
        })

        phi = (1 + math.sqrt(5)) / 2
        diff_pct = abs(k - 1/phi) / (1/phi) * 100 if k is not None else None

        print(f"[{label}]")
        print(f"  k={k:.6f}  A={A:.4f}  phase={phase:.4f}  RMS={rms:.6f}")
        if diff_pct is not None:
            print(f"  diff from 1/φ: {diff_pct:.2f}%")
        if k_log_invariant is not None:
            print(f"  k*log(N_max) = {k_log_invariant:.6f}")
        print()

    # Write results CSV
    with open(OUTPUT_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['parametrization','k','A','phase','RMS','k_log_invariant','n_points'])
        writer.writeheader()
        writer.writerows(results)

    print(f"Results written to {OUTPUT_PATH}")

    # Summary statistics
    ks = [r['k'] for r in results if r['k'] is not None]
    if ks:
        print(f"\n=== Summary ===")
        print(f"k range: [{min(ks):.4f}, {max(ks):.4f}]  spread={max(ks)-min(ks):.4f}")
        print(f"k mean: {np.mean(ks):.4f}  std: {np.std(ks):.4f}")
        phi_inv = 1 / ((1+math.sqrt(5))/2)
        diffs = [abs(k - phi_inv)/phi_inv*100 for k in ks]
        print(f"% diff from 1/φ: min={min(diffs):.2f}%  max={max(diffs):.2f}%  mean={np.mean(diffs):.2f}%")

        if np.std(ks) < 0.05:
            print("\n✓ k is ROBUST across coordinate parametrizations (std < 0.05)")
        else:
            print("\n✗ k varies significantly across parametrizations (std >= 0.05) — may be artifact")


if __name__ == "__main__":
    run_sensitivity()
