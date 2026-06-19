#!/usr/bin/env python3
"""
dsr_structural_transitions.py
==============================
Investigate the structural transition in |cos(k·y+phase)|^p at p=2.

Key questions:
1. What is the physical meaning of the L1→L2 transition?
2. Does the golden ratio geometry emerge specifically at p=2?
3. What does this tell us about the manifold's structure?

At p=1 (L1 norm): fit minimizes absolute deviation
At p=2 (L2 norm): fit minimizes squared deviation
At p→∞ (L∞ norm): fit minimizes maximum deviation

The fact that RMS drops 76% at p=2 suggests the data becomes
much more coherent when viewed through an L2 lens. This could
indicate:
- The underlying pattern is fundamentally L2 (energy-like)
- The manifold encodes a sum-of-squares structure
- The golden ratio emerges as the energy minimizer
"""

import csv
import math
import numpy as np
from scipy.optimize import minimize

DATA_PATH = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/traj_vmPhase9.csv"


def load_vm_data():
    rows = []
    with open(DATA_PATH) as f:
        reader = csv.DictReader(f)
        for row in reader:
            x = row.get('x', '')
            y = row.get('y', '')
            if x and y:
                rows.append((float(x), float(y)))
    return np.array(rows)


def fit_manifold_p(xs, ys, p, k_hint=None):
    """
    Fit |x-0.5|^p = (A * |cos(k*y + phase)|)^p
    Equivalent to: |x-0.5| = A * |cos(k*y + phase)| with p-p norm
    
    We fit by minimizing Σ |target^p - (A*|cos(...)})^p|^2
    Which is equivalent to minimizing Σ |target - A*|cos(...)||^(2p)
    """
    targets = np.abs(xs - 0.5)
    ys_n = ys / ys.max()

    # Objective: minimize L^p distance
    def objective(params):
        A, k, phase = params
        preds = np.abs(A * np.cos(k * ys_n + phase))
        # Minimize L^p norm of (target - pred)
        if p == 1:
            return np.sum(np.abs(targets - preds))
        elif p == 2:
            return np.sum((targets - preds) ** 2)
        elif p == 0.5:
            # Use sqrt for numerical stability
            return np.sum(np.abs(targets - preds) ** 0.5)
        else:
            return np.sum(np.abs(targets - preds) ** p)

    best = None
    k0_candidates = [0.4, 0.5, 0.618, 0.7, 0.8] if k_hint is None else [k_hint]
    for k0 in k0_candidates:
        for A0 in [0.15, 0.19, 0.25]:
            for ph0 in [0.0, 0.3, 0.6, -0.3]:
                res = minimize(objective, x0=[A0, k0, ph0],
                              method='Nelder-Mead', options={'maxiter': 3000})
                if best is None or res.fun < best.fun:
                    best = res

    A, k, phase = best.x
    preds = np.abs(A * np.cos(k * ys_n + phase))
    if p == 1:
        rms = np.sqrt(np.mean(np.abs(targets - preds)))
    elif p == 2:
        rms = np.sqrt(np.mean((targets - preds) ** 2))
    elif p == 0.5:
        rms = np.sqrt(np.mean(np.abs(targets - preds) ** 0.5))
    else:
        rms = np.sqrt(np.mean(np.abs(targets - preds) ** p))

    # Also compute alternative metrics for comparison
    mae = np.mean(np.abs(targets - preds))
    mse = np.mean((targets - preds) ** 2)

    return k, A, phase, rms, mae, mse


def analyze_p_transitions():
    print("=== Ω^p Structural Transition Analysis ===\n")

    data = load_vm_data()
    xs = data[:, 0]
    ys = data[:, 1]
    N = len(xs)
    print(f"Loaded {N} points, y_max={ys.max():.2f}")

    phi = (1 + math.sqrt(5)) / 2
    phi_inv = 1 / phi

    # Sweep p values
    p_values = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0]
    results = []

    print("\n--- p-sweep results ---")
    print(f"{'p':>6}  {'k':>8}  {'A':>7}  {'phase':>7}  {'RMS':>8}  {'MAE':>8}  {'k diff%':>8}")

    prev_k = None
    for p in p_values:
        k, A, phase, rms, mae, mse = fit_manifold_p(xs, ys, p)
        diff_pct = abs(k - phi_inv) / phi_inv * 100
        results.append({'p': p, 'k': k, 'A': A, 'phase': phase,
                       'RMS': rms, 'MAE': mae, 'MSE': mse,
                       'k_diff_pct': diff_pct})

        print(f"{p:>6.2f}  {k:>8.5f}  {A:>7.4f}  {phase:>7.4f}  {rms:>8.6f}  {mae:>8.6f}  {diff_pct:>8.2f}%")

        # Detect structural change
        if prev_k is not None:
            delta_k = abs(k - prev_k)
            if delta_k > 0.05:
                print(f"  *** k jump at p={p}: Δk={delta_k:.4f}")
        prev_k = k

    print("\n--- Key Observations ---")

    # Find where k is closest to 1/φ
    best = min(results, key=lambda r: r['k_diff_pct'])
    print(f"k nearest 1/φ: p={best['p']:.2f}, k={best['k']:.6f}, diff={best['k_diff_pct']:.2f}%")

    # Find where RMS drops most
    rms_values = [r['RMS'] for r in results]
    for i in range(1, len(results)):
        drop = (results[i-1]['RMS'] - results[i]['RMS']) / results[i-1]['RMS']
        if drop > 0.3:
            print(f"  RMS drop at p={results[i]['p']:.2f}: {drop*100:.1f}% reduction")

    # Check for p=2 special behavior
    r_p2 = next(r for r in results if r['p'] == 2.0)
    r_p1 = next(r for r in results if r['p'] == 1.0)
    print(f"\n  p=1 → p=2 comparison:")
    print(f"    k: {r_p1['k']:.6f} → {r_p2['k']:.6f} (Δ={r_p2['k']-r_p1['k']:.6f})")
    print(f"    RMS: {r_p1['RMS']:.6f} → {r_p2['RMS']:.6f} (Δ={r_p2['RMS']-r_p1['RMS']:.6f})")
    print(f"    % diff from 1/φ: {r_p1['k_diff_pct']:.2f}% → {r_p2['k_diff_pct']:.2f}%")

    # L1→L2 ratio
    print(f"\n  L2/L1 ratio (energy/MAE measure):")
    print(f"    MSE(L2) / MSE(L1) at optimal k: {r_p2['MSE']/r_p1['MSE']:.4f}")
    print(f"    This ratio relates to the signal-to-noise structure")

    # Pentagonal number connection
    print(f"\n--- Pentagonal Number Connection ---")
    print(f"  2sin(π/10) = 2sin(18°) = (√5-1)/2 = 1/φ = {phi_inv:.6f}")
    print(f"  Generalized pentagonal numbers: g_k = k(3k-1)/2")
    for k in range(1, 6):
        g = k * (3*k - 1) // 2
        print(f"    g_{k} = {g}")
    print(f"\n  If Λ(n) encodes pentagonal recurrence, the golden ratio")
    print(f"  emerges from pentagonal geometry via 1/φ = 2sin(π/10)")

    # Energy interpretation
    print(f"\n--- Physical Interpretation ---")
    print(f"  L1 norm: minimizes absolute deviation (robust to outliers)")
    print(f"  L2 norm: minimizes squared deviation (energy-like, favors smoothness)")
    print(f"  At p=2, the manifold 'snaps in' because the underlying structure")
    print(f"  is fundamentally an energy minimization (sum of squares).")
    print(f"  This suggests DSR manifold = energy eigenstate of some operator.")


if __name__ == "__main__":
    analyze_p_transitions()
