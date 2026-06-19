#!/usr/bin/env python3
"""
dsr_multi_start_optimization.py
==============================
Multi-start optimization to map ALL local minima in the manifold fitting.

The key question: Is k≈0.6146 (oscillating, local MSE min) a robust finding
that appears from multiple starting points, or does it only appear from
k₀≈0.63 initialization?

Run from many k₀ points across [0.01, 2.0] to map the landscape.
"""

import numpy as np
import pandas as pd
import json
from scipy.optimize import minimize

DATA = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/traj_vmPhase9.csv"
OUT_PATH = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase10_DualStrandToZeta/python_code/multistart_results.json"


def load_data():
    df = pd.read_csv(DATA)
    counts = df.groupby('traj_id')['frame'].transform('count')
    dfc = df[counts >= 3].copy().dropna(subset=['x', 'y'])
    y = dfc['y'].values.astype(float)
    xmag = np.abs(dfc['x'].values.astype(float) - 0.5)
    return y, xmag


def loss_abs(params, y, xmag):
    k, A, phase = params
    pred = np.abs(A * np.cos(k * y + phase))
    return np.mean((xmag - pred) ** 2)


def fit_point(y, xmag, k0, A0_std):
    """Fit from one starting point. Return None if failed."""
    try:
        res = minimize(
            lambda p: loss_abs(p, y, xmag),
            x0=[k0, A0_std * 1.2, 0.0],
            bounds=[(0.01, 2.0), (0.01, 0.6), (-np.pi, np.pi)],
            method='L-BFGS-B'
        )
        return res.x, res.fun
    except:
        return None, None


def run_multistart():
    print("=== Multi-Start Optimization ===\n")
    y, xmag = load_data()
    A0_std = float(np.std(xmag))
    print(f"N={len(y)}, A0_std={A0_std:.4f}")
    print(f"Searching k0 in [0.01, 2.0] step 0.05...\n")

    # Grid of starting points
    k0_grid = np.arange(0.01, 2.05, 0.05)
    results = []

    for k0 in k0_grid:
        params, loss = fit_point(y, xmag, k0, A0_std)
        if params is not None:
            k, A, phase = params
            rms = np.sqrt(loss)
            results.append({
                'k0': k0,
                'k': k,
                'A': A,
                'phase': phase,
                'rms': rms,
            })

    print(f"Completed {len(results)}/{len(k0_grid)} optimizations\n")

    # Cluster results by final k value
    ks_final = [r['k'] for r in results]
    print(f"Final k range: [{min(ks_final):.4f}, {max(ks_final):.4f}]")
    print(f"Unique k clusters (within 0.1):")

    # Group by clustering
    clusters = []
    for r in results:
        placed = False
        for c in clusters:
            if abs(r['k'] - c['center']) < 0.1:
                c['members'].append(r)
                c['count'] += 1
                c['rms_sum'] += r['rms']
                placed = True
                break
        if not placed:
            clusters.append({'center': r['k'], 'members': [r], 'count': 1, 'rms_sum': r['rms']})

    clusters.sort(key=lambda c: c['count'], reverse=True)
    for c in clusters:
        avg_rms = c['rms_sum'] / c['count']
        print(f"  k≈{c['center']:.4f}: {c['count']} solutions, avg RMS={avg_rms:.6f}")
        print(f"    k0 range: {[m['k0'] for m in c['members'][:5]]}...")

    # Key question: does k≈0.6146 appear from diverse starting points?
    near_0614 = [r for r in results if 0.55 < r['k'] < 0.68]
    near_0025 = [r for r in results if r['k'] < 0.05]
    print(f"\n=== Key Answer ===")
    print(f"k≈0.6146 found from {len(near_0614)} starting points")
    print(f"k≈0.025 found from {len(near_0025)} starting points")

    if len(near_0614) > 10:
        print("→ k≈0.6146 is a ROBUST local minimum (found from many starts)")
    else:
        print("→ k≈0.6146 may be initialization-sensitive (only from k₀≈0.6-0.7)")

    # Save results
    out = {
        'n_starts': len(k0_grid),
        'n_success': len(results),
        'all_results': results,
        'clusters': [{'center': c['center'], 'count': c['count'], 'avg_rms': c['rms_sum']/c['count']} for c in clusters],
        'k0614_count': len(near_0614),
        'k0025_count': len(near_0025),
    }
    with open(OUT_PATH, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nResults written to {OUT_PATH}")


if __name__ == "__main__":
    run_multistart()
