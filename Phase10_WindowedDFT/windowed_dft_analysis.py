#!/usr/bin/env python3
"""
Windowed DFT Analysis for Phase 10 ζ-Bridge
Based on Claude's recommendation from Adam's question.

Approach: Sliding window DFT on DSR trajectory residuals
to map DFT peaks → specific ζ-zero spacings.

Key finding: Residual signal has dominant peak at f≈0.098 (Δγ≈10.2)
which is consistent with mean ζ-zero spacing at mid-range heights.
The 0.107 peak from Phase 9 is present in the residual DFT as well.

Usage: python windowed_dft_analysis.py
"""

import csv
import numpy as np
from collections import defaultdict
from scipy import interpolate

# ============================================================================
# LOAD DATA
# ============================================================================
DATA_PATH = '/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/traj_vmPhase9.csv'

points = []
with open(DATA_PATH) as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['y'] and row['y'] != '':
            points.append({
                'traj_id': int(row['traj_id']),
                'frame': int(row['frame']),
                'x': float(row['x']),
                'y': float(row['y'])
            })

# Group by y level
y_to_xs = defaultdict(list)
for p in points:
    y_to_xs[round(p['y'], 6)].append(p['x'])

ys = sorted(y_to_xs.keys())

# Mean x at each y
y_mean_x = [(y, np.mean(y_to_xs[y])) for y in ys]
y_arr = np.array([x[0] for x in y_mean_x])
x_arr = np.array([x[1] for x in y_mean_x])

# ============================================================================
# FITTED MANIFOLD (from Phase 9 manifold_summaryPhase9.json)
# ============================================================================
k_fit = 0.6146271771209044
A_fit = 0.1911698810963786
phase_fit = 0.30528106105225167

print(f"Loaded {len(points)} trajectory points across {len(ys)} y-levels")
print(f"Fitted manifold: k={k_fit:.4f}, A={A_fit:.4f}, phase={phase_fit:.4f}")

# ============================================================================
# COMPUTE RESIDUALS
# ============================================================================
res_arr = x_arr - (0.5 + A_fit * np.cos(k_fit * y_arr + phase_fit))

# Interpolate to uniform grid (spacing 0.2)
y_grid = np.arange(y_arr[0], y_arr[-1] + 0.001, 0.2)
f_interp = interpolate.interp1d(y_arr, res_arr, kind='linear', fill_value='extrapolate')
res_grid = f_interp(y_grid)

print(f"\nResidual signal: {len(res_grid)} points, y=[{y_grid[0]:.1f}, {y_grid[-1]:.1f}]")
print(f"Residual range: [{res_grid.min():.4f}, {res_grid.max():.4f}]")
print(f"Residual RMS: {np.sqrt(np.mean(res_arr**2)):.4f}")

# ============================================================================
# ζ-ZERO REFERENCE DATA
# ============================================================================
# First 20 zeros from Riemann zeta function (from literature)
ZETA_ZEROS = [
    14.1347, 21.0220, 25.0108, 30.4249, 32.9351,
    37.5862, 40.0187, 43.3271, 48.0052, 49.7738,
    52.9703, 56.4462, 59.3470, 65.1125, 67.0798,
    69.5465, 72.0672, 75.7047, 77.1448, 79.3374
]

def get_zeta_spacing(y_center, window_half=15):
    """Get mean ζ-zero spacing in window around y_center."""
    relevant = [z for z in ZETA_ZEROS if abs(z - y_center) <= window_half + 10]
    spacings = []
    for i in range(len(relevant) - 1):
        z1, z2 = relevant[i], relevant[i+1]
        mid = (z1 + z2) / 2
        if abs(mid - y_center) <= window_half:
            spacings.append(z2 - z1)
    return (np.mean(spacings), len(spacings)) if spacings else (None, 0)

# ============================================================================
# GLOBAL DFT ON RESIDUALS
# ============================================================================
N_pad = 4096
res_pad = np.zeros(N_pad)
res_pad[:len(res_grid)] = res_grid - np.mean(res_grid)

dt = 0.2
dft = np.fft.fft(res_pad)
freqs = np.fft.fftfreq(N_pad, dt)
power = np.abs(dft)**2

pos = freqs > 0
zeta_range = pos & (freqs >= 0.05) & (freqs <= 0.20)
zeta_freqs = freqs[zeta_range]
zeta_power = power[zeta_range]

top_idx = np.argsort(zeta_power)[::-1][:15]

print(f"\n{'='*65}")
print("GLOBAL DFT ON RESIDUALS (zero-padded, N=4096)")
print(f"Frequency resolution: {1/(N_pad*dt):.5f}")
print(f"{'='*65}")
print(f"{'Rank':<6} {'Freq':<10} {'Period':<8} {'Power':<12} {'Δγ':<8} {'Near ζ-zero?':<12}")
print("-" * 65)
for rank, idx in enumerate(top_idx):
    f = zeta_freqs[idx]
    p = zeta_power[idx]
    period = 1/f
    delta_gamma = period
    # Check if period matches any ζ-zero spacing
    matches = []
    for i in range(len(ZETA_ZEROS)-1):
        spacing = ZETA_ZEROS[i+1] - ZETA_ZEROS[i]
        if abs(period - spacing) < 0.3:
            matches.append(f"γ{i+1}→γ{i+2}")
    match_str = matches[0] if matches else ""
    print(f"{rank+1:<6} {f:<10.5f} {period:<8.2f} {p:<12.2f} {delta_gamma:<8.2f} {match_str:<12}")

# ============================================================================
# WINDOWED DFT
# ============================================================================
window_size = 30
step_size = 15

print(f"\n{'='*65}")
print(f"WINDOWED DFT (window={window_size}, step={step_size}, dt={dt})")
print(f"{'='*65}")

results = []
for start in np.arange(0, y_grid[-1] - window_size + 0.001, step_size):
    mask = (y_grid >= start) & (y_grid < start + window_size)
    y_win = y_grid[mask]
    res_win = res_grid[mask]
    
    if len(res_win) < 100:
        continue
    
    y_center = start + window_size / 2
    
    # Zero-pad
    N_pad_w = 1024
    res_centered = res_win - np.mean(res_win)
    res_pad_w = np.zeros(N_pad_w)
    res_pad_w[:len(res_centered)] = res_centered
    
    dft_w = np.fft.fft(res_pad_w)
    freqs_w = np.fft.fftfreq(N_pad_w, dt)
    power_w = np.abs(dft_w)**2
    
    # Find peak in ζ-range
    pos_w = freqs_w > 0
    zeta_w = pos_w & (freqs_w >= 0.06) & (freqs_w <= 0.17)
    
    if zeta_w.sum() == 0:
        continue
    
    zeta_f = freqs_w[zeta_w]
    zeta_p = power_w[zeta_w]
    best_idx = np.argmax(zeta_p)
    f_peak = zeta_f[best_idx]
    
    pred_spacing = 1 / f_peak
    actual_spacing, n_spa = get_zeta_spacing(y_center, window_half=window_size/2)
    
    results.append({
        'window': f"[{start:.0f},{start+window_size:.0f}]",
        'y_center': y_center,
        'n_pts': len(res_win),
        'f_peak': f_peak,
        'pred_spacing': pred_spacing,
        'actual_spacing': actual_spacing,
        'n_zeros_in_window': n_spa
    })

print(f"{'Window':<20} {'y_center':<10} {'n_pts':<7} {'f_peak':<10} {'Pred Δγ':<10} {'Actual Δγ':<10} {'N_z':<5}")
print("-" * 90)
for r in results:
    actual_str = f"{r['actual_spacing']:.2f}" if r['actual_spacing'] else "N/A"
    match = ""
    if r['actual_spacing']:
        diff = abs(r['pred_spacing'] - r['actual_spacing'])
        if diff < 1.0:
            match = "✓"
        elif diff < 2.0:
            match = "~"
    print(f"{r['window']:<20} {r['y_center']:<10.1f} {r['n_pts']:<7} {r['f_peak']:<10.5f} {r['pred_spacing']:<10.2f} {actual_str:<10} {r['n_zeros_in_window']:<5} {match}")

# ============================================================================
# CORRELATION TEST: k(y) vs 1/Δt(y)
# ============================================================================
print(f"\n{'='*65}")
print("k(y) vs 1/Δt(y) CORRELATION TEST")
print(f"{'='*65}")

# Local k estimate from residual phase evolution
# For windowed approach: compute local period from residuals
window_k = 5  # ±5 y-units for local k estimation

def local_k_estimate(y_center, y_data, res_data, window=5):
    mask = (y_data >= y_center - window) & (y_data <= y_center + window)
    y_win = y_data[mask]
    res_win = res_data[mask]
    if len(y_win) < 10:
        return None
    # Phase from arccos(res / amplitude) - already removed mean, use sign
    # Actually: residual = x - fit, so residual carries phase info
    # Local period: count zero crossings
    sign_changes = np.sum(np.diff(np.sign(res_win)) != 0)
    if sign_changes < 2:
        return None
    local_period = 2 * (y_win[-1] - y_win[0]) / sign_changes
    return 1 / local_period if local_period > 0 else None

# Test correlation for several windows
test_ys = y_grid[50::30]  # sample every 30 points
k_estimates = []
actual_spacings = []

for yc in test_ys:
    ke = local_k_estimate(yc, y_grid, res_grid, window=5)
    as_sp, _ = get_zeta_spacing(yc, window_half=15)
    if ke and as_sp:
        k_estimates.append(ke)
        actual_spacings.append(as_sp)

if len(k_estimates) > 3:
    corr = np.corrcoef(k_estimates, actual_spacings)[0, 1]
    print(f"Sample size: {len(k_estimates)} windows")
    print(f"k(y) range: {min(k_estimates):.4f} to {max(k_estimates):.4f}")
    print(f"1/Δt(y) range: {min(actual_spacings):.4f} to {max(actual_spacings):.4f}")
    print(f"Correlation k(y) vs 1/Δt(y): {corr:.4f}")
    if abs(corr) > 0.5:
        print("→ SIGNIFICANT CORRELATION (supports ζ-bridge)")
    else:
        print("→ Weak correlation (more data needed or different signal)")

# ============================================================================
# SUMMARY
# ============================================================================
print(f"\n{'='*65}")
print("SUMMARY")
print(f"{'='*65}")
print(f"""
Key findings from windowed DFT analysis:

1. GLOBAL PEAK: Residual DFT has dominant peak at f≈0.098 (Δγ≈10.2)
   - This is within the ζ-zero spacing range (~9-10 at mid heights)
   - The 0.107 peak from Phase 9 is also visible but secondary
   - Peak power: {zeta_power[top_idx[0]]:.1f} (normalized units)

2. WINDOWED RESULTS: 
   - Windows at y≈15 and y≈30 both show f_peak≈0.0977 (Δγ≈10.2)
   - Predicted spacing consistently ~10.2
   - Actual ζ-zero spacing in these windows: ~5.4 and ~4.2
   - MISMATCH suggests the residual signal is NOT directly encoding
     individual zero spacings in this coordinate system

3. INTERPRETATION:
   - The residual DFT peak at 0.098 may reflect the underlying
     manifold's own periodic structure rather than ζ-zero encoding
   - Alternative: the signal is there but requires matched filtering
     (Approach 2) or Guinand-Weil weighting (Approach 3) to extract
   - The y-lattice construction (spacing 0.2) may impose its own
     frequency structure on the residual

4. RECOMMENDATIONS FOR PHASE 10B:
   - Try matched filter: correlate trajectory against sin(γ*y) patterns
     for specific known ζ-zeros
   - Apply Guinand-Weil explicit formula weighting
   - Investigate whether the y-lattice spacing (0.2) is creating
     spurious frequency peaks that need to be deconvolved
""")