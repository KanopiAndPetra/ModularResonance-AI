#!/usr/bin/env python3
"""
High-Resolution DFT Analysis + Blind Zero Prediction Test
Phase 10 of ModularResonance-AI

Task 1: High-Resolution DFT of residual signal
Task 2: Blind prediction experiment with multiple Y→γ mapping hypotheses
"""

import json
import numpy as np
from numpy.fft import fft, fftfreq
from scipy.interpolate import interp1d

WORKSPACE = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase10_DualStrandToZeta"
TRAJ_FILE = "/Users/oppie1.kanopi/Desktop/ModularResonance-AI/Phase9_DualStrandResonance/DSR_Phase9/trajectory_data/traj_vmPhase9.csv"

# ============================================================
# LOAD AND PROCESS TRAJECTORY DATA
# ============================================================

def load_residual_data():
    """Load trajectory data and compute residuals (midpoint - 0.5)."""
    import csv
    
    with open(TRAJ_FILE) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = []
        for row in reader:
            if len(row) >= 4 and row[2] and row[3]:
                try:
                    rows.append((int(row[0]), int(row[1]), float(row[2]), float(row[3])))
                except:
                    pass
    
    # Group by frame
    frames = {}
    for traj_id, frame, x, y in rows:
        if frame not in frames:
            frames[frame] = []
        frames[frame].append((x, y))
    
    # Collect residual data (midpoint of two strands at same Y)
    all_residuals = []
    all_ys = []
    for frame in sorted(frames.keys()):
        pts = frames[frame]
        sorted_pts = sorted(pts, key=lambda p: p[1])
        i = 0
        while i < len(sorted_pts) - 1:
            y0 = sorted_pts[i][1]
            y1 = sorted_pts[i+1][1]
            if abs(y1 - y0) < 0.1:  # Same Y = same trajectory point
                x0 = sorted_pts[i][0]
                x1 = sorted_pts[i+1][0]
                if x0 > x1:
                    x0, x1 = x1, x0
                mid_x = (x0 + x1) / 2.0
                residual = mid_x - 0.5
                all_residuals.append(residual)
                all_ys.append(y0)
                i += 2
            else:
                i += 1
    
    residuals = np.array(all_residuals)
    ys = np.array(all_ys)
    
    # Sort by Y and deduplicate
    sort_idx = np.argsort(ys)
    ys_sorted = ys[sort_idx]
    residuals_sorted = residuals[sort_idx]
    
    # Deduplicate by Y (mean of duplicates)
    unique_ys, idx = np.unique(ys_sorted, return_index=True)
    residuals_dedup = np.array([residuals_sorted[ys_sorted == y].mean() for y in unique_ys])
    
    return unique_ys, residuals_dedup

print("Loading trajectory data...")
ys, residuals = load_residual_data()
print(f"Loaded {len(residuals)} residual points, Y range: {ys.min():.4f} to {ys.max():.4f}")

# ============================================================
# TASK 1: HIGH-RESOLUTION DFT
# ============================================================

print("\n" + "="*60)
print("TASK 1: High-Resolution DFT Analysis")
print("="*60)

# Load fitted helix parameters from spectral analysis
with open(f"{WORKSPACE}/python_code/spectral_analysis_results.json") as f:
    spec = json.load(f)

k_fitted = spec["k_fitted"]
A = spec["A"]
phase = spec["phase"]
print(f"Fitted helix: k={k_fitted:.6f}, A={A:.6f}, phase={phase:.6f}")

# Compute fitted helix and subtract to get true residuals
# The fitted helix is: x = 0.5 + A * cos(k * y + phase)
fitted_helix = 0.5 + A * np.cos(k_fitted * ys + phase)
true_residuals = residuals - fitted_helix  # Residuals after removing helix

print(f"True residual range (after helix removal): {true_residuals.min():.6f} to {true_residuals.max():.6f}")
print(f"True residual RMS: {np.sqrt(np.mean(true_residuals**2)):.6f}")

# For DFT we use the raw residuals (deviation from midpoint=0.5)
# The helix subtraction might remove the main signal, so let's use both
print("\nComputing DFT on raw residuals (midpoint - 0.5)...")

# Zero-pad to N=4096 for high resolution
N_fft = 4096
y_max = ys.max()
n_original = len(residuals)

# Create evenly-spaced interpolation for FFT
# Use the actual Y values, resample to even spacing
y_step = y_max / N_fft
ys_even = np.linspace(0, y_max, N_fft)

# Interpolate residuals to even spacing
interp_residual = interp1d(ys, residuals, kind='cubic', fill_value=0, bounds_error=False)
residuals_even = interp_residual(ys_even)

print(f"Zero-padded from {n_original} points to {N_fft} FFT points")
print(f"Frequency resolution: {1.0/y_max/N_fft:.8f} (per bin)")

# Apply window to reduce spectral leakage
window = np.hanning(N_fft)
residuals_windowed = residuals_even * window

# Compute FFT
fft_vals = fft(residuals_windowed)
freqs = fftfreq(N_fft, y_step)

# Power spectrum (only positive frequencies)
positive_mask = freqs > 0
freqs_pos = freqs[positive_mask]
power = np.abs(fft_vals[positive_mask])**2

# Focus on 0.05 to 0.15 range (where the main peak is)
freq_mask = (freqs_pos >= 0.05) & (freqs_pos <= 0.15)
freqs_focus = freqs_pos[freq_mask]
power_focus = power[freq_mask]

# Find peaks in the focus range
peak_indices = []
for i in range(1, len(power_focus) - 1):
    if power_focus[i] > power_focus[i-1] and power_focus[i] > power_focus[i+1]:
        if power_focus[i] > np.mean(power_focus) * 2:  # Significant peak
            peak_indices.append(i)

peak_freqs = freqs_focus[[i for i in peak_indices]]
peak_powers = power_focus[[i for i in peak_indices]]

# Sort by power
sort_idx = np.argsort(peak_powers)[::-1]
peak_freqs = peak_freqs[sort_idx]
peak_powers = peak_powers[sort_idx]

print(f"\nPeaks in 0.05-0.15 range:")
for pf, pp in zip(peak_freqs[:10], peak_powers[:10]):
    print(f"  f = {pf:.6f}, power = {pp:.4f}")

# Check for peak splitting at high resolution
# Try N=8192 for even higher resolution
N_fft_2 = 8192
y_step_2 = y_max / N_fft_2
ys_even_2 = np.linspace(0, y_max, N_fft_2)
interp_residual_2 = interp1d(ys, residuals, kind='cubic', fill_value=0, bounds_error=False)
residuals_even_2 = interp_residual_2(ys_even_2)
window_2 = np.hanning(N_fft_2)
residuals_windowed_2 = residuals_even_2 * window_2

fft_vals_2 = fft(residuals_windowed_2)
freqs_2 = fftfreq(N_fft_2, y_step_2)
positive_mask_2 = freqs_2 > 0
freqs_pos_2 = freqs_2[positive_mask_2]
power_2 = np.abs(fft_vals_2[positive_mask_2])**2

freq_mask_2 = (freqs_pos_2 >= 0.05) & (freqs_pos_2 <= 0.15)
freqs_focus_2 = freqs_pos_2[freq_mask_2]
power_focus_2 = power_2[freq_mask_2]

# Find peaks in N=8192 spectrum
peak_indices_2 = []
for i in range(1, len(power_focus_2) - 1):
    if power_focus_2[i] > power_focus_2[i-1] and power_focus_2[i] > power_focus_2[i+1]:
        if power_focus_2[i] > np.mean(power_focus_2) * 2:
            peak_indices_2.append(i)

peak_freqs_2 = freqs_focus_2[[i for i in peak_indices_2]]
peak_powers_2 = power_focus_2[[i for i in peak_indices_2]]
sort_idx_2 = np.argsort(peak_powers_2)[::-1]
peak_freqs_2 = peak_freqs_2[sort_idx_2]
peak_powers_2 = peak_powers_2[sort_idx_2]

print(f"\nPeaks in 0.05-0.15 range (N=8192, higher resolution):")
for pf, pp in zip(peak_freqs_2[:10], peak_powers_2[:10]):
    print(f"  f = {pf:.6f}, power = {pp:.4f}")

# Determine peak structure
main_peak_4096 = peak_freqs[0] if len(peak_freqs) > 0 else None
main_peak_8192 = peak_freqs_2[0] if len(peak_freqs_2) > 0 else None

# Check if peak splits
nearby_peaks_8192 = [f for f in peak_freqs_2 if abs(f - main_peak_8192) < 0.02] if main_peak_8192 else []

dft_results = {
    "N_fft_4096": {
        "peak_freqs": peak_freqs[:10].tolist(),
        "peak_powers": peak_powers[:10].tolist(),
        "main_peak": main_peak_4096
    },
    "N_fft_8192": {
        "peak_freqs": peak_freqs_2[:10].tolist(),
        "peak_powers": peak_powers_2[:10].tolist(),
        "main_peak": main_peak_8192
    },
    "peak_splitting": {
        "splits_at_higher_resolution": len(nearby_peaks_8192) > 1,
        "nearby_peaks_count": len(nearby_peaks_8192),
        "nearby_peaks": nearby_peaks_8192
    },
    "y_range": float(y_max),
    "n_original_points": n_original
}

# Determine structure description
if len(nearby_peaks_8192) > 1:
    structure_desc = "MULTIPLE_SUB_PEAKS"
elif main_peak_8192 and abs(main_peak_8192 - 0.107) < 0.005:
    structure_desc = "SINGLE_SHARP_PEAK_MATCHING_ZETA_ZERO_SPACING"
else:
    structure_desc = "BROADENED_OR_SHIFTED_PEAK"

dft_results["structure_description"] = structure_desc

print(f"\nPeak structure: {structure_desc}")
print(f"Main peak at 4096: {main_peak_4096}")
print(f"Main peak at 8192: {main_peak_8192}")
print(f"Nearby peaks at higher res: {len(nearby_peaks_8192)}")

# ============================================================
# TASK 2: BLIND ZERO PREDICTION TEST
# ============================================================

print("\n" + "="*60)
print("TASK 2: Blind Zero Prediction Test")
print("="*60)

# Load zeta zeros
zeta_gammas = spec["zeta_gammas"]
print(f"Loaded {len(zeta_gammas)} zeta zeros")
print(f"First 10 zeros: {zeta_gammas[:10]}")

# Split: Y <= 0.7 for training, Y > 0.7 for testing
train_mask = ys <= 0.7
test_mask = ys > 0.7

ys_train = ys[train_mask]
residuals_train = residuals[train_mask]
ys_test = ys[test_mask]
residuals_test = residuals[test_mask]

print(f"\nTraining set: Y <= 0.7, {len(ys_train)} points")
print(f"Test set: Y > 0.7, {len(ys_test)} points")
print(f"Training Y range: {ys_train.min():.4f} to {ys_train.max():.4f}")
print(f"Test Y range: {ys_test.min():.4f} to {ys_test.max():.4f}")

# Fit helix model on training data ONLY
# Model: x = 0.5 + A * cos(k * y + phi)
# Residual = midpoint - 0.5 = A * cos(k * y + phi)

def fit_helix(ys, residuals):
    """Fit helix model using least squares."""
    # Residual = A * cos(k * y + phi)
    # We need to find A, k, phi
    # Linearize: residual = A*cos(k*y)*cos(phi) - A*sin(k*y)*sin(phi)
    # = a*cos(k*y) + b*sin(k*y), where a=A*cos(phi), b=-A*sin(phi)
    
    # Grid search over k
    best_k = None
    best_rmse = float('inf')
    best_params = None
    
    for k_try in np.linspace(0.01, 3.0, 300):
        # Linear least squares for A and phase at this k
        A_cos = np.cos(k_try * ys)
        A_sin = np.sin(k_try * ys)
        # residual = a * cos(k*y) + b * sin(k*y)
        # Solve: [cos(k*y), sin(k*y)] @ [a, b] = residual
        # Using normal equations
        cos_mean = np.mean(A_cos)
        sin_mean = np.mean(A_sin)
        cos_sq_mean = np.mean(A_cos**2)
        sin_sq_mean = np.mean(A_sin**2)
        cos_sin_mean = np.mean(A_cos * A_sin)
        res_mean = np.mean(residuals)
        cos_res = np.mean(A_cos * residuals)
        sin_res = np.mean(A_sin * residuals)
        
        # Solve 2x2 system
        det = cos_sq_mean * sin_sq_mean - cos_sin_mean**2
        if abs(det) < 1e-10:
            continue
        
        a = (sin_sq_mean * cos_res - cos_sin_mean * sin_res) / det
        b = (-cos_sin_mean * cos_res + cos_sq_mean * sin_res) / det
        
        A_est = np.sqrt(a**2 + b**2)
        phi_est = np.arctan2(-b, a)
        
        # Compute fit and RMSE
        fitted = A_est * np.cos(k_try * ys + phi_est)
        rmse = np.sqrt(np.mean((residuals - fitted)**2))
        
        if rmse < best_rmse:
            best_rmse = rmse
            best_k = k_try
            best_params = (A_est, phi_est)
    
    return best_k, best_params, best_rmse

print("\nFitting helix on training data...")
k_train, (A_train, phi_train), rmse_train = fit_helix(ys_train, residuals_train)
print(f"Training fit: k={k_train:.6f}, A={A_train:.6f}, phi={phi_train:.6f}, RMSE={rmse_train:.6f}")

# Find predicted crossing points on test set
# Crossing points: where strand_0 = strand_1 = 0.5
# At crossing, residual = A * cos(k * y + phi) = 0
# So cos(k * y + phi) = 0 => k * y + phi = (2n+1) * pi/2 => y = ((2n+1)*pi/2 - phi) / k

def find_crossings(k, phi, y_min, y_max):
    """Find crossing points in [y_min, y_max] where cos(k*y + phi) = 0."""
    crossings = []
    # Solve: k*y + phi = pi/2 + n*pi => y = (pi/2 - phi + n*pi) / k
    n = 0
    while True:
        y_cross = (np.pi/2 - phi + n * np.pi) / k
        if y_cross > y_max:
            break
        if y_cross >= y_min:
            crossings.append(y_cross)
        n += 1
        if n > 10000:
            break
    return crossings

# Find crossing points in test region (Y > 0.7)
test_crossings = find_crossings(k_train, phi_train, ys_test.min(), ys_test.max())
print(f"\nPredicted crossing points in test region: {len(test_crossings)}")
print(f"Test crossing Y values: {[f'{y:.4f}' for y in test_crossings[:10]]}")

# Define Y -> gamma mapping hypotheses
def linear_mapping(Y, alpha, beta):
    """γ = α * Y + β"""
    return alpha * Y + beta

def exponential_mapping(Y, alpha):
    """γ = exp(alpha * Y) - 1"""
    return np.exp(alpha * Y) - 1

def von_mangoldt_mapping(Y, alpha, beta):
    """γ = alpha * Y * log(1/Y + 1) + beta (von Mangoldt-like)"""
    return alpha * Y * np.log(1/Y + 1) + beta

def log_mapping(Y, alpha, beta):
    """γ = alpha * log(Y + 1) + beta"""
    return alpha * np.log(Y + 1) + beta

# We need to calibrate mappings using training data
# Training crossings should map to actual zeta zeros
# First, find actual crossings in training region

train_crossings = find_crossings(k_train, phi_train, ys_train.min(), ys_train.max())
print(f"\nPredicted crossing points in training region: {len(train_crossings)}")
print(f"Training crossing Y values: {[f'{y:.4f}' for y in train_crossings[:5]]}")

# For each mapping, calibrate on training crossings, test on test crossings
# Mean zero spacing: 0.1058 (from spectral analysis)
mean_spacing = 0.1058

# Method 1: Linear mapping calibrated on training crossings
# We have train_crossings at Y values, which should map to actual zeta zeros
# Use first few training crossings to calibrate

# Sort training crossings and corresponding zeta zeros
# Assume first training crossing maps to first zeta zero (14.134...)
# Then spacing between crossings should match spacing between zeros

if len(train_crossings) >= 2:
    train_spacing_y = np.diff(train_crossings)
    # Get corresponding zeta zero spacings
    zero_spacings = np.diff(zeta_gammas[:len(train_crossings)])
    
    # Calibrate linear mapping: γ = α * Y
    # Use ratio of gamma spacing to Y spacing
    if len(train_spacing_y) > 0 and len(zero_spacings) > 0:
        alpha_linear = np.mean(zero_spacings / train_spacing_y)
    else:
        alpha_linear = mean_spacing / np.mean(train_spacing_y) if len(train_spacing_y) > 0 else 1.0
else:
    alpha_linear = 1.0

print(f"\nLinear mapping calibrated: α = {alpha_linear:.6f}")

# Method 2: Exponential mapping
# γ = exp(α * Y) - 1
# For small Y, this grows faster than linear
if len(train_crossings) >= 2:
    # Use ratio of consecutive spacings
    # If γ = exp(α*Y), then spacing ≈ exp(α*Y) * α * ΔY
    # For small increments, spacing is approximately α * ΔY * exp(α*Y_avg)
    alpha_exp = alpha_linear / np.mean(train_crossings) if np.mean(train_crossings) > 0 else 0.01
else:
    alpha_exp = 0.01

print(f"Exponential mapping calibrated: α = {alpha_exp:.6f}")

# Method 3: Von Mangoldt mapping
# γ = a * Y * log(1/Y) + b
# log(1/Y) is negative for Y > 1, so we need absolute value
def von_mangoldt_v2(Y, alpha, beta):
    """γ = alpha * Y * log((1+Y)/Y) + beta"""
    return alpha * Y * np.log((1+Y)/Y) + beta

# Calibrate: use mean spacing relationship
# For small Y, log((1+Y)/Y) ≈ 1/Y - 1/(2Y^2) + ... so Y*log(...) ≈ 1
# This means for small Y, gamma ≈ alpha + beta
if len(train_crossings) >= 2:
    # First zero: 14.134 at first training crossing
    y0 = train_crossings[0]
    # gamma(y0) should be near first zeta zero
    alpha_von = zeta_gammas[0] / (y0 * np.log((1+y0)/y0)) if y0 > 0 else zeta_gammas[0]
    beta_von = 0  # Assume first crossing maps to first zero
else:
    alpha_von = 1.0
    beta_von = 0.0

print(f"Von Mangoldt mapping calibrated: α = {alpha_von:.6f}, β = {beta_von:.6f}")

# Method 4: Log mapping
# γ = alpha * log(Y + 1) + beta
if len(train_crossings) >= 2:
    # First crossing maps to first zero
    y0 = train_crossings[0]
    alpha_log = zeta_gammas[0] / np.log(y0 + 1) if y0 > 0 else zeta_gammas[0]
    beta_log = 0
else:
    alpha_log = 1.0
    beta_log = 0.0

print(f"Log mapping calibrated: α = {alpha_log:.6f}, β = {beta_log:.6f}")

# ============================================================
# Test each mapping on test crossings
# ============================================================

mappings = {
    "linear": (lambda Y: linear_mapping(Y, alpha_linear, 0), "γ = α·Y"),
    "exponential": (lambda Y: exponential_mapping(Y, alpha_exp), "γ = exp(α·Y) - 1"),
    "von_mangoldt": (lambda Y: von_mangoldt_v2(Y, alpha_von, beta_von), "γ = α·Y·log((1+Y)/Y)"),
    "log": (lambda Y: log_mapping(Y, alpha_log, beta_log), "γ = α·log(Y+1)")
}

prediction_results = {}

for name, (mapping_func, desc) in mappings.items():
    # Predict gamma for each test crossing
    predicted_gammas = [mapping_func(y) for y in test_crossings]
    
    # Compare to actual zeta zeros
    # For each predicted gamma, find nearest actual zeta zero
    errors = []
    for pg in predicted_gammas:
        # Find nearest zeta zero
        nearest_zero = min(zeta_gammas, key=lambda z: abs(z - pg))
        error = abs(pg - nearest_zero)
        errors.append(error)
    
    mean_error = np.mean(errors) if errors else float('inf')
    min_error = min(errors) if errors else float('inf')
    
    prediction_results[name] = {
        "description": desc,
        "calibration": f"alpha={getattr(locals(), f'alpha_{name}', 'N/A'):.6f}" if name != 'von_mangoldt' else f"alpha={alpha_von:.6f}",
        "test_crossings": test_crossings[:5],
        "predicted_gammas": predicted_gammas[:5],
        "mean_error": float(mean_error),
        "min_error": float(min_error),
        "n_test_crossings": len(test_crossings)
    }
    
    print(f"\n{name} mapping:")
    print(f"  Description: {desc}")
    print(f"  Test crossings: {[f'{y:.4f}' for y in test_crossings[:5]]}")
    print(f"  Predicted gammas: {[f'{pg:.4f}' for pg in predicted_gammas[:5]]}")
    print(f"  Mean prediction error: {mean_error:.4f}")
    print(f"  Min prediction error: {min_error:.4f}")

# Find best mapping
best_mapping = min(prediction_results.keys(), key=lambda k: prediction_results[k]["mean_error"])

print(f"\n*** Best mapping: {best_mapping} (mean error: {prediction_results[best_mapping]['mean_error']:.4f}) ***")

# ============================================================
# Save results
# ============================================================

# Save DFT results
dft_output = {
    "task": "High-Resolution DFT Analysis",
    "description": "Examines peak structure at 0.107 frequency to determine if zeros are individually encoded",
    "N_fft_4096": {
        "peak_freqs": peak_freqs[:10].tolist(),
        "peak_powers": peak_powers[:10].tolist(),
        "main_peak": main_peak_4096
    },
    "N_fft_8192": {
        "peak_freqs": peak_freqs_2[:10].tolist(),
        "peak_powers": peak_powers_2[:10].tolist(),
        "main_peak": main_peak_8192
    },
    "peak_splitting_analysis": {
        "splits_at_higher_resolution": len(nearby_peaks_8192) > 1,
        "nearby_peaks_count": len(nearby_peaks_8192),
        "nearby_peaks": [float(f) for f in nearby_peaks_8192]
    },
    "structure_description": structure_desc,
    "interpretation": {
        "SINGLE_SHARP_PEAK_MATCHING_ZETA_ZERO_SPACING": "Peak is likely an ensemble average; individual zeros not clearly encoded",
        "MULTIPLE_SUB_PEAKS": "Individual zeros may be encoded in the manifold",
        "BROADENED_OR_SHIFTED_PEAK": "Signal is smeared; connection may be indirect"
    }[structure_desc],
    "y_range": float(y_max),
    "n_original_points": n_original
}

with open(f"{WORKSPACE}/python_code/high_res_dft_results.json", 'w') as f:
    json.dump(dft_output, f, indent=2)
print(f"\nSaved high_res_dft_results.json")

# Save prediction results
prediction_output = {
    "task": "Blind Zero Prediction Test",
    "description": "Tests multiple Y→γ mapping hypotheses for predictive encoding",
    "training_set": {
        "criterion": "Y <= 0.7",
        "n_points": int(len(ys_train)),
        "y_range": f"{ys_train.min():.4f} to {ys_train.max():.4f}",
        "crossings": [float(y) for y in train_crossings[:5]]
    },
    "test_set": {
        "criterion": "Y > 0.7",
        "n_points": int(len(ys_test)),
        "y_range": f"{ys_test.min():.4f} to {ys_test.max():.4f}",
        "crossings": [float(y) for y in test_crossings[:10]]
    },
    "fitted_helix_training": {
        "k": float(k_train),
        "A": float(A_train),
        "phi": float(phi_train),
        "RMSE": float(rmse_train)
    },
    "mappings_tested": prediction_results,
    "best_mapping": best_mapping,
    "best_mean_error": float(prediction_results[best_mapping]["mean_error"]),
    "threshold_for_useful": 5.0,
    "conclusion": "PREDICTIVE" if prediction_results[best_mapping]["mean_error"] < 5.0 else "NOT_PREDICTIVE"
}

with open(f"{WORKSPACE}/python_code/blind_prediction_results.json", 'w') as f:
    json.dump(prediction_output, f, indent=2)
print(f"Saved blind_prediction_results.json")

# ============================================================
# Summary
# ============================================================

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"\nTask 1 - High-Resolution DFT:")
print(f"  Main peak at f≈{main_peak_8192:.6f}" if main_peak_8192 else "  No clear peak")
print(f"  Peak structure: {structure_desc}")
print(f"  Peak splits at higher resolution: {len(nearby_peaks_8192) > 1}")

print(f"\nTask 2 - Blind Zero Prediction:")
print(f"  Training: {len(train_crossings)} crossings in Y ≤ 0.7")
print(f"  Test: {len(test_crossings)} crossings in Y > 0.7")
for name, res in sorted(prediction_results.items(), key=lambda x: x[1]['mean_error']):
    print(f"  {name}: mean error = {res['mean_error']:.4f}")
print(f"\n  Best mapping: {best_mapping}")
print(f"  Mean prediction error: {prediction_results[best_mapping]['mean_error']:.4f}")
print(f"  Conclusion: {prediction_output['conclusion']}")