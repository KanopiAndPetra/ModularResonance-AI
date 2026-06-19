#!/usr/bin/env python3
"""
build_full_von_mangoldt.py
Von Mangoldt animation pipeline.
"""

import os, math, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import imageio
from scipy import ndimage

FRAME_DIR_VM = "frames_vm"
FRAME_DIR_LOG = "frames_log"
OUT_VM_VIDEO = "von_mangoldt_animation.mp4"
OUT_LOG_VIDEO = "logdriver_animation.mp4"
TRAJ_VM_CSV = "traj_vm.csv"
TRAJ_LOG_CSV = "traj_log.csv"

NUM_FRAMES = 10
DELTA_T = 1.0
RES = 300
Y_MAX = 60.0
FPS = 3
LOCAL_MAX_PCT = 99.5
MIN_REGION_AREA = 1
TRACK_MAX_DIST = 0.03

os.makedirs(FRAME_DIR_VM, exist_ok=True)
os.makedirs(FRAME_DIR_LOG, exist_ok=True)

def S_314_von_mangoldt_vec(Y):
    T = np.abs(Y)
    out = np.ones_like(T, dtype=float)
    mask = T >= 1.0
    if np.any(mask):
        TT = T[mask]
        idx = (TT / (2.0 * np.pi)) * np.log(np.maximum(TT / (2.0 * np.pi), 1e-12))
        out[mask] = 2.5 + 1.5 * np.sin(2.0 * np.pi * idx)
    return out

def Psi_vm_vec(X, Y):
    s_mag = S_314_von_mangoldt_vec(Y)
    Zmag = np.sqrt(X**2 + Y**2)
    return s_mag * Zmag * np.sin(np.pi * X)

def S_314_logdriver_vec(Y):
    T = np.abs(Y)
    out = np.ones_like(T, dtype=float)
    mask = T >= 1.0
    if np.any(mask):
        out[mask] = 2.5 + 1.5 * np.sin(2.0 * np.pi * np.log(T[mask] + 1.0))
    return out

def Psi_log_vec(X, Y):
    s_mag = S_314_logdriver_vec(Y)
    Zmag = np.sqrt(X**2 + Y**2)
    return s_mag * Zmag * np.sin(np.pi * X)

def generate_heatmap(Psi_vec_func, T_shift, y_max=Y_MAX, res=RES):
    x = np.linspace(0.0, 1.0, res)
    y = np.linspace(0.0 + T_shift, y_max + T_shift, res)
    X, Y = np.meshgrid(x, y)
    Energy_Pi = np.pi * np.sqrt(X**2 + Y**2)
    Energy_Psi = Psi_vec_func(X, Y)
    Diff = np.abs(Energy_Pi - Energy_Psi)
    InvDiff = 1.0 / (Diff + 0.1)
    return X, Y, InvDiff

def local_maxima_coords(arr, pct=LOCAL_MAX_PCT):
    k = 3
    pad = k // 2
    padded = np.pad(arr, pad, mode='edge')
    sm = np.zeros_like(arr)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            sm[i, j] = padded[i:i+k, j:j+k].mean()
    thr = np.percentile(sm, pct)
    padded2 = np.pad(sm, 1, mode='constant', constant_values=-np.inf)
    is_max = np.ones_like(sm, dtype=bool)
    M, N = sm.shape
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            neighbor = padded2[1+di:M+1+di, 1+dj:N+1+dj]
            is_max &= (sm > neighbor)
    ys, xs = np.where(is_max & (sm >= thr))
    return xs, ys, sm[ys, xs]

def track_trajectories(peaks_list, max_dist=TRACK_MAX_DIST):
    trajectories = []
    if len(peaks_list) == 0:
        return trajectories
    first = peaks_list[0]
    for p in first:
        trajectories.append([p])
    for frame_idx in range(1, len(peaks_list)):
        current = peaks_list[frame_idx]
        used = [False] * len(current)
        for traj in trajectories:
            last = traj[-1]
            best_i = -1; best_d = None
            if last is None:
                traj.append(None)
                continue
            for i, c in enumerate(current):
                if used[i]:
                    continue
                d = math.hypot(last[0] - c[0], last[1] - c[1])
                if best_d is None or d < best_d:
                    best_d = d; best_i = i
            if best_d is not None and best_d <= max_dist:
                traj.append(current[best_i])
                used[best_i] = True
            else:
                traj.append(None)
        for i, c in enumerate(current):
            if not used[i]:
                new = [None] * frame_idx
                new.append(c)
                trajectories.append(new)
    return trajectories

def save_trajs(trajs, path):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['traj_id', 'frame', 'x', 'y'])
        for tid, t in enumerate(trajs):
            for fi, p in enumerate(t):
                if p is None:
                    w.writerow([tid, fi, '', ''])
                else:
                    w.writerow([tid, fi, p[0], p[1]])

def compute_locking_metrics(trajectories, x_center=0.5, x_tol=0.08):
    fractions = []
    x_vars = []
    for traj in trajectories:
        xs = [p[0] for p in traj if p is not None]
        if len(xs) == 0:
            continue
        frac = sum(1 for x in xs if abs(x - x_center) <= x_tol) / len(xs)
        fractions.append(frac)
        x_vars.append(np.var(xs))
    return {
        'num_traj': len(trajectories),
        'mean_fraction_locked': float(np.mean(fractions)) if fractions else 0.0,
        'median_fraction_locked': float(np.median(fractions)) if fractions else 0.0,
        'mean_x_variance': float(np.mean(x_vars)) if x_vars else float('nan')
    }

def run_pipeline():
    peaks_vm = []
    peaks_log = []
    for i in range(NUM_FRAMES):
        T_shift = i * DELTA_T

        X, Y, Heat_vm = generate_heatmap(Psi_vm_vec, T_shift, y_max=Y_MAX, res=RES)
        xs, ys, vals = local_maxima_coords(Heat_vm, pct=LOCAL_MAX_PCT)
        x_coords = X[ys, xs]
        y_coords = Y[ys, xs] - T_shift
        peaks_vm.append(list(zip(x_coords.tolist(), y_coords.tolist())))

        plt.figure(figsize=(5, 8))
        plt.pcolormesh(X, Y - T_shift, Heat_vm, shading='auto', cmap='plasma')
        plt.scatter(x_coords, y_coords, s=10, facecolors='none', edgecolors='white')
        plt.axvline(x=0.5, color='cyan', linewidth=1)
        plt.title(f"VM Frame {i}")
        plt.tight_layout()
        plt.savefig(f"{FRAME_DIR_VM}/vm_frame_{i:03d}.png", dpi=120)
        plt.close()

        Xl, Yl, Heat_log = generate_heatmap(Psi_log_vec, T_shift, y_max=Y_MAX, res=RES)
        xs2, ys2, vals2 = local_maxima_coords(Heat_log, pct=LOCAL_MAX_PCT)
        x_coords2 = Xl[ys2, xs2]
        y_coords2 = Yl[ys2, xs2] - T_shift
        peaks_log.append(list(zip(x_coords2.tolist(), y_coords2.tolist())))

        plt.figure(figsize=(5, 8))
        plt.pcolormesh(Xl, Yl - T_shift, Heat_log, shading='auto', cmap='plasma')
        plt.scatter(x_coords2, y_coords2, s=10, facecolors='none', edgecolors='white')
        plt.axvline(x=0.5, color='cyan', linewidth=1)
        plt.title(f"LOG Frame {i}")
        plt.tight_layout()
        plt.savefig(f"{FRAME_DIR_LOG}/log_frame_{i:03d}.png", dpi=120)
        plt.close()

    traj_vm = track_trajectories(peaks_vm, max_dist=TRACK_MAX_DIST)
    traj_log = track_trajectories(peaks_log, max_dist=TRACK_MAX_DIST)

    save_trajs(traj_vm, TRAJ_VM_CSV)
    save_trajs(traj_log, TRAJ_LOG_CSV)

    files_vm = sorted(f for f in os.listdir(FRAME_DIR_VM) if f.endswith(".png"))
    frames_vm = [imageio.v3.imread(os.path.join(FRAME_DIR_VM, f)) for f in files_vm]
    imageio.mimsave(OUT_VM_VIDEO, frames_vm, fps=FPS)

    files_log = sorted(f for f in os.listdir(FRAME_DIR_LOG) if f.endswith(".png"))
    frames_log = [imageio.v3.imread(os.path.join(FRAME_DIR_LOG, f)) for f in files_log]
    imageio.mimsave(OUT_LOG_VIDEO, frames_log, fps=FPS)

    metrics_vm = compute_locking_metrics(traj_vm, x_center=0.5, x_tol=0.08)
    metrics_log = compute_locking_metrics(traj_log, x_center=0.5, x_tol=0.08)

    print("Done. Outputs:")
    print(" -", OUT_VM_VIDEO)
    print(" -", OUT_LOG_VIDEO)
    print(" -", FRAME_DIR_VM, "and", FRAME_DIR_LOG, "(png frames)")
    print(" -", TRAJ_VM_CSV, ",", TRAJ_LOG_CSV)
    print("Metrics (VM):", metrics_vm)
    print("Metrics (LOG):", metrics_log)

if __name__ == "__main__":
    run_pipeline()
