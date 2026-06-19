#!/usr/bin/env python3
import numpy as np, pandas as pd, json, os, sys
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_and_filter(path, minlen=3):
    df = pd.read_csv(path)
    counts = df.groupby('traj_id')['frame'].transform('count')
    dfc = df[counts >= minlen].copy()
    dfc = dfc.dropna(subset=['x','y'])
    return dfc

def loss_abs(params, y, xmag):
    k, A, phase = params
    pred = np.abs(A * np.cos(k * y + phase))
    return np.mean((xmag - pred)**2)

def fit_abs_model(y, xmag, init=None, bounds=None):
    if bounds is None:
        bounds = [(0.01,2.0),(0.01,0.6),(-np.pi,np.pi)]
    if init is None:
        init = [0.63, np.clip(np.std(xmag)*1.2, 0.01, 0.6), 0.0]
    res = minimize(lambda p: loss_abs(p, y, xmag), x0=init, bounds=bounds, method='L-BFGS-B')
    return res

def main(csv_path, out_dir, B=200, N=500):
    os.makedirs(out_dir, exist_ok=True)
    df = load_and_filter(csv_path)
    if df.empty:
        print('No points after filtering (min track len=3)'); return 1
    y = df['y'].values.astype(float)
    x = df['x'].values.astype(float)
    xmag = np.abs(x - 0.5)
    res = fit_abs_model(y, xmag)
    k_opt, A_opt, ph_opt = res.x
    rms_H1 = float(np.sqrt(res.fun))
    rms_H0 = float(np.sqrt(np.mean(xmag**2)))
    rng = np.random.default_rng(12345)
    boot = []
    for b in range(B):
        idx = rng.integers(0, len(y), size=len(y))
        yb, xm = y[idx], xmag[idx]
        try:
            rb = fit_abs_model(yb, xm)
            boot.append(rb.x)
        except:
            continue
    boot = np.array(boot)
    if boot.size>0:
        k_ci = np.percentile(boot[:,0],[2.5,97.5]).tolist()
        A_ci = np.percentile(boot[:,1],[2.5,97.5]).tolist()
        ph_ci = np.percentile(boot[:,2],[2.5,97.5]).tolist()
    else:
        k_ci=A_ci=ph_ci=[None,None]
    rms_perm = []
    for i in range(N):
        xm_sh = rng.permutation(xmag)
        r = fit_abs_model(y, xm_sh)
        rms_perm.append(float(np.sqrt(r.fun)))
    pval = float(np.mean(np.array(rms_perm) <= rms_H1))
    out = {"n_points": int(len(y)), "fitted":{"k":float(k_opt),"A":float(A_opt),"phase":float(ph_opt),"rms_H1":rms_H1}, "null_vertical":{"rms_H0":rms_H0}, "bootstrap_CI":{"k_CI":k_ci,"A_CI":A_ci,"phase_CI":ph_ci}, "permutation":{"n_perm":int(N),"p_rms_leq_obs":pval}}
    with open(os.path.join(out_dir,'manifold_summary.json'),'w') as fh:
        json.dump(out, fh, indent=2)
    # plots
    try:
        plt.figure(figsize=(6,4))
        plt.scatter(xmag,y,s=6,alpha=0.6)
        ygrid = np.linspace(y.min(),y.max(),1000)
        plt.plot(np.abs(A_opt*np.cos(k_opt*ygrid+ph_opt)), ygrid, color='red', linewidth=2)
        plt.gca().invert_yaxis()
        plt.xlabel('|x-0.5|'); plt.ylabel('y')
        plt.tight_layout(); plt.savefig(os.path.join(out_dir,'manifold_fit_mag.png'))
        plt.close()
        plt.figure(figsize=(5,3))
        plt.hist(rms_perm,bins=50,alpha=0.7); plt.axvline(rms_H1,color='red',linestyle='--')
        plt.tight_layout(); plt.savefig(os.path.join(out_dir,'manifold_perm_rms.png'))
        plt.close()
    except Exception as e:
        print('plot error', e)
    print(json.dumps(out, indent=2))
    return 0

if __name__=='__main__':
    if len(sys.argv)<3:
        print('usage: dsl_manifold_test.py <csv> <outdir>'); sys.exit(1)
    main(sys.argv[1], sys.argv[2])
