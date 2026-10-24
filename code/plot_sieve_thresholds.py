"""
Figure 2 of Part II: the best values of 1/lam + 2/beta1 (which bounds Omega(n(N-n))) and of 1/lam + 1/beta1
(which bounds max(Omega(n), Omega(N-n))) attainable with L - lam R > 0, as functions of the level of distribution
theta, for the vector sieve with Richert weights of Part I (floating-point model of modular_sieve_scan.py).
For each theta the two quantities are minimised separately over alpha, beta1 (lam = L/R, the largest admissible
value); the integer thresholds k + 1 give the bounds Omega(n(N-n)) <= k and max Omega <= k'.

Floating point, not certified; the certified points are those of Part I (theta = 2/5, 1/2) and of this paper
(theta = 2/3 and 19/20).  Usage:  python plot_sieve_thresholds.py   (theta values missing from the cache
results/sieve_thresholds.txt are computed and appended; add --recompute to redo the whole scan)
"""
import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from modular_sieve_scan import best_L, R_value

DATA = "../results/sieve_thresholds.txt"
thetas = np.round(np.arange(0.38, 0.9901, 0.01), 4)
alphas = np.arange(0.020, 0.1001, 0.004)
done = {}                                                     # theta -> (joint, side), from the cache
if os.path.exists(DATA) and "--recompute" not in sys.argv:
    for th, bj, bs in np.loadtxt(DATA, ndmin=2):
        done[round(float(th), 4)] = (bj, bs)
todo = [th for th in thetas if float(th) not in done]          # only the theta values missing from the cache
t0 = time.time()
for th in todo:
    bj, bs = np.inf, np.inf
    for a in alphas:
        L, _ = best_L(a, th)
        if L <= 0:
            continue
        for b1 in np.arange(0.20, th - 0.005, 0.01):
            if b1 <= a:
                continue
            R = R_value(a, th, b1, nu=100)
            lam = L / R
            bj = min(bj, 1 / lam + 2 / b1)
            bs = min(bs, 1 / lam + 1 / b1)
    done[float(th)] = (bj, bs)
    print(f"theta={th:.2f}: min 1/lam+2/beta1 = {bj:.3f}, min 1/lam+1/beta1 = {bs:.3f}  ({time.time() - t0:.0f}s)", flush=True)
cached = len(todo) == 0
joint = np.array([done[float(th)][0] for th in thetas])
side = np.array([done[float(th)][1] for th in thetas])
if not cached:
    np.savetxt(DATA, np.column_stack([thetas, joint, side]), fmt="%.4f",
               header="theta  min(1/lam+2/beta1)  min(1/lam+1/beta1)   (floating point, modular_sieve_scan model)")

fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot(thetas, joint, color="#2a78d6", lw=1.8, label=r"$\min\,(1/\lambda+2/\beta_1)$:  $\Omega(n(N-n))\leq k$ if below $k+1$")
ax.plot(thetas, side, color="#eb6834", lw=1.8, label=r"$\min\,(1/\lambda+1/\beta_1)$:  $\max\Omega\leq k'$ if below $k'+1$")
for k in (10, 9, 8, 7, 6, 5, 4):
    ax.axhline(k, color="#999", lw=0.7, ls=":")
for th, lab in ((2 / 5, r"$2/5$: $(9,6)$"), (1 / 2, r"$1/2$: $(8,5)$"),
                (2 / 3, r"$2/3$: $(6,4)$"), (19 / 20, r"$19/20$: $(5,4)$")):
    ax.axvline(th, color="#bbb", lw=0.8)
    ax.text(th + 0.004, 3.1, lab, fontsize=7.5, rotation=90, va="bottom", color="#333")
ax.set_xlim(0.38, 1.0)
ax.set_ylim(3.0, 10)
ax.set_yticks(range(3, 11))
ax.set_xlabel(r"level of distribution $\theta$")
ax.set_ylabel("best attainable value (dotted: integer thresholds)")
ax.grid(True, axis="x", alpha=0.2)
ax.legend(fontsize=7.5, loc="upper right")
fig.tight_layout()
fig.savefig("../figures/sieve_thresholds.pdf")
fig.savefig("../figures/sieve_thresholds.png", dpi=150)
print("wrote ../figures/sieve_thresholds.pdf/.png" + (" (scan read from cache)" if cached else " and " + DATA))
