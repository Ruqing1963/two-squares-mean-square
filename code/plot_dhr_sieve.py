r"""
Figure 4 of Part II: the two-dimensional sieve against the vector sieve (Section 7).
(a) F_2, f_2 (dhr_sieve.py) with alpha_2, beta_2, the point s = 9 of Theorem 7.1 and the range [1.8, 8] of the
    integral R_2^*; dashed: the symmetric vector sieve, upper F_1(s/2)^2 and lower 2 f_1(s/2) F_1(s/2) - F_1(s/2)^2,
    positive only for s > 2 + 2 sqrt(e).
(b) The best attainable 1/lam + 2/b and 1/lam + 1/b against the level theta, for both sieves (dimension 2 computed
    here, cached in results/dhr_thresholds.txt; vector sieve from results/sieve_thresholds.txt).
Usage: python plot_dhr_sieve.py [--recompute]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import exp, sqrt
from dhr_sieve import DHR, scan
from sieve_core import F_lin, f_lin

BLUE, ORANGE, AQUA, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#898781", "#e1e0d9"
D2 = DHR(2, verbose=False)

# ---- (b): thresholds of the dimension-2 sieve against theta (cached) ----
cache = "../results/dhr_thresholds.txt"
done = {}
if os.path.exists(cache) and "--recompute" not in sys.argv:
    for line in open(cache):
        if line[0].isdigit():
            th, jt, sd = (float(x) for x in line.split()[:3])
            done[round(th, 4)] = (jt, sd)
thetas = [round(t, 4) for t in np.arange(0.40, 0.991, 0.03)] + [0.5, 2 / 3, 0.95]
thetas = sorted(set(round(t, 4) for t in thetas))
with open(cache, "a" if done else "w") as fh:
    if not done:
        fh.write("theta  min(1/lam+2/b)  min(1/lam+1/b)   (dimension-2 sieve, floating point, dhr_sieve.scan)\n")
    for th in thetas:
        if th in done:
            continue
        rows = scan(D2, th)
        jt = min(r[0] for r in rows)
        sd = min(r[1] for r in rows)
        done[th] = (jt, sd)
        fh.write(f"{th:.4f} {jt:.4f} {sd:.4f}\n")
        fh.flush()
        print(f"theta = {th:.4f}: joint {jt:.4f}, side {sd:.4f}", flush=True)
vec = {}
for line in open("../results/sieve_thresholds.txt"):
    parts = line.split()
    if parts and parts[0][0].isdigit():
        vec[round(float(parts[0]), 4)] = (float(parts[1]), float(parts[2]))

fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.3))
# ---- (a) ----
s = np.linspace(1.5, 10, 900)
F2, f2 = D2.Fs(s), D2.fs(s)
ax.axvspan(1.8, 8, color=GRID, alpha=0.55, lw=0, label=r"range of the integral $R_2^*$")
ax.plot(s, F2, color=BLUE, lw=2, label=r"$F_2(s)$")
ax.plot(s, f2, color=ORANGE, lw=2, label=r"$f_2(s)$")
Fv = F_lin(s / 2) ** 2
fv = 2 * f_lin(s / 2) * F_lin(s / 2) - F_lin(s / 2) ** 2
ax.plot(s, Fv, color=BLUE, lw=1.3, ls="--", label=r"vector sieve: $F_1(s/2)^2$")
ax.plot(s, np.maximum(fv, 0), color=ORANGE, lw=1.3, ls="--", label=r"vector sieve: $2f_1F_1-F_1^2$ at $s/2$")
ax.axhline(1, color=MUTED, lw=0.8)
for xv, lab, dy in ((D2.beta, r"$\beta_2=4.27$", 0.06), (2 + 2 * sqrt(exp(1)), r"$2+2\sqrt{e}=5.30$", 0.06),
                    (D2.alpha, r"$\alpha_2=5.36$", 1.42), (9, r"$s=9$", 0.06)):
    ax.axvline(xv, color=MUTED, lw=0.8, ls=":")
    ax.text(xv + 0.05, dy, lab, fontsize=7.5, color=INK, rotation=90, va="bottom")
ax.plot([9], [D2.fs(9)], "o", color=ORANGE, ms=5)
ax.text(9.1, 0.93, r"$f_2(9)=0.9997$", fontsize=7.5, ha="right")
ax.set_xlim(1.5, 10)
ax.set_ylim(0, 2.2)
ax.set_xlabel(r"$s$")
ax.set_title("(a) sieve functions of dimension 2 and the symmetric vector sieve", fontsize=10)
ax.legend(fontsize=7.5, loc="upper right", framealpha=0.95)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
# ---- (b) ----
tv = sorted(vec)
td = sorted(done)
bx.plot(tv, [vec[t][0] for t in tv], color=BLUE, lw=1.3, ls="--", label=r"vector sieve, $\min(1/\lambda+2/b)$")
bx.plot(td, [done[t][0] for t in td], color=BLUE, lw=2, marker="o", ms=3, label=r"dimension 2, $\min(1/\lambda+2/b)$")
bx.plot(tv, [vec[t][1] for t in tv], color=ORANGE, lw=1.3, ls="--", label=r"vector sieve, $\min(1/\lambda+1/b)$")
bx.plot(td, [done[t][1] for t in td], color=ORANGE, lw=2, marker="o", ms=3, label=r"dimension 2, $\min(1/\lambda+1/b)$")
for k in (4, 5, 6, 7, 8, 9):
    bx.axhline(k, color=GRID, lw=0.8)
    bx.text(0.995, k + 0.05, f"{k}", fontsize=7, color=MUTED, ha="right")
for th, lab in ((0.5, r"$\theta=1/2$"), (2 / 3, r"$2/3$"), (0.95, r"$19/20$")):
    bx.axvline(th, color=MUTED, lw=0.8, ls=":")
    bx.text(th + 0.005, 3.08, lab, fontsize=7.5, color=INK)
bx.annotate(r"8.03 $\rightarrow$ 7.56: $\Omega\leq7$", (0.5, done[0.5][0]), xytext=(0.56, 8.6), fontsize=7.5,
            arrowprops=dict(arrowstyle="-", color=INK, lw=0.6))
bx.set_xlim(0.38, 1.0)
bx.set_ylim(3, 9.7)
bx.set_xlabel(r"level of distribution $\theta$")
bx.set_ylabel("threshold")
bx.set_title(r"(b) best thresholds against the level ($\Omega\leq k$ needs $\leq k+1$)", fontsize=10)
bx.legend(fontsize=7.5, loc="upper right", framealpha=0.95)
for sp in ("top", "right"):
    bx.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig("../figures/dhr_sieve.pdf")
fig.savefig("../figures/dhr_sieve.png", dpi=160)
print("written figures/dhr_sieve.pdf/.png")
