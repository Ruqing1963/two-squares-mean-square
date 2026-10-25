r"""
Figure 5 of Part II: cancellation in the sum over the moduli (Section 8.2), for N = 3000002.
(a) A/N, Q/N, |S_1|/N, |S_mu|/N as functions of theta (sums over the coprime odd square-free pairs with de <= N^theta).
(b) The partial sums sum_{de <= M} mu(de) a_{d,e} / Q(M) against M, with the band |.| <= 1.5.
Recomputes the coefficients as in modulus_sum_scan.py (about 10 s).   Usage: python plot_modulus_sums.py [N]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import pi, gcd, log
from sup_error_scan import r_table, factor, delta2_table


def mobius_sqfree_odd(D):
    mu = np.ones(D + 1, dtype=np.int64)
    sqf = np.ones(D + 1, dtype=bool)
    isp = np.ones(D + 1, dtype=bool)
    isp[:2] = False
    for p in range(2, D + 1):
        if not isp[p]:
            continue
        isp[p * p::p] = False
        mu[p::p] *= -1
        sqf[p * p::p * p] = False
    mu[~sqf] = 0
    mu[::2] = 0
    sqf[::2] = False
    sqf[0] = False
    return mu, sqf


def beta_odd_primes(d, e, N):
    prod = 8 / pi ** 2
    for p, k in factor(N).items():
        if p == 2:
            continue
        if (d * e) % p:
            prod *= (1 - p ** -(k + 1.0)) / (1 - 1.0 / p)
        else:
            prod /= 1 - p ** -2.0
            prod *= (1 + 1 / p) * p ** -k if p % 4 == 3 else 2 - (1 + 1 / p) * p ** -k
    for p in factor(d * e):
        if N % p:
            prod /= 1 - p ** -2.0
            prod *= 1 - (1 if p % 4 == 1 else -1) / p
    return prod


BLUE, ORANGE, AQUA, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#898781", "#e1e0d9"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000002
TH_MAX = 0.75
Dmax = int(N ** TH_MAX) + 1
r = r_table(N)
mu, sqf = mobius_sqfree_odd(Dmax)
d2 = {(a, b): delta2_table(a, b) for a in (1, 3) for b in (1, 3)}
pairs = []
for d in range(1, Dmax + 1, 2):
    if not sqf[d]:
        continue
    for e in range(1, Dmax // d + 1, 2):
        if not sqf[e] or gcd(d, e) != 1:
            continue
        al = np.arange((N * pow(d, -1, e)) % e if e > 1 else 0, N // d + 1, e)
        T = float(np.sum(r[al] * r[(N - d * al) // e]))
        M = pi ** 2 * N / (d * e) * beta_odd_primes(d, e, N) * d2[(d % 4, e % 4)][N % 512]
        pairs.append((d * e, mu[d * e], T - M))
pairs.sort()
k = np.array([p[0] for p in pairs], dtype=float)
m = np.array([p[1] for p in pairs], dtype=float)
a = np.array([p[2] for p in pairs])
cA, cQ2, c1, cmu = np.cumsum(np.abs(a)), np.cumsum(a ** 2), np.cumsum(a), np.cumsum(m * a)
thetas = np.linspace(0.2, TH_MAX, 120)
idx = np.searchsorted(k, N ** thetas, side="right") - 1

fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.2))
ax.plot(thetas, cA[idx] / N, color=ORANGE, lw=2, label=r"$A/N=\sum|a_{d,e}|/N$")
ax.plot(thetas, np.sqrt(cQ2[idx]) / N, color=INK, lw=2, label=r"$Q/N=(\sum a_{d,e}^2)^{1/2}/N$")
ax.plot(thetas, np.abs(c1[idx]) / N, color=BLUE, lw=1.4, label=r"$|S_1|/N=|\sum a_{d,e}|/N$")
ax.plot(thetas, np.abs(cmu[idx]) / N, color=AQUA, lw=1.4, label=r"$|S_\mu|/N=|\sum\mu(de)a_{d,e}|/N$")
ax.axhline(1, color=MUTED, lw=0.8, ls="--")
ax.text(0.745, 1.15, "sum = N", fontsize=7.5, color=MUTED, ha="right")
ax.axvline(0.5, color=MUTED, lw=0.8, ls=":")
ax.text(0.505, 2e-4, r"$\theta=1/2$", fontsize=7.5, color=INK)
ax.set_yscale("log")
ax.set_xlim(0.2, TH_MAX)
ax.set_ylim(1e-4, 30)
ax.set_xlabel(r"$\theta$  ($de\leq N^\theta$)")
ax.set_ylabel("sum / N")
ax.set_title(f"(a) sums of $a_{{G_{{d,e}}}}(N)$ over the moduli, $N={N}$", fontsize=10)
ax.legend(fontsize=7.8, loc="upper left", framealpha=0.95)
ax.grid(True, which="major", color=GRID, lw=0.6)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
# (b)
Q = np.sqrt(cQ2)
bx.fill_between(k, -1.5, 1.5, color=GRID, alpha=0.5, lw=0)
bx.plot(k, cmu / Q, color=AQUA, lw=0.8, label=r"$\sum_{de\leq M}\mu(de)\,a_{d,e}\,/\,Q(M)$")
bx.plot(k, c1 / Q, color=BLUE, lw=0.8, alpha=0.8, label=r"$\sum_{de\leq M}a_{d,e}\,/\,Q(M)$")
bx.axhline(0, color=INK, lw=0.6)
bx.set_xscale("log")
bx.set_xlim(10, k[-1])
bx.set_ylim(-4, 4)
bx.set_xlabel(r"$M$")
bx.set_ylabel("partial sum / Q(M)")
bx.set_title(r"(b) normalised partial sums; band $|\cdot|\leq1.5$", fontsize=10)
bx.legend(fontsize=7.8, loc="upper left", framealpha=0.95)
for sp in ("top", "right"):
    bx.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig("../figures/modulus_sums.pdf")
fig.savefig("../figures/modulus_sums.png", dpi=160)
print("written figures/modulus_sums.pdf/.png;", len(pairs), "pairs")
