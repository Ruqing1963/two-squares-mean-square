r"""
Moments and extreme values of the normalised coefficients u(N) = a_G(N)/sqrt(N), G = G_{m,m'}, over the admissible
N <= X of Part I (N coprime to mm', v_2(N) <= 3), with the exact main term of local_densities.py (Part II, Section 8.3).

Prints, for each pair: the root mean square, the kurtosis m_4/m_2^2 (3 for a Gaussian), the ratio sup/rms and the
value sqrt(2 log n_samples) a Gaussian sample of the same size would give, the share of the fourth moment carried by
the N with d(N) >= 24, the conditional mean of u(N)^2/m_2 given d(N) in dyadic ranges, and the ten largest |u(N)|
with d(N) and omega(N).   Floating point.   Usage: python coefficient_moments.py [X] [m,m' ...]
"""
import sys
import time
import numpy as np
from math import gcd, log, sqrt
from sup_error_scan import r_table, sigma_over_n_odd, factor
from local_densities import siegel_all

args = [a for a in sys.argv[1:] if not a.startswith('--')]
X = int(args[0]) if args and ',' not in args[0] else 2 ** 21
pairs = [tuple(int(t) for t in a.split(',')) for a in args if ',' in a] or \
    [(1, 5), (3, 5), (1, 101), (13, 29), (1, 1009), (61, 67), (101, 103)]

t0 = time.time()
r = r_table(X)
f_sigma, v2 = sigma_over_n_odd(X)
# divisor function and omega by sieving
d = np.ones(X + 1, dtype=np.int32)
om = np.zeros(X + 1, dtype=np.int32)
for k in range(2, X + 1):
    d[k::k] += 1
sieve = np.ones(X + 1, dtype=bool)
sieve[:2] = False
for p in range(2, int(X ** 0.5) + 1):
    if sieve[p]:
        sieve[p * p::p] = False
for p in np.nonzero(sieve)[0]:
    om[p::p] += 1
print(f"tables to {X} in {time.time() - t0:.1f}s\n")

for (m, mp) in pairs:
    assert gcd(m, mp) == 1 and m % 2 == 1 and mp % 2 == 1
    t1 = time.time()
    T, M = siegel_all(m, mp, X, r, f_sigma, v2)
    n = np.arange(X + 1)
    ok = (n >= 1) & (v2 <= 3)
    for p in factor(m * mp):
        ok &= n % p != 0
    ok[0] = False
    u = (T[ok] - M[ok]) / np.sqrt(n[ok])
    dN, oN, NN = d[ok], om[ok], n[ok]
    m2 = np.mean(u ** 2)
    m4 = np.mean(u ** 4)
    rms = sqrt(m2)
    ns = len(u)
    isup = np.argmax(np.abs(u))
    big = dN >= 24
    share4 = np.sum(u[big] ** 4) / np.sum(u ** 4)
    share2 = np.sum(u[big] ** 2) / np.sum(u ** 2)
    frac_big = np.mean(big)
    print(f"({m},{mp}): mm' = {m * mp}, {ns} admissible N <= {X}:  rms = {rms:.4f},  kurtosis m4/m2^2 = {m4 / m2 ** 2:.3f},"
          f"  sup/rms = {abs(u[isup]) / rms:.2f}  (Gaussian sample of this size: {sqrt(2 * log(ns)):.2f}),"
          f"  share of N with d(N) >= 24: {frac_big:.3f} of the N, {share2:.3f} of m2, {share4:.3f} of m4")
    print("   E[u^2]/m2 given d(N) in:", end="")
    for lo, hi in ((1, 4), (4, 8), (8, 16), (16, 32), (32, 64), (64, 1 << 30)):
        sel = (dN >= lo) & (dN < hi)
        if sel.any():
            print(f"  [{lo},{hi}): {np.mean(u[sel] ** 2) / m2:.2f} ({sel.sum()})", end="")
    print()
    idx = np.argsort(-np.abs(u))[:10]
    print("   largest |u(N)|/rms:", ", ".join(f"N={NN[i]} ({abs(u[i]) / rms:.1f}; d={dN[i]}, omega={oN[i]})" for i in idx))
    print(f"   [{time.time() - t1:.0f}s]")
    sys.stdout.flush()
