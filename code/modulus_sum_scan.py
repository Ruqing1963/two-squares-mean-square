r"""
Sums over the moduli of the cuspidal coefficients a_{G_{d,e}}(N) = T(d,e;N) - M(d,e;N) at a FIXED N (Part II,
Section 8).  The remainder of the sieve at the level theta is sum_{de <= N^theta} lambda_d lambda'_e R(d,e;N), and
R(d,e;N) is, up to the Moebius inversion of (n,N) = 1, the coefficient a_{G_{d,e}}(N).  Part I bounds the sum of the
absolute values, which limits theta to 1/2; going beyond needs cancellation between the signs of a_{G_{d,e}}(N) as
(d,e) varies.  This script measures, for several N = 2 mod 4 and D = N^theta, over the coprime odd square-free pairs
(d,e) with de <= D:
    A = sum |a|,   S_1 = sum a,   S_mu = sum mu(de) a,   S_rnd = sum (random sign) a,   Q = (sum a^2)^(1/2),
so that |S|/Q = 1 is square-root cancellation and A/Q the maximum possible.  Main term exact for every N
(local_densities.py formulae; beta_2 depends on (d, e) mod 4 only).  Floating point.
Usage: python modulus_sum_scan.py [N ...]
"""
import sys
import time
import numpy as np
from math import pi, gcd, log
from sup_error_scan import r_table, factor, delta2_table

Ns = [int(a) for a in sys.argv[1:]] or [2000002, 2000006, 2000018, 3000002]
THETAS = (0.5, 0.6, 0.7, 0.75)
rng = np.random.default_rng(1)


def mobius_sqfree_odd(D):
    """mu(k) for k <= D, and a flag for odd square-free k."""
    mu = np.ones(D + 1, dtype=np.int64)
    sqf = np.ones(D + 1, dtype=bool)
    p = 2
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
    """prod over odd p of beta_p(Q_{d,e}; N), with the exact formulae at p | de."""
    prod = 8 / pi ** 2        # prod_{p odd} (1 - p^-2)
    fN = factor(N)
    for p, k in fN.items():
        if p == 2:
            continue
        if (d * e) % p:
            prod *= (1 - p ** -(k + 1.0)) / (1 - 1.0 / p)      # sigma_{-1}(p^k); the (1-p^-2) is in 8/pi^2
        else:
            prod /= 1 - p ** -2.0
            prod *= (1 + 1 / p) * p ** -k if p % 4 == 3 else 2 - (1 + 1 / p) * p ** -k
    for p in factor(d * e):
        if N % p:
            prod /= 1 - p ** -2.0
            prod *= 1 - (1 if p % 4 == 1 else -1) / p
    return prod


for N in Ns:
    assert N % 4 == 2
    t0 = time.time()
    Dmax = int(N ** THETAS[-1]) + 1
    r = r_table(N)
    mu, sqf = mobius_sqfree_odd(Dmax)
    d2 = {(a, b): delta2_table(a, b) for a in (1, 3) for b in (1, 3)}
    bp_cache = {}
    pairs = []
    for d in range(1, Dmax + 1, 2):
        if not sqf[d]:
            continue
        for e in range(1, Dmax // d + 1, 2):
            if not sqf[e] or gcd(d, e) != 1:
                continue
            # T(d,e;N) = sum_{alpha = alpha_0 mod e} r(alpha) r((N - d alpha)/e),  alpha_0 = N d^{-1} mod e
            al = np.arange((N * pow(d, -1, e)) % e if e > 1 else 0, N // d + 1, e)
            T = float(np.sum(r[al] * r[(N - d * al) // e]))
            M = pi ** 2 * N / (d * e) * beta_odd_primes(d, e, N) * d2[(d % 4, e % 4)][N % 512]
            pairs.append((d * e, mu[d * e], T - M))
    pairs.sort()
    k = np.array([p[0] for p in pairs], dtype=float)
    m = np.array([p[1] for p in pairs], dtype=float)
    a = np.array([p[2] for p in pairs])
    sgn = rng.choice([-1.0, 1.0], size=len(a))
    print(f"N = {N} = {'*'.join(f'{p}^{v}' if v > 1 else str(p) for p, v in sorted(factor(N).items()))}, "
          f"{len(pairs)} pairs with de <= N^{THETAS[-1]} = {Dmax}  [{time.time() - t0:.0f}s]")
    print(f"   rms of a (de)^(1/2) / N^(1/2) over all pairs: {np.sqrt(np.mean(a ** 2 * k)) / N ** 0.5:.3f}")
    print(f"   {'theta':>6} {'D':>7} {'#pairs':>7} {'A=sum|a|':>10} {'Q=(sum a^2)^1/2':>16} {'S_1/Q':>8} {'S_mu/Q':>8} {'S_rnd/Q':>8} {'A/Q':>7} {'A/N':>8}")
    for th in THETAS:
        D = N ** th
        sel = k <= D
        A = np.sum(np.abs(a[sel]))
        Q = np.sqrt(np.sum(a[sel] ** 2))
        S1 = np.sum(a[sel])
        Smu = np.sum(m[sel] * a[sel])
        Sr = np.sum(sgn[sel] * a[sel])
        print(f"   {th:>6.2f} {D:>7.0f} {sel.sum():>7} {A:>10.0f} {Q:>16.1f} {S1 / Q:>8.2f} {Smu / Q:>8.2f} {Sr / Q:>8.2f} {A / Q:>7.1f} {A / N:>8.4f}")
    sys.stdout.flush()
