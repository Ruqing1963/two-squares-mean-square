r"""
The sandwich for the norm (Part II, Section 6.7):

    Lambda_1(m,m')  <=  ||G_{m,m'}||^2  <=  Lambda_{sqrt3/2}(m,m'),

where Lambda_h is the sum over ALL cusps of Gamma_0(4mm') of the Parseval integral over the strip {y >= h} in the
width coordinates of the cusp, i.e. the formula of cusp_lower_bound.py with every exponent multiplied by h:

   Lambda_h = sum_{delta | mm'} [ (1/4pi) sum_n |a_{A,B}(n)|^2 n^-1 (e^{-h pi delta n/mm'} + e^{-4 h pi delta n/mm'})
                                 + (1/pi) sum_n |a^odd_{A,B}(n)|^2 n^-1 e^{-h pi delta n/mm'} ].

h = 1: the disjoint regions above the cusps (lower bound, Theorem 1.5(i) and Section 6.6).
h = sqrt3/2: the covering by the translates of the standard fundamental domain (Part I, Lemma 2.8; upper bound).
The script prints Lambda_1, Lambda_{sqrt3/2}, the Rankin-Selberg estimate of ||G||^2 (all n <= X), the part of the
norm below height 1 (estimate minus Lambda_1) and its bound Lambda_{sqrt3/2} - Lambda_1, both divided by log mm'.
For (1,p) the bound is also split into the four Atkin-Lehner cusps and the two cusps of denominator 2.
Floating point; not certified.   Usage:  python norm_sandwich.py [X] [m,m' ...]
"""
import sys
import time
import numpy as np
from math import pi, gcd, log, sqrt
from sup_error_scan import r_table, sigma_over_n_odd, factor
from local_densities import siegel_all

args = [a for a in sys.argv[1:] if not a.startswith('--')]
X = int(args[0]) if args and ',' not in args[0] else 2 ** 21
pairs = [tuple(int(t) for t in a.split(',')) for a in args if ',' in a] or \
    [(1, 5), (3, 5), (3, 7), (5, 7), (1, 101), (13, 29), (15, 77), (1, 1009), (61, 67), (101, 103), (1, 15015),
     (1, 10007), (1, 100003), (1, 1000003)]
CUT = 14
H = sqrt(3) / 2

Xtab = max(X, max(CUT * m * mp for m, mp in pairs) + 1)
t0 = time.time()
r = r_table(Xtab)
f_sigma, v2 = sigma_over_n_odd(Xtab)
print(f"tables to {Xtab} in {time.time() - t0:.1f}s;  ||G||^2 estimated from n <= X = 2^{X.bit_length() - 1}\n")


def divisors(n):
    ds = [1]
    for p in factor(n):
        ds += [d * p for d in ds]
    return sorted(ds)


def norm_estimate(m, mp):
    L = 4 * m * mp
    T, M = siegel_all(m, mp, X, r[:X + 1], f_sigma[:X + 1], v2[:X + 1])
    n = np.arange(X + 1)
    a2 = (T - M) ** 2
    vol = (pi / 3) * L * np.prod([1 + 1 / p for p in factor(L)])
    return vol * np.mean(a2[1:] / n[1:]) / (16 * pi ** 2)


print(f"{'(m,m)':>12} {'mm':>8} {'||G||^2 est':>11} {'Lam_1':>8} {'Lam_h':>8} {'est-Lam_1':>9} {'Lam_h-Lam_1':>11} "
      f"{'(Lam_h-Lam_1)/log':>17} {'AL part':>8} {'odd part':>8} {'Lam_1 AL(1,n)':>13}")
for (m, mp) in pairs:
    assert gcd(m, mp) == 1 and m % 2 == 1 and mp % 2 == 1
    t1 = time.time()
    mm = m * mp
    est = norm_estimate(m, mp) if mm <= 20000 else float('nan')
    lam1 = lamh = 0.0
    al_diff = odd_diff = 0.0
    lam1_four = 0.0
    for delta in divisors(mm):
        dm, dmp = gcd(delta, m), gcd(delta, mp)
        A, B = dm * (mp // dmp), dmp * (m // dm)
        Nd = CUT * mm // delta + 1
        n = np.arange(Nd + 1, dtype=np.float64)
        T, M = siegel_all(A, B, Nd, r[:Nd + 1], f_sigma[:Nd + 1], v2[:Nd + 1])
        a2 = (T[1:] - M[1:]) ** 2 / n[1:]
        To, Mo = siegel_all(A, B, Nd, r[:Nd + 1], f_sigma[:Nd + 1], v2[:Nd + 1], odd=True)
        o2 = (To[1:] - Mo[1:]) ** 2 / n[1:]
        parts = {}
        for h in (1.0, H):
            w1 = np.exp(-h * pi * delta * n[1:] / mm)
            w4 = np.exp(-4 * h * pi * delta * n[1:] / mm)
            c_al = (np.sum(a2 * w1) + np.sum(a2 * w4)) / (4 * pi)
            c_odd = np.sum(o2 * w1) / pi
            parts[h] = (c_al, c_odd)
        lam1 += sum(parts[1.0])
        lamh += sum(parts[H])
        al_diff += parts[H][0] - parts[1.0][0]
        odd_diff += parts[H][1] - parts[1.0][1]
        if delta in (m, mp):
            lam1_four += parts[1.0][0]
    diff = lamh - lam1
    print(f"{str((m, mp)):>12} {mm:>8} {est:>11.4f} {lam1:>8.4f} {lamh:>8.4f} {est - lam1:>9.4f} {diff:>11.4f} "
          f"{diff / log(mm):>17.4f} {al_diff:>8.4f} {odd_diff:>8.4f} {lam1_four:>13.4f}   [{time.time() - t1:.0f}s]")
    sys.stdout.flush()

print("\nLam_1 <= ||G||^2 <= Lam_h (h = sqrt3/2).  est-Lam_1: the part of the norm below height 1 (floating point, "
      "needs the estimate to have converged: mm' <~ 10^4);\nLam_h-Lam_1: its proven upper bound, split into the "
      "Atkin-Lehner cusps (d_2 = 1, 4) and the cusps of denominator 2.  Lam_1 AL(1,n): the four cusps of Thm 1.5(i).")
