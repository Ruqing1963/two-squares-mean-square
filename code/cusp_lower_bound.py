r"""
Lower bound for ||G_{m,m'}||^2 from the Siegel regions {y >= 1} of the cusps of Gamma_0(4mm') (Part II, Theorem 1.4
and Section 6), for coprime odd SQUARE-FREE m, m'.

For every delta | mm' put delta_m = (delta, m), delta_m' = (delta, m'), A = delta_m m'/delta_m', B = delta_m' m/delta_m
(so AB = mm', (A, B) = 1).  The three cusps of Gamma_0(L), L = 4mm', with denominator d = d_2 delta, d_2 in {1, 2, 4},
have widths 4mm'/delta, mm'/delta, mm'/delta, and at them G_{m,m'} expands (up to a unimodular constant and the scaling
of the variable) as G_{A,B}, G^odd_{A,B}, G_{A,B}, where G^odd is the cusp part of the theta series of the odd vectors.
Integrating |G|^2 y^2 dmu over the regions {y >= 1} of all cusps (disjoint in Gamma_0(L)\H) gives

   ||G_{m,m'}||^2  >=  sum_{delta | mm'} [ (1/4pi) sum_n |a_{A,B}(n)|^2 n^-1 (e^{-pi delta n/mm'} + e^{-4 pi delta n/mm'})
                                           + (1/pi) sum_n |a^odd_{A,B}(n)|^2 n^-1 e^{-pi delta n/mm'} ],

a_{A,B}(n) = T(A,B;n) - M(A,B;n).  The two Atkin-Lehner terms with delta in {m, m'} (pair (1, mm')) are Theorem 1.4(i);
their asymptotic value is (1/pi)((log m)^2 + (log m')^2) (Theorem 1.4(ii)).  The script prints, for each pair, the
Rankin-Selberg estimate of ||G||^2 over ALL n <= X (local_densities.py makes the main term exact for every n), the
estimate over the admissible n of Part I (Table 1), the three parts of the lower bound, and the ratio.
Floating point; not certified.   Usage:  python cusp_lower_bound.py [X] [m,m' ...]  [--detail]
"""
import sys
import time
import numpy as np
from math import pi, gcd, log
from sup_error_scan import r_table, sigma_over_n_odd, factor
from local_densities import siegel_all

args = [a for a in sys.argv[1:] if not a.startswith('--')]
detail = '--detail' in sys.argv
X = int(args[0]) if args and ',' not in args[0] else 2 ** 21
pairs = [tuple(int(t) for t in a.split(',')) for a in args if ',' in a] or \
    [(1, 5), (3, 5), (3, 7), (5, 7), (1, 101), (13, 29), (15, 77), (1, 1009), (61, 67), (105, 143), (101, 103),
     (1, 15015), (1009, 1013)]
CUT = 12  # e^{-pi delta n/mm'} is kept for n <= CUT mm'/delta (e^{-12 pi} < 1e-16)

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


def norm_estimates(m, mp):
    """Rankin-Selberg estimate of ||G||^2 from n <= X: over all n, and over the admissible n of Part I."""
    L = 4 * m * mp
    T, M = siegel_all(m, mp, X, r[:X + 1], f_sigma[:X + 1], v2[:X + 1])
    n = np.arange(X + 1)
    a2 = (T - M) ** 2
    vol = (pi / 3) * L * np.prod([1 + 1 / p for p in factor(L)])
    all_n = vol * np.mean(a2[1:] / n[1:]) / (16 * pi ** 2)
    ok = (n >= 1) & (v2[:X + 1] <= 3)
    for p in factor(m * mp):
        ok &= n % p != 0
    adm = vol * np.mean(a2[ok] / n[ok]) / (16 * pi ** 2)
    return all_n, adm


rows = []
print(f"{'(m,m)':>12} {'mm':>8} {'||G||^2 all n':>13} {'adm. (Tab.1)':>12} {'LB_AL':>8} {'LB_odd':>8} {'LB':>8} "
      f"{'LB/||G||^2':>10} {'Thm(i)':>8} {'Thm(ii) asympt.':>15}")
for (m, mp) in pairs:
    assert gcd(m, mp) == 1 and m % 2 == 1 and mp % 2 == 1
    assert all(e == 1 for e in factor(m * mp).values()), "square-free only"
    t1 = time.time()
    G2_all, G2_adm = norm_estimates(m, mp)
    mm = m * mp
    lb_al = lb_odd = thm = 0.0
    det = []
    for delta in divisors(mm):
        dm, dmp = gcd(delta, m), gcd(delta, mp)
        A, B = dm * (mp // dmp), dmp * (m // dm)
        assert A * B == mm and gcd(A, B) == 1
        Nd = CUT * mm // delta + 1
        n = np.arange(Nd + 1, dtype=np.float64)
        T, M = siegel_all(A, B, Nd, r[:Nd + 1], f_sigma[:Nd + 1], v2[:Nd + 1])
        a2 = (T[1:] - M[1:]) ** 2 / n[1:]
        w1 = np.exp(-pi * delta * n[1:] / mm)
        w4 = np.exp(-4 * pi * delta * n[1:] / mm)
        c_al = (np.sum(a2 * w1) + np.sum(a2 * w4)) / (4 * pi)
        To, Mo = siegel_all(A, B, Nd, r[:Nd + 1], f_sigma[:Nd + 1], v2[:Nd + 1], odd=True)
        c_odd = np.sum((To[1:] - Mo[1:]) ** 2 / n[1:] * w1) / pi
        lb_al += c_al
        lb_odd += c_odd
        if delta in (m, mp):
            thm += c_al
        det.append((delta, A, B, c_al, c_odd))
    lb = lb_al + lb_odd
    asym = (log(m) ** 2 + log(mp) ** 2) / pi
    print(f"{str((m, mp)):>12} {mm:>8} {G2_all:>13.4f} {G2_adm:>12.4f} {lb_al:>8.4f} {lb_odd:>8.4f} {lb:>8.4f} "
          f"{lb / G2_all:>10.3f} {thm:>8.4f} {asym:>15.4f}   [{time.time() - t1:.0f}s]")
    sys.stdout.flush()
    rows.append((m, mp, G2_all, G2_adm, lb_al, lb_odd, lb, thm, asym, det))
    if detail:
        print(f"      {'delta':>8} {'(A,B)':>18} {'AL (d_2=1,4)':>13} {'odd (d_2=2)':>12}")
        for delta, A, B, c_al, c_odd in det:
            print(f"      {delta:>8} {str((A, B)):>18} {c_al:>13.4f} {c_odd:>12.4f}")

print("\nLB_AL: Atkin-Lehner cusps (d_2 = 1, 4), LB_odd: cusps of denominator 2 delta; Thm(i): the four cusps with the "
      "pair (1, mm') (delta = m, m');\nThm(ii) asympt.: (1/pi)((log m)^2 + (log m')^2), the leading term of Thm(i).  "
      "The Siegel regions cover the fraction 3/pi = 0.955 of the area of Gamma_0(L)\\H.")
