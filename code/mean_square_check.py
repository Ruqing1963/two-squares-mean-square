r"""
Mean square of the cuspidal error a_G(N) = T(m,m';N) - M(m,m';N) over N <= Y, and the Petersson norm of G
estimated from it by the Rankin-Selberg asymptotic

    sum_{N <= Y} |a_G(N)|^2 / N  ~  (16 pi^2 ||G||^2 / vol(Gamma_0(L)\H)) * Y,      L = 4 m m',

vol = (pi/3) L prod_{p | L} (1 + 1/p) (measure dx dy / y^2).  Theorem A of the supplement states, rigorously,

    sum_{N <= Y} |a_G(N)|^2  <=  e^{8 pi} Y (3 + 17 Y/L) ||G||^2,      ||G||^2 << (mm')^{1/2 + eps},

so that the root mean square of (T - M)/sqrt(N) over N <= Y is << (mm')^{-1/4 + eps} once Y >= mm'.

The sums below run over the admissible N ((N, mm') = 1, v_2(N) <= 3) only, where the main term of
sup_error_scan.py is exact; they are therefore lower bounds for the full sums, and the norm estimate assumes the
mean square per N to be the same on the remaining N.  Floating point, not certified.
Usage:  python mean_square_check.py [X] [m,m' ...]
"""
import sys
import time
import numpy as np
from sup_error_scan import r_table, sigma_over_n_odd, delta2_table, T_all, M_all, factor

args = sys.argv[1:]
X = int(args[0]) if args and ',' not in args[0] else 2 ** 21
pairs = [tuple(int(t) for t in a.split(',')) for a in args if ',' in a] or \
    [(1, 5), (9, 1), (3, 5), (1, 101), (13, 29), (1, 1009), (61, 67), (101, 103), (1, 15015), (1009, 1013)]

t0 = time.time()
r = r_table(X)
f_sigma, v2 = sigma_over_n_odd(X)
print(f"X = 2^{X.bit_length() - 1}; tables in {time.time() - t0:.1f}s\n")
print(f"{'(m,m)':>13} {'mm':>8} {'L':>8} {'#N':>8} {'S2(Y)/S2(Y/2)':>14} {'rms/(mm)^-1/4':>14} "
      f"{'||G||^2 est':>12} {'||G||^2/(mm)^1/2':>16} {'||G||^2/(mm)^1/4':>16}")
n_all = np.arange(X + 1)
for (m, mp) in pairs:
    L = 4 * m * mp
    dens2 = delta2_table(m, mp)
    T = T_all(r, m, mp, X)
    M = M_all(m, mp, X, f_sigma, v2, dens2)
    ok = (n_all >= 1) & (v2 <= 3)
    for p in factor(m * mp):
        ok &= n_all % p != 0
    a2 = (T - M) ** 2
    S2 = a2[ok].sum()
    S2_half = a2[ok & (n_all <= X // 2)].sum()
    mean_a2_over_N = np.mean(a2[ok] / n_all[ok])
    rms = np.sqrt(mean_a2_over_N)
    vol = (np.pi / 3) * L * np.prod([1 + 1 / p for p in factor(L)])
    G2 = vol * mean_a2_over_N / (16 * np.pi ** 2)
    print(f"{str((m, mp)):>13} {m * mp:>8} {L:>8} {ok.sum():>8} {S2 / S2_half:>14.3f} {rms * (m * mp) ** 0.25:>14.3f} "
          f"{G2:>12.4f} {G2 / (m * mp) ** 0.5:>16.4f} {G2 / (m * mp) ** 0.25:>16.4f}")
print("\nS2(Y)/S2(Y/2) ~ 4 is the Y^2 growth of Theorem A; ||G||^2/(mm')^{1/2} bounded is the sharpness of Lemma 2.1.")
