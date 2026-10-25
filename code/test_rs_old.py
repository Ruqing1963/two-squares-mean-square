r"""
Numerical test of the oldform case of the Rankin-Selberg formula (Part II, Theorem 6.7(ii)):

    <Theta_{1,p}, f(dz)> = (4/pi^2) chi(p) eta_p(f) L(1,f) L(1,f x chi_{-4}) / d        (f newform of level p or 2p, d | 4p/M)

for p prime.  The T_3-eigencomponent v of G_{1,p} inside the oldspace of f is v = c_1 f + c_2 f(2z) + c_4 f(4z)
(c_4 = 0 if the level of f is 2p); it is the orthogonal projection of G onto that space, so <Theta, v> = ||v||^2,
and the formula predicts ||v||^2 = (4/pi^2) chi(p) eta_p(f) L L (c_1 + c_2/2 + c_4/4).  The c_i and a_f(2) are read
off from a_v(1), a_v(2), a_v(4), a_v(8); ||v||^2 from Rankin-Selberg; L-values from the functional equation.
Usage: python test_rs_old.py p [p ...]
"""
import sys
import numpy as np
from math import pi, sqrt
from sup_error_scan import r_table, sigma_over_n_odd, factor
from local_densities import siegel_all

ps = [int(a) for a in sys.argv[1:]] or [11, 13, 29]
X = 2 ** 16


def T3(G, n_max):
    n = np.arange(n_max + 1)
    H = np.zeros(n_max + 1)
    H[1:] = G[3 * n[1:]]
    H[3::3] += 3 * G[n[3::3] // 3]
    return H


def L1(coef, Nl, nn):
    S = lambda t: np.sum(coef[1:] / nn[1:] * np.exp(-2 * pi * nn[1:] * t / sqrt(Nl)))
    best = None
    for w in (1, -1):
        vals = [S(t) + w * S(1 / t) for t in (0.6, 1.0, 1.5)]
        spread = max(vals) - min(vals)
        if best is None or spread < best[0]:
            best = (spread, w, float(np.mean(vals)))
    return best


for p in ps:
    N = 4 * p
    r = r_table(81 * X)
    fs, v2 = sigma_over_n_odd(81 * X)
    T, M = siegel_all(1, p, 81 * X, r, fs, v2)
    G = T - M
    seqs = [G[: 3 ** 10 + 1].copy()]
    for k in range(9):
        seqs.append(T3(seqs[-1], len(seqs[-1]) // 3))
    s = np.array([v[1] for v in seqs])
    lam = None
    for dd in range(1, 5):
        A = np.array([s[i:i + dd] for i in range(len(s) - dd)])
        coef = np.linalg.lstsq(A, s[dd:], rcond=None)[0]
        if np.abs(A @ coef - s[dd:]).max() < 1e-6 * np.abs(s).max():
            lam = np.sort(np.roots(np.r_[1, -coef[::-1]]).real)
            break
    print(f"p = {p}: T_3 eigenvalues on the span of G_{{1,{p}}}: {np.round(lam, 6)}")
    nn = np.arange(X + 1, dtype=float)
    chi = np.zeros(X + 1)
    chi[1::4] = 1
    chi[3::4] = -1
    for i, li in enumerate(lam):
        v = G[: (3 ** (len(lam) - 1)) * X + 1].copy()
        for j, lj in enumerate(lam):
            if j != i:
                v = T3(v, len(v) // 3) - lj * v[: len(v) // 3 + 1]
        v = v[: X + 1] / np.prod([li - lj for j, lj in enumerate(lam) if j != i])
        if abs(v[1]) < 1e-9:
            continue
        c1 = v[1]
        a_odd = v / c1
        # hypothesis B: f of level 2p (a_f(2) = +-1, c_4 = 0)
        a2B = v[4] / v[2] if abs(v[2]) > 1e-9 else None
        hypB = a2B is not None and abs(abs(a2B) - 1) < 1e-6 and abs(v[8] - a2B * v[4]) < 1e-6 * max(1, abs(v[8]))
        # hypothesis N: newform of level 4p (a(2) = 0)
        hypN = abs(v[2]) < 1e-9
        if hypN:
            print(f"   eigenvalue {li:.4f}: newform of level {N} (tested in test_rs_formula.py)")
            continue
        if hypB:
            a2, c2, c4, Mlev = a2B, v[2] - c1 * a2B, 0.0, 2 * p
            a_f2 = lambda j: a2 ** j
        else:
            # hypothesis A: f of level p; unknown a2: solve a_v(8) = c1 a_f(8) + c2 a_f(4) + c4 a_f(2)
            def resid(a2):
                c2 = v[2] - c1 * a2
                c4 = v[4] - c1 * (a2 ** 2 - 2) - c2 * a2
                return c1 * (a2 ** 3 - 4 * a2) + c2 * (a2 ** 2 - 2) + c4 * a2 - v[8], c2, c4
            grid = np.linspace(-2.9, 2.9, 5801)
            vals = np.array([resid(a)[0] for a in grid])
            roots = [grid[k] for k in range(len(grid) - 1) if vals[k] * vals[k + 1] <= 0]
            # refine each root by bisection and keep the one consistent with a_v(16)
            best = None
            for r0 in roots:
                lo, hi = r0, r0 + grid[1] - grid[0]
                for _ in range(60):
                    mid = (lo + hi) / 2
                    if resid(lo)[0] * resid(mid)[0] <= 0:
                        hi = mid
                    else:
                        lo = mid
                a2 = (lo + hi) / 2
                _, c2, c4 = resid(a2)
                a4, a8, a16 = a2 ** 2 - 2, a2 ** 3 - 4 * a2, a2 ** 4 - 6 * a2 ** 2 + 4
                err = abs(c1 * a16 + c2 * a8 + c4 * a4 - v[16])
                if best is None or err < best[0]:
                    best = (err, a2, c2, c4)
            err, a2, c2, c4 = best
            Mlev = p
            def a_f2(j, a2=a2):
                x, y = 1.0, a2
                for _ in range(j):
                    x, y = y, a2 * y - 2 * x
                return x
        # coefficients of f: a_f(2^j k) = a_f(2^j) a_f(k) for k odd, a_f(k) = a_v(k)/c1 for k odd
        af = np.zeros(X + 1)
        for k in range(1, X + 1, 2):
            af[k] = a_odd[k]
        j = 1
        while 2 ** j <= X:
            kk = np.arange(1, X // 2 ** j + 1, 2)
            af[kk * 2 ** j] = a_f2(j) * af[kk]
            j += 1
        vol = (pi / 3) * N * 1.5 * (1 + 1 / p)
        normv = np.sum(v[1:] ** 2 / nn[1:]) * vol / (16 * pi ** 2 * X)
        normv2 = np.sum(v[1:X // 2 + 1] ** 2 / nn[1:X // 2 + 1]) * vol / (16 * pi ** 2 * (X // 2))
        sp, w, Lf = L1(af, Mlev, nn)
        cand = [(Nl,) + L1(chi * af, Nl, nn) for Nl in (Mlev, 2 * Mlev, 4 * Mlev, 8 * Mlev, 16 * Mlev)]
        Nl, sp2, w2, Lfc = min(cand, key=lambda t: t[1])
        eta = -af[p]
        chip = 1 if p % 4 == 1 else -1
        pred = 4 / pi ** 2 * chip * eta * Lf * Lfc * (c1 + c2 / 2 + c4 / 4)
        print(f"   eigenvalue {li:.4f}: oldforms of f of level {Mlev}, a_f(2) = {a2:.6f}, a_f(3) = {af[3]:.4f}, "
              f"a_f({p}) = {af[p]:.4f}; v = {c1:.4f} f + {c2:.4f} f(2z) + {c4:.4f} f(4z)")
        print(f"      ||v||^2 (RS at X/2, X) = {normv2:.6f}, {normv:.6f};  L(1,f) = {Lf:.8f} (w={w:+d}, spread {sp:.1e}), "
              f"L(1,f x chi) = {Lfc:.8f} (level {Nl}, w={w2:+d}, spread {sp2:.1e})")
        print(f"      formula (4/pi^2) chi(p) eta_p L L (c1 + c2/2 + c4/4) = {pred:.6f};  ratio ||v||^2/formula = {normv / pred:.5f}")
