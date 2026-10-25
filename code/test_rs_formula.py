r"""
Numerical test of the Rankin-Selberg formula (Part II, Theorem 6.7)

    <Theta_{1,p}, phi> = (4/pi^2) chi(p) eta_p(phi) L(1,phi) L(1, phi x chi_{-4})

for phi a newform of weight 2 and level 4p.  The cuspidal part G_{1,p} = T - M (exact main term) is decomposed with
the Hecke operator T_3 (Prony on the sequence (T_3^k G)(1) gives the distinct eigenvalues on the span of G); the
eigencomponent with a(2) = 0 and a(9) = a(3)^2 - 3 is a newform phi of level 4p; its norm comes from Rankin-Selberg
(sum_{k<=y} a(k)^2/k ~ 16 pi^2 ||phi||^2 y / vol), its central values from the smoothed sums with the functional
equation (the root numbers and the level of the twist are detected by t-independence).  eta_p = -a_phi(p).
Usage: python test_rs_formula.py p [p ...]
"""
import sys
import numpy as np
from math import pi, sqrt
from sup_error_scan import r_table, sigma_over_n_odd
from local_densities import siegel_all

ps = [int(a) for a in sys.argv[1:]] or [5, 13, 19]
X = 2 ** 16


def T3(G, n_max):
    n = np.arange(n_max + 1)
    H = np.zeros(n_max + 1)
    H[1:] = G[3 * n[1:]]
    H[3::3] += 3 * G[n[3::3] // 3]
    return H


for p in ps:
    N = 4 * p
    r = r_table(81 * X)
    fs, v2 = sigma_over_n_odd(81 * X)
    T, M = siegel_all(1, p, 81 * X, r, fs, v2)
    G = T - M
    # Prony on s_k = (T_3^k G)(1), k = 0..9 (G needed up to 3^9 < 81 X)
    seqs = [G[: 3 ** 10 + 1].copy()]
    for k in range(9):
        seqs.append(T3(seqs[-1], len(seqs[-1]) // 3))
    s = np.array([v[1] for v in seqs])
    # find the smallest order d with a linear recurrence; eigenvalues = roots
    lam = None
    for d in range(1, 5):
        A = np.array([s[i:i + d] for i in range(len(s) - d)])
        b = s[d:]
        coef, res, rk, sv = np.linalg.lstsq(A, b, rcond=None)
        resid = np.abs(A @ coef - b).max()
        if resid < 1e-6 * np.abs(s).max():
            roots = np.roots(np.r_[1, -coef[::-1]])
            lam = np.sort(roots.real)
            break
    print(f"p = {p}, level {N}: T_3 eigenvalues on the span of G_{{1,{p}}}: {np.round(lam, 6)}")
    # eigencomponents v_i = prod_{j != i} (T_3 - lam_j) G / prod (lam_i - lam_j), on n <= X
    for i, li in enumerate(lam):
        v = G[: (3 ** (len(lam) - 1)) * X + 1].copy()
        for j, lj in enumerate(lam):
            if j != i:
                v = T3(v, len(v) // 3) - lj * v[: len(v) // 3 + 1]
        v = v[: X + 1] / np.prod([li - lj for j, lj in enumerate(lam) if j != i])
        if abs(v[1]) < 1e-9:
            continue
        a = v / v[1]
        new = abs(a[2]) < 1e-6 and abs(a[9] - (a[3] ** 2 - 3)) < 1e-6 and abs(a[25] - (a[5] ** 2 - 5)) < 1e-6
        print(f"   eigenvalue {li:.4f}: a(2) = {a[2]:.4f}, a(3) = {a[3]:.4f}, a(5) = {a[5]:.4f}, a({p}) = {a[p]:.4f},"
              f" a(9)-(a(3)^2-3) = {a[9] - (a[3] ** 2 - 3):.1e}  -> {'NEWFORM of level %d' % N if new else 'old'}")
        if not new:
            continue
        c = v[1]
        nn = np.arange(X + 1, dtype=float)
        vol = (pi / 3) * N * 1.5 * (1 + 1 / p)
        norms = [np.sum(a[1:y + 1] ** 2 / nn[1:y + 1]) * vol / (16 * pi ** 2 * y) for y in (X // 4, X // 2, X)]
        norm = norms[-1]
        print(f"      ||phi||^2 (RS at y = X/4, X/2, X): {norms[0]:.6f} {norms[1]:.6f} {norms[2]:.6f};"
              f"  c = {c:.6f};  <Theta,phi> = c||phi||^2 = {c * norm:.6f}")
        chi = np.zeros(X + 1)
        chi[1::4] = 1
        chi[3::4] = -1

        def L1(coef, Nl):
            S = lambda t: np.sum(coef[1:] / nn[1:] * np.exp(-2 * pi * nn[1:] * t / sqrt(Nl)))
            best = None
            for w in (1, -1):
                vals = [S(t) + w * S(1 / t) for t in (0.6, 1.0, 1.5)]
                spread = max(vals) - min(vals)
                if best is None or spread < best[0]:
                    best = (spread, w, np.mean(vals))
            return best

        sp, w, Lf = L1(a, N)
        print(f"      L(1,phi) = {Lf:.8f} (root number {w:+d}, spread {sp:.1e})")
        cand = [(Nl,) + L1(chi * a, Nl) for Nl in (N, 4 * N, 16 * N)]
        Nl, sp2, w2, Lfc = min(cand, key=lambda t: t[1])
        print(f"      L(1,phi x chi) = {Lfc:.8f} (twist level {Nl}, root number {w2:+d}, spread {sp2:.1e})")
        eta = -a[p]
        chip = 1 if p % 4 == 1 else -1
        pred = 4 / pi ** 2 * chip * eta * Lf * Lfc
        print(f"      formula (4/pi^2) chi(p) eta_p L L = {pred:.6f}   (chi(p) = {chip:+d}, eta_p = {eta:+.0f});"
              f"   ratio computed/formula = {c * norm / pred:.5f}")
