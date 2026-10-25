r"""
Spectral support of G_{1,q}, q prime (Part II, Section 6.8): which Hecke eigencomponents does the cuspidal part of
theta^2(z) theta^2(qz) have?  The T_3-eigencomponents of G_{1,q} are found by Prony's method on (T_3^k G)(1) and
classified by a(2), a(4), a(8) (normalised by a(1)):
   newform of level 4q:         a(2) = a(4) = 0;
   oldforms of a newform f of level q or 2q:  a(2) = a_f(2) + c_2/c_1 etc.
Theorem 6.7 and the Atkin-Lehner symmetries say that a newform phi of level 4q occurs iff eta_q(phi) = chi(q),
eta_4(phi) = -chi(q) and L(1,phi) L(1,phi x chi) != 0; since eta_4 = -1 for every newform of level 4q, no newform
occurs when q = 3 (mod 4).  The script prints, for each q, the components with a(2), a(3), a(5), a(q), the share of
||G||^2 they carry (Rankin-Selberg over n <= X_norm), and whether any newform of level 4q occurs.
Usage: python newform_support.py q [q ...]
"""
import sys
import numpy as np
from math import pi
from sup_error_scan import r_table, sigma_over_n_odd
from local_densities import siegel_all

qs = [int(a) for a in sys.argv[1:]] or [11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 5, 13, 17, 29, 37, 41, 53, 61]
KMAX = 14            # Prony on s_0..s_KMAX, recurrence order up to KMAX//2
X = 2 ** 12          # components computed on n <= X (needs G up to 3^(dim-1) X)
Xn = 2 ** 16         # Rankin-Selberg window for the norm shares (only for dim <= 5)


def T3(G, n_max):
    n = np.arange(n_max + 1)
    H = np.zeros(n_max + 1)
    H[1:] = G[3 * n[1:]]
    H[3::3] += 3 * G[n[3::3] // 3]
    return H


Gmax = 3 ** KMAX + 1
r = r_table(Gmax)
fs, v2 = sigma_over_n_odd(Gmax)
for q in qs:
    N = 4 * q
    T, M = siegel_all(1, q, Gmax, r, fs, v2)
    G = T - M
    seqs = [G[: 3 ** KMAX + 1].copy()]
    for k in range(KMAX):
        seqs.append(T3(seqs[-1], len(seqs[-1]) // 3))
    s = np.array([v[1] for v in seqs])
    lam = None
    for dd in range(1, KMAX // 2 + 1):
        A = np.array([s[i:i + dd] for i in range(len(s) - dd)])
        coef = np.linalg.lstsq(A, s[dd:], rcond=None)[0]
        if np.abs(A @ coef - s[dd:]).max() < 1e-7 * np.abs(s).max():
            lam = np.sort(np.roots(np.r_[1, -coef[::-1]]).real)
            break
    if lam is None:
        print(f"q = {q}: more than {KMAX // 2} distinct T_3-eigenvalues; skipped")
        continue
    dim = len(lam)
    Xc = X if 3 ** (dim - 1) * X < Gmax else (Gmax - 1) // 3 ** (dim - 1)
    chi_q = 1 if q % 4 == 1 else -1
    print(f"q = {q} (chi(q) = {chi_q:+d}), level {N}: {dim} T_3-eigencomponents, eigenvalues {np.round(lam, 4)}")
    any_new = False
    unresolved = False
    vol = (pi / 3) * N * 1.5 * (1 + 1 / q)
    nn = np.arange(Xn + 1, dtype=float)
    Gn = None
    if dim <= 5:
        Yn = 3 ** (dim - 1) * Xn + 1
        rn = r_table(Yn)
        fsn, v2n = sigma_over_n_odd(Yn)
        Tn, Mn = siegel_all(1, q, Yn, rn, fsn, v2n)
        Gn = Tn - Mn
        normG = np.sum(Gn[1:Xn + 1] ** 2 / nn[1:]) * vol / (16 * pi ** 2 * Xn)
    def T5(v, n_max):
        n = np.arange(n_max + 1)
        H = np.zeros(n_max + 1)
        H[1:] = v[5 * n[1:]]
        H[5::5] += 5 * v[n[5::5] // 5]
        return H

    def classify(v, label):
        """v: a T_3-eigencomponent (or a T_5-subcomponent) on n <= len(v)-1; prints its type."""
        nonlocal_any = False
        a = v / v[1]
        new = abs(a[2]) < 1e-6 and abs(a[4]) < 1e-6
        kind = f"NEWFORM of level {N}" if new else "old (levels q, 2q)"
        share = ""
        if Gn is not None and len(v) == Xn + 1:
            share = f", share of ||G||^2 = {np.sum(v[1:] ** 2 / nn[1:]) * vol / (16 * pi ** 2 * Xn) / normG:.3f}"
        print(f"   {label}: a(2) = {a[2]:+.4f}, a(4) = {a[4]:+.4f}, a(5) = {a[5]:+.4f}, a(q) = {a[q]:+.4f}"
              f"  -> {kind}{share}")
        return new

    for i, li in enumerate(lam):
        src = Gn if Gn is not None else G
        Xi = Xn if Gn is not None else Xc
        v = src[: 3 ** (dim - 1) * Xi + 1].copy()
        for j, lj in enumerate(lam):
            if j != i:
                v = T3(v, len(v) // 3) - lj * v[: len(v) // 3 + 1]
        v = v[: Xi + 1] / np.prod([li - lj for j, lj in enumerate(lam) if j != i])
        if abs(v[1]) < 1e-9 * max(1, np.abs(v).max()):
            continue
        a = v / v[1]
        # a single eigenform (times oldform shifts) has a(15) = a(3) a(5); a mixture of two systems with the same
        # a(3) does not: split it further with T_5
        if q != 5 and abs(a[15] - a[3] * a[5]) > 1e-6 and len(v) > 5 ** 5:
            s5 = [v[1]]
            w = v.copy()
            for k in range(4):
                w = T5(w, (len(w) - 1) // 5)
                s5.append(w[1])
            s5 = np.array(s5)
            lam5 = None
            for dd in (1, 2):
                A5 = np.array([s5[t:t + dd] for t in range(len(s5) - dd)])
                coef = np.linalg.lstsq(A5, s5[dd:], rcond=None)[0]
                if np.abs(A5 @ coef - s5[dd:]).max() < 1e-7 * np.abs(s5).max():
                    lam5 = np.sort(np.roots(np.r_[1, -coef[::-1]]).real)
                    break
            if lam5 is not None and len(lam5) == 2:
                for ii, mu in enumerate(lam5):
                    u = v.copy()
                    nu = (len(u) - 1) // 5
                    u = (T5(u, nu) - lam5[1 - ii] * u[: nu + 1]) / (mu - lam5[1 - ii])
                    if abs(u[1]) > 1e-9:
                        any_new |= classify(u, f"a(3) = {li:+.4f}, a(5) = {mu:+.4f} (split by T_5)")
                continue
            print(f"   a(3) = {li:+.4f}: mixture or inaccurate component (a(15) != a(3)a(5)), not separated by T_5")
            unresolved = True
            continue
        any_new |= classify(v, f"a(3) = {li:+.4f}")
    verdict = 'YES' if any_new else ('INCONCLUSIVE (unresolved components)' if unresolved else 'none')
    print(f"   newform component of level {N}: {verdict}"
          f"   (Theorem: none when chi(q) = -1; when chi(q) = +1 the newforms with eta_q = +1 and L L != 0 occur)")
    sys.stdout.flush()
