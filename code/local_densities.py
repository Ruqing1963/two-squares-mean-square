r"""
Exact Siegel main term of theta^2(Az) theta^2(Bz) for ALL n (not only the admissible ones), for coprime odd
SQUARE-FREE A, B, and the same for the theta series of the odd vectors (the expansion at the cusps of
denominator 2 delta).  Used by cusp_lower_bound.py (Part II, Section 6).

    M(A,B;n) = pi^2 n / (AB) * beta_2(n) * prod_{p odd} beta_p(n),

    beta_p(n) = (1 - p^-2) sigma_{-1}(p^{v_p(n)})                         (p nmid 2AB)
    beta_p(n) = (1 + 1/p) p^{-k}                 (p | AB, p = 3 mod 4, k = v_p(n))
    beta_p(n) = 2 - (1 + 1/p) p^{-k}             (p | AB, p = 1 mod 4, k = v_p(n))
    beta_2(n): counts modulo 2^nu for v_2(n) <= 3 (sup_error_scan.delta2_table), and for v_2(n) >= 2 the
               recursion  beta_2(n) = beta_2(n/4)/4 + [A = B mod 4][v_2(n) = 2]/2
                                     + [A != B mod 4] (1 + [v_2(n) >= 3]/2),
               which is checked against the table at v_2(n) = 2, 3 and against direct counts up to v_2(n) = 8.

The recursion for p | AB: solutions of A(x^2+y^2) + B(z^2+w^2) = p n' with p | A are either non-singular
((z, w) != 0 mod p, which requires p = 1 mod 4 and contributes 2(1 - 1/p)) or have z = w = 0 mod p, which
reduces to the form Q_{A/p, pB} and n' with the factor 1/p.  The 2-adic recursion sorts the solutions by the
parities of (x, y) and (z, w): all even (factor 1/4 and n/4), all odd (density 1/2 exactly when
n = 2(A + B) mod 8), one odd in each pair (density 1 exactly when n = A + B mod 4).

Odd vectors: r_odd(alpha) = r(alpha) [alpha = 2 mod 4], T_odd(A,B;n) = sum_{A alpha + B beta = n} r_odd r_odd,
M_odd(A,B;n) = pi^2 n/(AB) * (1/2)[n = 2(A+B) mod 8] * prod_{p odd} beta_p(n)   (Siegel's theorem for the
coset (1+2Z)^4; for (A,B) = (1,1) this is Jacobi's 16 sigma(n/4)).
"""
import numpy as np
from math import pi
from scipy.signal import fftconvolve
from sup_error_scan import r_table, sigma_over_n_odd, delta2_table, factor, chi4, counts_mod


def vp_table(X, p):
    """v_p(n) for 0 <= n <= X (v_p(0) := 0)."""
    v = np.zeros(X + 1, dtype=np.int64)
    q = p
    while q <= X:
        v[q::q] += 1
        q *= p
    return v


def beta2_all(A, B, X, v2, check_to=8):
    """beta_2(Q_{A,B}; n) for 0 <= n <= X."""
    n = np.arange(X + 1)
    tab = delta2_table(A, B)                      # exact for v_2(n) <= 3, indexed by n mod 512
    b = tab[n % 512].copy()
    same = (A - B) % 4 == 0
    vmax = int(v2.max())
    # recursion for v >= 2 (apply from v = 2 upwards so that b[n/4] is final when used)
    for v in range(2, vmax + 1):
        idx = np.nonzero(v2 == v)[0]
        if len(idx) == 0:
            continue
        val = b[idx // 4] / 4                      # b[n/4] is final: table for v <= 3, recursion before
        if same:
            val = val + (0.5 if v == 2 else 0.0)
        else:
            val = val + 1.0 + (0.5 if v >= 3 else 0.0)
        if v <= 3:
            assert np.allclose(val, b[idx], rtol=1e-12, atol=1e-12), (A, B, v)
        else:
            b[idx] = val
    # direct check of the recursion against counts modulo 2^(v+6) for v <= check_to
    for v in range(4, min(check_to, vmax) + 1):
        nu = v + 6
        cnt = counts_mod(A, B, 2 ** nu)
        idx = np.nonzero(v2 == v)[0][:64]
        direct = cnt[idx % 2 ** nu] / 8.0 ** nu
        assert np.allclose(direct, b[idx], rtol=1e-9, atol=1e-12), (A, B, v, direct[:3], b[idx][:3])
    return b


def siegel_all(A, B, X, r=None, f_sigma=None, v2=None, odd=False):
    """T(A,B;n) and M(A,B;n) for 0 <= n <= X (exact for every n), or the odd-vector versions."""
    if r is None:
        r = r_table(X)
    if f_sigma is None:
        f_sigma, v2 = sigma_over_n_odd(X)
    n = np.arange(X + 1, dtype=np.float64)
    bad = [p for p in factor(A * B)]
    for p in bad:
        assert (A * B) % (p * p) != 0, "square-free A, B only"
    # product over p nmid 2AB of (1 - p^-2) sigma_{-1}(p^{v_p}) = (8/pi^2) prod_{p|AB} (1-p^-2)^{-1} * sigma(n_odd)/n_odd
    #   * prod_{p|AB} sigma_{-1}(p^{v_p(n)})^{-1}
    good = (8 / pi ** 2) / np.prod([1 - p ** -2.0 for p in bad]) if bad else 8 / pi ** 2
    loc = f_sigma.copy()
    for p in bad:
        k = vp_table(X, p)
        sig = (1 - p ** -(k + 1.0)) / (1 - 1.0 / p)          # sigma_{-1}(p^k)
        if p % 4 == 3:
            bp = (1 + 1.0 / p) * p ** (-k.astype(float))
        else:
            bp = 2 - (1 + 1.0 / p) * p ** (-k.astype(float))
        loc *= bp / sig
    if odd:
        b2 = 0.5 * ((np.arange(X + 1) - 2 * (A + B)) % 8 == 0)
        rr = r * ((np.arange(X + 1) % 4) == 2)
        rr[0] = 0.0
    else:
        b2 = beta2_all(A, B, X, v2)
        rr = r
    M = pi ** 2 * n / (A * B) * good * loc * b2
    a = np.zeros(X + 1)
    a[::A] = rr[: X // A + 1]
    b = np.zeros(X + 1)
    b[::B] = rr[: X // B + 1]
    T = np.round(fftconvolve(a, b)[: X + 1])
    return T, M


if __name__ == "__main__":
    import time
    X = 2 ** 17
    t0 = time.time()
    r = r_table(X)
    f_sigma, v2 = sigma_over_n_odd(X)
    print(f"tables to {X} in {time.time() - t0:.1f}s")
    # levels 4 and 12 have no cusp forms: T = M exactly for every n, also for the odd vectors
    for (A, B) in ((1, 1), (1, 3)):
        for odd in (False, True):
            T, M = siegel_all(A, B, X, r, f_sigma, v2, odd=odd)
            err = np.abs(T[1:] - M[1:]).max()
            print(f"(A,B) = ({A},{B}), odd = {odd}: max |T - M| over n <= {X} is {err:.2e}")
            assert err < 1e-6
    # level 20: T - M = (8/3) a_f(n) for every n, f_20 = eta(2z)^2 eta(10z)^2 = q P(q^2)^2 P(q^10)^2
    def euler(Xmax, step):
        s = np.zeros(Xmax + 1)
        k = 0
        while True:
            e1, e2 = k * (3 * k - 1) // 2 * step, k * (3 * k + 1) // 2 * step
            if e1 > Xmax and e2 > Xmax:
                break
            if e1 <= Xmax:
                s[e1] += (-1) ** k
            if e2 <= Xmax and k > 0:
                s[e2] += (-1) ** k
            k += 1
        return s

    def sq(s):
        return np.round(fftconvolve(s, s)[: len(s)])

    F = np.round(fftconvolve(sq(euler(X, 2)), sq(euler(X, 10)))[: X + 1])
    af = np.zeros(X + 1)
    af[1:] = F[:X]
    T, M = siegel_all(1, 5, X, r, f_sigma, v2)
    err = np.abs(T[1:] - M[1:] - 8 / 3 * af[1:]).max()
    print(f"(A,B) = (1,5): max |T - M - (8/3) a_f(n)| over ALL n <= {X} is {err:.2e}")
    assert err < 1e-6
    T, M = siegel_all(1, 5, X, r, f_sigma, v2, odd=True)
    # the odd-vector cusp part must be a cusp form: its mean square grows like Y^2, not Y^3
    a = T - M
    nn = np.arange(X + 1)
    for y in (X // 8, X // 4, X // 2, X):
        print(f"   odd vectors, (1,5): sum_(n<=y) a_odd(n)^2 / y^2 at y = {y}: {np.sum(a[1:y + 1] ** 2) / y ** 2:.4f}")
    print("all checks passed")
