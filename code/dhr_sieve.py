r"""
Diamond-Halberstam-Richert sieve functions F_kappa, f_kappa of dimension kappa (kappa = 2 for n(N-n); kappa = 1
reproduces the linear sieve as a check), and the weighted two-dimensional sieve for the product n(N-n) under the
weight r(n) r(N-n) (Part II, Section 7 / Remark 3.2).

Sieve functions (DHR 1988; Diamond-Halberstam, A Higher-Dimensional Sieve Method, 2008):
    sigma(s) = A s^k (0 < s <= 2),  (s^-k sigma(s))' = -k s^(-k-1) sigma(s-2)  (s > 2),  A = (2e^gamma)^-k / Gamma(k+1);
    F(s) = 1/sigma(s) (0 < s <= alpha),  f(s) = 0 (0 < s <= beta),
    (s^k F)' = k s^(k-1) f(s-1)  (s > alpha),   (s^k f)' = k s^(k-1) F(s-1)  (s > beta),
and alpha = alpha_k, beta = beta_k are the values for which F, f -> 1.  We determine them through the adjoint
functions p, q:  s^k p'(s) = -k (s+1)^(k-1) p(s+1),  s^k q'(s) = k (s+1)^(k-1) q(s+1),  p ~ s^-k, q ~ s^k, for which
    Phi_P(s) = s^k P(s) p(s) + k int_{s-1}^{s} (u+1)^(k-1) P(u) p(u+1) du    (P = F + f)
    Phi_Q(s) = s^k Q(s) q(s) - k int_{s-1}^{s} (u+1)^(k-1) Q(u) q(u+1) du    (Q = F - f)
are constant in s >= max(alpha, beta) + 1; F, f -> 1 iff Phi_P = 2 and Phi_Q = 0.  Check: kappa = 1 gives
alpha = beta = 2 and the linear sieve functions.

Weighted sieve for n(N-n) (dimension 2, level D = N^theta, elements <= N^2): sift by p < N^a, Richert weights
1 - lam sum_{N^a <= p < N^b, p | n(N-n)} (1 - log p / (b log N)); positivity needs
    L2 - lam R2 > 0,   L2 = f_2(theta/a),   R2 = 2 int_a^b (1 - t/b) F_2((theta - t)/a) dt/t,
and then Omega(n(N-n)) < 1/lam + 2/b and max(Omega(n), Omega(N-n)) < 1/lam + 1/b.  The script minimises
1/lam + 2/b = R2/L2 + 2/b and R2/L2 + 1/b over (a, b) and compares with the vector sieve of Part I.
Floating point; not certified.   Usage: python dhr_sieve.py [theta ...]
"""
import sys
import numpy as np
from math import gamma as Gamma, exp, log, comb
from scipy.optimize import fsolve

EG = exp(np.euler_gamma)
H = 1.0 / 2000
S_TOP = 48.0
grid = np.arange(0.0, S_TOP + H / 2, H)
N1 = int(round(1 / H))
N2 = 2 * N1


def sigma_table(k):
    """g(s) = s^-k sigma(s) on the grid (g = A on (0,2])."""
    A = (2 * EG) ** (-k) / Gamma(k + 1)
    g = np.full_like(grid, A)
    s = grid
    i2 = N2
    # g'(s) = -k s^(-k-1) (s-2)^k g(s-2)
    for i in range(i2 + 1, len(s)):
        d0 = -k * s[i - 1] ** (-k - 1) * s[i - 1 - N2] ** k * g[i - 1 - N2]
        d1 = -k * s[i] ** (-k - 1) * s[i - N2] ** k * g[i - N2]
        g[i] = g[i - 1] + H * (d0 + d1) / 2
    return g, A


def gbinom(x, i):
    out = 1.0
    for j in range(i):
        out *= (x - j) / (j + 1)
    return out


def adjoint_tables(k, nterms=14):
    """p and q on the grid (p ~ s^-k, q ~ s^k), by backward integration from an asymptotic series."""
    a = [1.0]
    for n in range(1, nterms):
        a.append(k * sum(a[j] * gbinom(-j - 1, n - j) for j in range(n)) / n)
    b = [1.0]
    for n in range(1, nterms):
        b.append(-k * sum(b[j] * gbinom(2 * k - 1 - j, n - j) for j in range(n)) / n)
    s = grid
    p = np.zeros_like(s)
    q = np.zeros_like(s)
    top = len(s) - 1
    start = top - N1
    ss = s[start:]
    p[start:] = ss ** (-k) * sum(a[j] * ss ** (-j) for j in range(nterms))
    q[start:] = ss ** k * sum(b[j] * ss ** (-j) for j in range(nterms))
    for i in range(start - 1, N1 - 1, -1):   # down to s = 1/H... stop at s = 1
        dp1 = -k * (s[i + 1] + 1) ** (k - 1) * p[i + 1 + N1] / s[i + 1] ** k
        dp0 = -k * (s[i] + 1) ** (k - 1) * p[i + N1] / s[i] ** k
        p[i] = p[i + 1] - H * (dp0 + dp1) / 2
        dq1 = k * (s[i + 1] + 1) ** (k - 1) * q[i + 1 + N1] / s[i + 1] ** k
        dq0 = k * (s[i] + 1) ** (k - 1) * q[i + N1] / s[i] ** k
        q[i] = q[i + 1] - H * (dq0 + dq1) / 2
    return p, q


def sieve_functions(k, alpha, beta, g):
    """F, f on the grid for given alpha, beta (F = 1/sigma up to alpha, f = 0 up to beta)."""
    s = grid
    F = np.zeros_like(s)
    f = np.zeros_like(s)
    with np.errstate(divide='ignore', invalid='ignore'):
        F[1:] = 1.0 / (s[1:] ** k * g[1:])
    ia = int(np.searchsorted(s, alpha, side='right'))   # first index with s > alpha
    ib = int(np.searchsorted(s, beta, side='right'))
    F[ia:] = 0.0
    u = s ** k * F                                      # u = s^k F, v = s^k f
    v = np.zeros_like(s)

    def interp(arr, x):
        return float(np.interp(x, s, arr))

    # partial first steps from alpha and beta
    ua = 1.0 / interp(g, alpha)                         # u(alpha) = alpha^k F(alpha) = alpha^k / sigma(alpha) = 1/g(alpha)
    u_prev_s, u_prev = alpha, ua
    v_prev_s, v_prev = beta, 0.0
    for i in range(min(ia, ib), len(s)):
        si = s[i]
        if i >= ia:
            if i == ia:
                d0 = k * alpha ** (k - 1) * interp(f, alpha - 1)
                d1 = k * si ** (k - 1) * f[i - N1]
                u[i] = ua + (si - alpha) * (d0 + d1) / 2
            else:
                d0 = k * s[i - 1] ** (k - 1) * f[i - 1 - N1]
                d1 = k * si ** (k - 1) * f[i - N1]
                u[i] = u[i - 1] + H * (d0 + d1) / 2
            F[i] = u[i] / si ** k
        if i >= ib:
            if i == ib:
                d0 = k * beta ** (k - 1) * interp(F, beta - 1)
                d1 = k * si ** (k - 1) * F[i - N1]
                v[i] = (si - beta) * (d0 + d1) / 2
            else:
                d0 = k * s[i - 1] ** (k - 1) * F[i - 1 - N1]
                d1 = k * si ** (k - 1) * F[i - N1]
                v[i] = v[i - 1] + H * (d0 + d1) / 2
            f[i] = v[i] / si ** k
    return F, f


def phis(k, F, f, p, q, s0):
    s = grid
    i0 = int(round(s0 / H))
    P = F + f
    Q = F - f
    sl = slice(i0 - N1, i0 + 1)
    w = (s[sl] + 1) ** (k - 1)
    phiP = s0 ** k * P[i0] * p[i0] + k * np.trapz(w * P[sl] * p[i0 - N1 + N1: i0 + 1 + N1], dx=H)
    phiQ = s0 ** k * Q[i0] * q[i0] - k * np.trapz(w * Q[sl] * q[i0 - N1 + N1: i0 + 1 + N1], dx=H)
    return phiP, phiQ


def solve_alpha_beta(k, guess, g, p, q, verbose=True):
    def resid(x):
        al, be = x
        F, f = sieve_functions(k, al, be, g)
        s0 = round((max(al, be) + 1.5) / H) * H
        P, Q = phis(k, F, f, p, q, s0)
        return [P - 2.0, Q]
    sol, info, ier, msg = fsolve(resid, guess, full_output=True, xtol=1e-10)
    if verbose:
        print(f"  kappa = {k}: alpha = {sol[0]:.6f}, beta = {sol[1]:.6f}  (residual {np.abs(info['fvec']).max():.1e}, {msg.strip()})")
    return sol


class DHR:
    def __init__(self, k, guess=None, verbose=True):
        self.k = k
        self.g, self.A = sigma_table(k)
        self.p, self.q = adjoint_tables(k)
        if guess is None:
            guess = (2.0, 2.0) if k == 1 else (5.35, 4.28)   # kappa = 2: the root located by a scan of (Phi_P - 2, Phi_Q)
        self.alpha, self.beta = solve_alpha_beta(k, guess, self.g, self.p, self.q, verbose)
        self.F, self.f = sieve_functions(k, self.alpha, self.beta, self.g)
        if verbose:
            s0 = round((max(self.alpha, self.beta) + 1.5) / H) * H
            for s in (s0, s0 + 3, s0 + 8):
                P, Q = phis(k, self.F, self.f, self.p, self.q, s)
                print(f"     Phi_P({s:.2f}) = {P:.8f}, Phi_Q = {Q:.2e};  F({s:.2f}) = {self.Fs(s):.8f}, f = {self.fs(s):.8f}")

    def Fs(self, s):
        s = np.asarray(s, dtype=float)
        return np.where(s >= grid[-1], 1.0, np.interp(s, grid, self.F))

    def fs(self, s):
        s = np.asarray(s, dtype=float)
        return np.where(s >= grid[-1], 1.0, np.interp(s, grid, self.f))


# ---------------------------------------------------------------- the weighted sieve for n(N-n)
def dhr_LR(D, theta, a, b, nu=600):
    L = float(D.fs(theta / a))
    t = np.linspace(a, b, nu + 1)
    R = 2 * float(np.trapz((1 - t / b) * D.Fs((theta - t) / a) / t, t))
    return L, R


def scan(D, theta):
    rows = []
    for a in np.arange(0.02, theta / D.beta - 1e-9, 0.002):
        for b in np.arange(0.20, theta - 0.0049, 0.005):
            L, R = dhr_LR(D, theta, a, b)
            if L <= 0:
                continue
            rows.append((R / L + 2 / b, R / L + 1 / b, a, b, L, R))
    return rows


C2H = 1.66485   # the constant c_{2h} of Part I, Proposition 4.3


def target(D, theta, k, kp, rational=None):
    """Largest margin L2 - lam R2 over (a, b) for the pair (k, k'): lam >= max(1/(k+1-2/b), 1/(k'+1-1/b))."""
    best = None
    for a in np.arange(0.02, theta / D.beta - 1e-9, 0.002):
        for b in np.arange(0.20, theta - 0.0049, 0.005):
            if k + 1 - 2 / b <= 0 or kp + 1 - 1 / b <= 0:
                continue
            lam = max(1 / (k + 1 - 2 / b), 1 / (kp + 1 - 1 / b))
            L, R = dhr_LR(D, theta, a, b)
            net = L - lam * R
            if best is None or net > best[0]:
                best = (net, a, b, lam, L, R)
    net, a, b, lam, L, R = best
    print(f"   target ({k},{kp}): max margin L2 - lam R2 = {net:+.4f} at a = {a:.3f}, b = {b:.3f}, lam = {lam:.4f} "
          f"(L2 = {L:.4f}, R2 = {R:.4f}); weight w(n) <= 2^{k + 4}, constant c = net 2^-(k+4) 8 c2h / a^2 = "
          f"{net * 2.0 ** -(k + 4) * 8 * C2H / a ** 2:.3f}" if net > 0 else
          f"   target ({k},{kp}): NOT reachable, best margin {net:+.4f} at a = {a:.3f}, b = {b:.3f}")
    if rational:
        a, b, lam = rational
        L, R = dhr_LR(D, theta, a, b, nu=2000)
        net = L - lam * R
        print(f"      rational parameters a = {a:.4f}, b = {b:.4f}, lam = {lam:.4f}: 1/lam + 2/b = {1 / lam + 2 / b:.4f}, "
              f"1/lam + 1/b = {1 / lam + 1 / b:.4f}, L2 = {L:.5f}, R2 = {R:.5f}, net = {net:+.5f}, "
              f"c = {net * 2.0 ** -(k + 4) * 8 * C2H / a ** 2:.3f}")


if __name__ == "__main__":
    thetas = [float(eval(t)) for t in sys.argv[1:]] or [0.5, 2 / 3, 0.95]
    print("linear sieve check (kappa = 1): expect alpha = beta = 2, F(3) = 2e^g/3 = %.6f, f(4) = 2e^g log 3/4 = %.6f"
          % (2 * EG / 3, 2 * EG * log(3) / 4))
    D1 = DHR(1, guess=(2.0, 2.0))
    print(f"     F(3) = {D1.Fs(3.0):.6f}, f(4) = {D1.fs(4.0):.6f}, f(2.5) = {D1.fs(2.5):.6f} (exact {2 * EG * log(1.5) / 2.5:.6f})")
    print("dimension 2:")
    D2 = DHR(2)
    print("     F_2, f_2 at s = 2, 3, 4, 4.5, 5, 6, 8, 10:")
    for s in (2, 3, 4, 4.5, 5, 6, 8, 10):
        print(f"       s = {s:4}: F_2 = {float(D2.Fs(s)):.6f}   f_2 = {float(D2.fs(s)):.6f}")
    for theta in thetas:
        rows = scan(D2, theta)
        rows.sort()
        print(f"\ntheta = {theta:.4f}: dimension-2 weighted sieve, best 1/lam + 2/b (joint) and 1/lam + 1/b (side):")
        for jt, sd, a, b, L, R in rows[:6]:
            print(f"   joint {jt:7.4f}  side {sd:7.4f}   a = {a:.3f}  b = {b:.3f}  L2 = {L:.4f}  R2 = {R:.4f}")
        rows.sort(key=lambda r: r[1])
        for jt, sd, a, b, L, R in rows[:4]:
            print(f"   joint {jt:7.4f}  side {sd:7.4f}   a = {a:.3f}  b = {b:.3f}  L2 = {L:.4f}  R2 = {R:.4f}   (best side)")
        if abs(theta - 0.5) < 1e-9:
            target(D2, theta, 7, 5, rational=(3 / 50, 37 / 100, 2 / 5))
            target(D2, theta, 7, 4)
            target(D2, theta, 6, 5)
        elif abs(theta - 2 / 3) < 1e-9:
            target(D2, theta, 6, 4)
            target(D2, theta, 5, 4)
            target(D2, theta, 6, 3)
        elif abs(theta - 0.95) < 1e-9:
            target(D2, theta, 5, 3)
            target(D2, theta, 4, 3)
            target(D2, theta, 5, 2)
