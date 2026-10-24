r"""
Exact check at level 20 of the whole numerical chain of Part II (Section 5).

S_2(Gamma_0(20)) is spanned by the newform f_20 = eta(2z)^2 eta(10z)^2 of the elliptic curve 20a1,
y^2 = x^3 + x^2 + 4x + 4 (modular degree 1, Manin constant 1), so G_{1,5} = c f_20 with
c = a_G(1) = T(1,5;1) - M(1,5;1) = 4 - 4/3 = 8/3.  The script
  (1) computes a_f(n), n <= 2^20, from the eta product;
  (2) estimates ||f_20||^2 by the Rankin-Selberg asymptotic  sum_{n<=Y} a_f(n)^2/n ~ 16 pi^2 ||f||^2 Y / vol,
      vol(Gamma_0(20)\H) = (pi/3) * 20 * (3/2)(6/5) = 12 pi;
  (3) computes ||f_20||^2 = covol(Lambda) / (4 pi^2) exactly from the period lattice Lambda of 20a1
      (phi^* omega = 2 pi i f(z) dz, deg phi = 1), the lattice being validated through g_2, g_3;
  (4) checks that T - M of sup_error_scan.py equals (8/3) a_f(N) on all admissible N.
The three values of ||f||^2 agree to four decimals: 0.08135.  Hence ||G_{1,5}||^2 = (64/9)||f_20||^2 = 0.5785.
Usage:  python norm_level20.py
"""
import numpy as np
import mpmath as mp
from scipy.signal import fftconvolve
from sup_error_scan import r_table, sigma_over_n_odd, delta2_table, T_all, M_all

Y = 2 ** 20


def euler(Xmax, step):
    """prod_{n>=1} (1 - q^{step n}) as a power series (Euler's pentagonal theorem)."""
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


# (1) f_20 = q P(q^2)^2 P(q^10)^2
F = np.round(fftconvolve(sq(euler(Y, 2)), sq(euler(Y, 10)))[: Y + 1])
a = np.zeros(Y + 1)
a[1:] = F[:Y]
print("a_f(n), n = 1..13:", a[1:14].astype(int), " (20a1: a2=0, a3=-2, a5=-1, a7=2, a11=0, a13=2)")

# (2) Rankin-Selberg estimate
n = np.arange(Y + 1)
vol = 12 * np.pi
for y in (2 ** 16, 2 ** 18, 2 ** 20):
    S = np.sum(a[1:y + 1] ** 2 / n[1:y + 1])
    print(f"Rankin-Selberg, y = 2^{y.bit_length() - 1}: ||f_20||^2 = {S * vol / (16 * np.pi ** 2 * y):.6f}")

# (3) period lattice of y^2 = (x+1)(x^2+4), omega = dx/(2y)
mp.mp.dps = 30
cubic = lambda x: (x + 1) * (x * x + 4)
w1 = 2 * mp.quad(lambda x: 1 / (2 * mp.sqrt(cubic(x))), [-1, mp.inf])                       # real period
w2 = 2 * mp.quad(lambda t: 1j / (2 * mp.sqrt(4 - t * t) * mp.sqrt(1 + 1j * t)), [-2, 2])    # x = i t, -2 <= t <= 2
p, q = 4 - mp.mpf(1) / 3, 4 - mp.mpf(4) / 3 + mp.mpf(2) / 27    # x = X - 1/3:  X^3 + p X + q, so g2 = -4p, g3 = -4q


def eisenstein(tau):
    qq = mp.exp(2j * mp.pi * tau)
    E4, E6 = mp.mpf(1), mp.mpf(1)
    for k in range(1, 400):
        s3 = sum(d ** 3 for d in range(1, k + 1) if k % d == 0)
        s5 = sum(d ** 5 for d in range(1, k + 1) if k % d == 0)
        E4 += 240 * s3 * qq ** k
        E6 -= 504 * s5 * qq ** k
    return E4, E6


print(f"periods: omega1 = {mp.nstr(w1, 12)}, omega2 = {mp.nstr(w2, 12)};  curve: g2 = {mp.nstr(-4 * p, 10)}, g3 = {mp.nstr(-4 * q, 10)}")
for name, om1, om2 in (("Z w1 + Z w2", w1, w2), ("Z w1 + Z (w1+w2)/2", w1, (w1 + w2) / 2)):
    tau = om2 / om1
    E4, E6 = eisenstein(tau)
    g2 = (2 * mp.pi / om1) ** 4 * E4 / 12
    g3 = -(2 * mp.pi / om1) ** 6 * E6 / 216
    cov = abs(mp.im(om1 * mp.conj(om2)))
    print(f"  lattice {name}: g2 = {mp.nstr(mp.re(g2), 10)}, g3 = {mp.nstr(mp.re(g3), 10)}, covolume = {mp.nstr(cov, 10)},"
          f" ||f_20||^2 = covol/(4 pi^2) = {mp.nstr(cov / (4 * mp.pi ** 2), 10)}")
print("  (the lattice with g2, |g3| equal to those of the curve is the period lattice up to the sign of g3, i.e. up to the twist by -1, which has the same covolume)")

# (4) the pipeline of Part I against (8/3) a_f(N)
r = r_table(Y)
f_sigma, v2 = sigma_over_n_odd(Y)
T = T_all(r, 1, 5, Y)
M = M_all(1, 5, Y, f_sigma, v2, delta2_table(1, 5))
ok = (n >= 1) & (v2 <= 3) & (n % 5 != 0)
print("pipeline: max |(T - M) - (8/3) a_f(N)| over admissible N =", np.max(np.abs((T - M)[ok] - (8 / 3) * a[ok])))
print("||G_{1,5}||^2 = (64/9) ||f_20||^2 = ", 64 / 9 * 0.0813495, "; estimate over admissible N only (biased, Table 1):",
      vol * np.mean((T - M)[ok] ** 2 / n[ok]) / (16 * np.pi ** 2))
