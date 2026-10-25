r"""
Rigorous (outward-rounded) enclosures of the Diamond-Halberstam-Richert sieve functions F_2, f_2 and certification of
the two inequalities of Part II, Theorem 7.1:

    f_2(9) >= 0.9995,     R_2 = 2 int_{1/18}^{2/5} (1 - 5t/2) F_2(9 - 18t) dt/t = 2 int_{1.8}^{8} (s - 1.8)/(7.2 (9 - s)) F_2(s) ds <= 2.625.

Inputs taken as rigorous:
  * alpha_2, beta_2: the values of Booker-Browning (arXiv:1511.00601, ancillary dhr.html, computed in ball arithmetic),
    "correctly truncated at the 20th decimal place": alpha_2 in [5.35772744559446184227, +1e-20],
    beta_2 in [4.26645028414864191641, +1e-20];
  * sigma_2(s) = A s^2 on (0,2], A = 1/(2 (2e^gamma)^2); (s^-2 sigma)' = -2 s^-3 sigma(s-2) for s > 2, which gives the
    closed form g(s) := s^-2 sigma(s) = A (4 - 2 log(s/2) - 8/s + 4/s^2) on [2,4]; on [4, alpha_2] one more
    integration, g(s) = g(4) - 2 int_4^s u^-3 (u-2)^2 g(u-2) du, enclosed cell by cell (u^-3 (u-2)^2 increases for
    u < 6, g decreases);
  * F_2 = 1/sigma_2 on (0, alpha_2], f_2 = 0 on (0, beta_2], (s^2 F_2)' = 2 s f_2(s-1), (s^2 f_2)' = 2 s F_2(s-1);
    F_2 is decreasing and f_2 increasing (Diamond-Halberstam, Lemma 6.2), which bounds each cell integral by its
    endpoint values, as in Part I, Appendix A (certify_modular.py).
Closed forms in mpmath interval arithmetic (80 bits); the recursions in IEEE doubles with one outward nextafter step
after each operation.  Grid step 1/H_DEN.   Usage: python certify_dhr.py [H_DEN]
"""
import math
import sys
import numpy as np
from mpmath import iv
iv.prec = 80

dn = lambda x: math.nextafter(x, -math.inf)
up = lambda x: math.nextafter(x, math.inf)
H_DEN = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
H = 1.0 / H_DEN
S_TOP = 9

ALPHA = iv.mpf(['5.35772744559446184227', '5.35772744559446184228'])
BETA = iv.mpf(['4.26645028414864191641', '4.26645028414864191642'])
EG = iv.exp(iv.euler)
A = 1 / (2 * (2 * EG) ** 2)


def g_closed(s):
    """g(s) = s^-2 sigma_2(s) for s in (0, 4], s an interval."""
    if s.b <= 2:
        return A
    assert s.a >= 2 and s.b <= 4
    return A * (4 - 2 * iv.log(s / 2) - 8 / s + 4 / s ** 2)


def lo(v):
    return dn(float(v.a))


def hi(v):
    return up(float(v.b))


n = S_TOP * H_DEN
sg = [iv.mpf([i, i]) / H_DEN for i in range(n + 1)]              # exact grid points as intervals
pl = [lo(v) for v in sg]
ph = [hi(v) for v in sg]

# ---- sigma_2 on (0, alpha_2]: g_lo[i] <= g(s_i) <= g_hi[i] for s_i <= alpha_2 ----
g_lo = np.zeros(n + 1)
g_hi = np.zeros(n + 1)
i4 = 4 * H_DEN
for i in range(1, i4 + 1):
    v = g_closed(sg[i])
    g_lo[i], g_hi[i] = lo(v), hi(v)
ia = int(math.floor(float(ALPHA.a) * H_DEN))                      # last grid index with s_i <= alpha_2
# g(s) = g(4) - int_4^s phi(u) du, phi(u) = 2 u^-3 (u-2)^2 g(u-2): phi increases in the first factor, g decreases
phi_fac = lambda u: 2 * (u - 2) ** 2 / u ** 3                       # u interval
I_lo, I_hi = g_lo[i4], g_hi[i4]
for i in range(i4, ia):
    fa, fb = phi_fac(sg[i]), phi_fac(sg[i + 1])                   # increasing in u on [4, 6]
    j = i - 2 * H_DEN                                              # index of u - 2
    cell_lo = dn(dn(lo(fa) * g_lo[j + 1]) * dn(pl[i + 1] - ph[i]))   # smallest value of phi times the width
    cell_hi = up(up(hi(fb) * g_hi[j]) * up(ph[i + 1] - pl[i]))
    I_lo, I_hi = dn(I_lo - cell_hi), up(I_hi - cell_lo)
    g_lo[i + 1], g_hi[i + 1] = I_lo, I_hi
# g at alpha_2 itself (partial cell from s_ia to alpha_2)
fa, fb = phi_fac(sg[ia]), phi_fac(ALPHA)
j = ia - 2 * H_DEN
w_lo, w_hi = dn(float(ALPHA.a) - ph[ia]), up(float(ALPHA.b) - pl[ia])
gA_lo = dn(g_lo[ia] - up(up(hi(fb) * g_hi[j]) * w_hi))
gA_hi = up(g_hi[ia] - dn(dn(lo(fa) * g_lo[j + 1]) * w_lo))
print(f"sigma_2: g(4) in [{g_lo[i4]:.10f}, {g_hi[i4]:.10f}], sigma_2(4) in [{16 * g_lo[i4]:.8f}, {16 * g_hi[i4]:.8f}];"
      f" sigma_2(alpha_2) in [{dn(gA_lo * dn(float(ALPHA.a) ** 2)):.8f}, {up(gA_hi * up(float(ALPHA.b) ** 2)):.8f}]")

# ---- F_2 and f_2 on the grid ----
F_lo = np.zeros(n + 1)
F_hi = np.full(n + 1, math.inf)
f_lo = np.zeros(n + 1)
f_hi = np.zeros(n + 1)
for i in range(1, ia + 1):
    F_lo[i] = dn(1 / up(up(ph[i] ** 2) * g_hi[i]))
    F_hi[i] = up(1 / dn(dn(pl[i] ** 2) * g_lo[i]))
# u = s^2 F:  u(alpha_2) = 1/g(alpha_2)
u_lo, u_hi = dn(1 / gA_hi), up(1 / gA_lo)
ib = int(math.floor(float(BETA.a) * H_DEN))                       # last grid index with s_i <= beta_2
v_lo, v_hi = 0.0, 0.0                                              # v = s^2 f, v(beta_2) = 0
# first partial cells: f on [beta_2, s_{ib+1}] uses F(x-1) with x-1 <= s_{ib+1}-1 < alpha_2 (grid or sigma);
# F on [alpha_2, s_{ia+1}] uses f(x-1), x-1 in [alpha_2-1, s_{ia+1}-1], already available (beta_2 < alpha_2 - 1).
def f_lo_at(x_lo_index):
    return f_lo[x_lo_index]          # f increasing: value at the grid point below is a lower bound


# process in increasing s; f first (its cells start earlier)
# -- f: partial cell [beta_2, s_{ib+1}]
wb_lo, wb_hi = dn(pl[ib + 1] - hi(BETA)), up(ph[ib + 1] - lo(BETA))
jb = ib + 1 - H_DEN                                                # index of s_{ib+1} - 1
v_lo = dn(dn(dn(2 * lo(BETA)) * F_lo[jb]) * wb_lo)                 # F decreasing: min at the right end x = s_{ib+1}
v_hi = up(up(up(2 * ph[ib + 1]) * F_hi[jb - 1]) * wb_hi)           # max at the left end x = beta_2 >= s_{ib}: F(beta_2-1) <= F(s_{ib}-1)
f_lo[ib + 1], f_hi[ib + 1] = dn(v_lo / up(ph[ib + 1] ** 2)), up(v_hi / dn(pl[ib + 1] ** 2))
for i in range(ib + 1, n):
    j = i - H_DEN                                                  # s_i - 1
    w_lo, w_hi = dn(pl[i + 1] - ph[i]), up(ph[i + 1] - pl[i])
    v_lo = dn(v_lo + dn(dn(dn(2 * pl[i]) * F_lo[j + 1]) * w_lo))
    v_hi = up(v_hi + up(up(up(2 * ph[i + 1]) * F_hi[j]) * w_hi))
    f_lo[i + 1], f_hi[i + 1] = dn(v_lo / up(ph[i + 1] ** 2)), up(v_hi / dn(pl[i + 1] ** 2))
    if i + 1 == ia + 1:
        # -- F: partial cell [alpha_2, s_{ia+1}], then the regular cells
        wa_lo, wa_hi = dn(pl[ia + 1] - hi(ALPHA)), up(ph[ia + 1] - lo(ALPHA))
        ja = ia - H_DEN                                            # grid point below alpha_2 - 1
        u_lo = dn(u_lo + dn(dn(dn(2 * lo(ALPHA)) * f_lo[ja]) * wa_lo))        # f increasing: min at x = alpha_2, f(alpha_2-1) >= f(s_ia - 1)
        u_hi = up(u_hi + up(up(up(2 * ph[ia + 1]) * f_hi[ja + 1]) * wa_hi))   # max at x = s_{ia+1}
        F_lo[ia + 1], F_hi[ia + 1] = dn(u_lo / up(ph[ia + 1] ** 2)), up(u_hi / dn(pl[ia + 1] ** 2))
    elif i + 1 > ia + 1:
        u_lo = dn(u_lo + dn(dn(dn(2 * pl[i]) * f_lo[j]) * w_lo))
        u_hi = up(u_hi + up(up(up(2 * ph[i + 1]) * f_hi[j + 1]) * w_hi))
        F_lo[i + 1], F_hi[i + 1] = dn(u_lo / up(ph[i + 1] ** 2)), up(u_hi / dn(pl[i + 1] ** 2))

for s in (2, 3, 4, 4.5, 5, 6, 8, 9):
    i = int(round(s * H_DEN))
    print(f"s = {s:>4}: F_2 in [{F_lo[i]:.7f}, {F_hi[i]:.7f}]   f_2 in [{f_lo[i]:.7f}, {f_hi[i]:.7f}]")

# ---- the two inequalities ----
i9 = 9 * H_DEN
f9_lo = f_lo[i9]
# R_2 = 2 int_{1.8}^{8} w(s) F_2(s) ds, w(s) = (s - 1.8)/(7.2 (9 - s)) increasing, F_2 decreasing
R_hi = 0.0
i_lo, i_hi = int(round(1.8 * H_DEN)), 8 * H_DEN
for i in range(i_lo, i_hi):
    w = (sg[i + 1] - iv.mpf('1.8')) / (iv.mpf('7.2') * (9 - sg[i + 1]))
    R_hi = up(R_hi + up(up(hi(w) * F_hi[i]) * up(ph[i + 1] - pl[i])))
R_hi = up(2 * R_hi)
net_lo = dn(f9_lo - up(R_hi / dn(3.0)))
print(f"\ncertified (grid 1/{H_DEN}):  f_2(9) >= {f9_lo:.6f}   R_2 <= {R_hi:.6f}   "
      f"L_2 - R_2/3 >= {net_lo:.6f}   constant >= {dn(dn(dn(net_lo * 2.0 ** -8) * 1.66485) * 324):.4f}")
print("Lemma 7.2 (Theorem 7.1) needs f_2(9) >= 0.9995 and R_2 <= 2.625:",
      "CERTIFIED" if f9_lo >= 0.9995 and R_hi <= 2.625 else "NOT certified")
