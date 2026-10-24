r"""
Figure 1 of Part II: (a) the running Rankin-Selberg estimate of ||G_{m,m'}||^2, normalised by (log mm')^2, as a
function of y/L (it reaches a plateau once y >> L, which is the Y^2 growth of Theorem A); (b) ||G_{m,m'}||^2_est
at Y = 2^21 against mm' on log-log axes, with the lines 0.365 (log mm')^2 (the fit of Section 5, consistent with
Theorem 1: ||G||^2 << (mm')^eps), 3 (mm')^{1/4} and 0.79 (mm')^{1/2} (the growth allowed by Waibel's general bound,
through the point (1,101)), and the exact value (64/9)||f_20||^2 at (1,5) (norm_level20.py).

Same tables and FFT convolution as mean_square_check.py (sup_error_scan.py of Part I); admissible N only
((N, mm') = 1, v_2(N) <= 3).  Floating point, not certified.
Usage:  python plot_mean_square.py [X]
"""
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sup_error_scan import r_table, sigma_over_n_odd, delta2_table, T_all, M_all, factor

X = int(sys.argv[1]) if len(sys.argv) > 1 else 2 ** 21
reliable = [(1, 101), (13, 29), (1, 1009), (61, 67), (101, 103)]            # no small prime factors, Y >= 50 L
biased = [(1, 5), (9, 1), (3, 5), (1, 15015)]                                # 25%-64% of the N excluded
transition = [(1009, 1013)]                                                 # Y < L
curves = [(1, 101), (13, 29), (1, 1009), (101, 103), (1009, 1013)]
colors = ["#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#c9a400"]

t0 = time.time()
r = r_table(X)
f_sigma, v2 = sigma_over_n_odd(X)
n_all = np.arange(X + 1)
print(f"X = 2^{X.bit_length() - 1}; tables in {time.time() - t0:.1f}s")


def norm_estimates(m, mp):
    """Returns (L, mask, cumulative estimate ||G||^2_est(y) on the admissible N <= y, as arrays over y = n)."""
    L = 4 * m * mp
    T = T_all(r, m, mp, X)
    M = M_all(m, mp, X, f_sigma, v2, delta2_table(m, mp))
    ok = (n_all >= 1) & (v2 <= 3)
    for p in factor(m * mp):
        ok &= n_all % p != 0
    vol = (np.pi / 3) * L * np.prod([1 + 1 / p for p in factor(L)])
    q = np.where(ok, (T - M) ** 2 / np.maximum(n_all, 1), 0.0)
    cum = np.cumsum(q) / np.maximum(np.cumsum(ok), 1)
    return L, ok, vol * cum / (16 * np.pi ** 2)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.0))

# (a) running estimate against y/L
ys = np.unique(np.round(2 ** np.arange(10, np.log2(X) + 1e-9, 0.125)).astype(int))
for (m, mp), col in zip(curves, colors):
    L, ok, G2 = norm_estimates(m, mp)
    ax1.plot(ys / L, G2[ys] / np.log(m * mp) ** 2, color=col, lw=1.6, label=rf"$(m,m')=({m},{mp})$")
    print(f"({m},{mp}): L={L}, ||G||^2_est(Y)/(log mm')^2 = {G2[X] / np.log(m * mp) ** 2:.3f}")
ax1.axhline(0.365, color="#eb6834", lw=0.8, ls="--")
ax1.text(2e-4, 0.375, r"$0.365$", fontsize=8, color="#eb6834")
ax1.axvline(1.0, color="#777", lw=0.8, ls=":")
ax1.text(1.05, 0.03, r"$y=L$", fontsize=8, color="#555")
ax1.set_xscale("log")
ax1.set_xlabel(r"$y/L$")
ax1.set_ylabel(r"$\|G\|^2_{\mathrm{est}}(y)\,/\,(\log mm')^{2}$")
ax1.set_ylim(0, 0.6)
ax1.grid(True, alpha=0.25)
ax1.legend(fontsize=7.5, loc="lower right")
ax1.set_title("(a) running estimate of the norm", fontsize=9)

# (b) norm estimate against mm'
for group, marker, face, lab in ((reliable, "o", "#2a78d6", r"$mm'$ without small prime factors, $Y\geq 50L$"),
                                 (biased, "s", "none", r"$25\%$–$64\%$ of the $N$ excluded"),
                                 (transition, "^", "none", r"$Y<L$")):
    xs, gs = [], []
    for (m, mp) in group:
        L, ok, G2 = norm_estimates(m, mp)
        xs.append(m * mp); gs.append(G2[X])
    ax2.scatter(xs, gs, marker=marker, s=38, facecolors=face, edgecolors="#2a78d6", linewidths=1.2, label=lab, zorder=3)
mm = np.logspace(np.log10(3), np.log10(2e6), 50)
ax2.plot(mm, 0.365 * np.log(mm) ** 2, color="#eb6834", lw=1.4, label=r"$0.365\,(\log mm')^2$")
ax2.plot(mm, 3 * mm ** 0.25, color="#1baf7a", lw=1.0, ls="-.", label=r"$3\,(mm')^{1/4}$")
ax2.plot(mm, 0.79 * mm ** 0.5, color="#777", lw=1.2, ls="--", label=r"$0.79\,(mm')^{1/2}$")
ax2.scatter([5], [64 / 9 * 0.0813495], marker="*", s=90, color="#eb6834", zorder=4, label=r"$(1,5)$ exact: $\frac{64}{9}\|f_{20}\|^2=0.5785$")
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xlabel(r"$mm'$")
ax2.set_ylabel(r"$\|G_{m,m'}\|^2_{\mathrm{est}}$  ($Y=2^{21}$)")
ax2.grid(True, which="both", alpha=0.25)
ax2.legend(fontsize=7.5, loc="upper left")
ax2.set_title("(b) the norm against the level", fontsize=9)

fig.tight_layout()
fig.savefig("../figures/mean_square.pdf")
fig.savefig("../figures/mean_square.png", dpi=150)
print("wrote ../figures/mean_square.pdf/.png")
