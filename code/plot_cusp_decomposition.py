r"""
Figure 3 of Part II: the Petersson norm decomposed over the cusps (Section 6.4, Table 3).
(a) For each pair (m,m') with Y >= 35 L, the parts of the lower bound Lambda(m,m') as fractions of the Rankin-Selberg
    estimate ||G||^2: the four cusps of Theorem 1.5(i), the other Atkin-Lehner cusps, the cusps of denominator 2,
    and the rest ||G||^2 - Lambda (the part of the fundamental domain below height 1); dashed: 3/pi, the share of the
    area.  (b) Lambda and the leading term (1/pi)(log^2 m + log^2 m') + (5/pi^2) log m log m' against ||G||^2.
Reads results/cusp_lower_bound.txt.   Usage: python plot_cusp_decomposition.py
"""
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import pi, log

BLUE, BLUE_L, ORANGE, GREY, INK, MUTED = "#2a78d6", "#86b6ef", "#eb6834", "#c3c2b7", "#0b0b0b", "#898781"
rows = []
pat = re.compile(r"^\s*\((\d+), (\d+)\)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)")
with open("../results/cusp_lower_bound.txt") as fh:
    for line in fh:
        mt = pat.match(line)
        if mt:
            m, mp, mm, g2, adm, lb_al, lb_odd, lb, ratio, thm, asym = (float(x) for x in mt.groups())
            rows.append((int(m), int(mp), g2, lb_al, lb_odd, lb, thm))
rows.sort(key=lambda r: r[0] * r[1])

fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.3), gridspec_kw=dict(width_ratios=[1.45, 1]))
labels = [f"({m},{mp})" for m, mp, *_ in rows]
x = np.arange(len(rows))
thm = np.array([r[6] / r[2] for r in rows])
oth = np.array([(r[3] - r[6]) / r[2] for r in rows])
odd = np.array([r[4] / r[2] for r in rows])
rest = np.array([1 - r[5] / r[2] for r in rows])
w = 0.72
ax.bar(x, thm, w, color=BLUE, label="four cusps of Theorem 1.5(i)")
ax.bar(x, oth, w, bottom=thm, color=BLUE_L, label="other Atkin-Lehner cusps", edgecolor="white", linewidth=0.6)
ax.bar(x, odd, w, bottom=thm + oth, color=ORANGE, label="cusps of denominator 2", edgecolor="white", linewidth=0.6)
ax.bar(x, rest, w, bottom=thm + oth + odd, color=GREY, label=r"rest: $\|G\|^2-\Lambda$ (below height 1)",
       edgecolor="white", linewidth=0.6)
ax.axhline(3 / pi, color=INK, linestyle="--", linewidth=1, label=r"$3/\pi=95.5\%$, the share of the area")
for i, r in enumerate(rows):
    ax.text(i, 1.012, f"{100 * r[5] / r[2]:.0f}%", ha="center", va="bottom", fontsize=7.5, color=INK)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
ax.set_ylim(0, 1.07)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
ax.set_ylabel(r"fraction of $\|G_{m,m'}\|^2$")
ax.set_title(r"(a) the norm, cusp by cusp (pairs ordered by $mm'$)", fontsize=10)
ax.legend(fontsize=7.8, loc="lower left", framealpha=0.95)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

g2 = np.array([r[2] for r in rows])
lb = np.array([r[5] for r in rows])
lead = np.array([(log(r[0]) ** 2 + log(r[1]) ** 2) / pi + 5 / pi ** 2 * log(r[0]) * log(r[1]) for r in rows])
fit = np.array([0.365 * log(r[0] * r[1]) ** 2 for r in rows])
lim = [0.3, 60]
bx.plot(lim, lim, color=INK, linewidth=1, linestyle="--", label="equality")
bx.scatter(g2, lb, s=34, color=BLUE, zorder=3, label=r"$\Lambda(m,m')$, all cusps")
bx.scatter(g2, lead, s=34, color=ORANGE, marker="s", zorder=3,
           label=r"$\frac{1}{\pi}(\log^2 m+\log^2 m')+\frac{5}{\pi^2}\log m\,\log m'$")
bx.scatter(g2, fit, s=30, facecolors="none", edgecolors=MUTED, marker="o", zorder=2, label=r"$0.365\,(\log mm')^2$")
for r, gg, ll in zip(rows, g2, lb):
    if r[0] * r[1] in (5, 377, 1155, 1009, 10403, 15015):
        bx.annotate(f"({r[0]},{r[1]})", (gg, ll), xytext=(4, -9), textcoords="offset points", fontsize=7, color=INK)
bx.set_xscale("log")
bx.set_yscale("log")
bx.set_xlim(*lim)
bx.set_ylim(*lim)
bx.set_xlabel(r"$\|G_{m,m'}\|^2$ (Rankin-Selberg estimate, all $n\leq2^{21}$)")
bx.set_ylabel("lower bound / formula")
bx.set_title("(b) the lower bound and the leading term", fontsize=10)
bx.legend(fontsize=7.8, loc="upper left", framealpha=0.95)
bx.grid(True, which="major", color="#e1e0d9", linewidth=0.6)
for s in ("top", "right"):
    bx.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("../figures/cusp_decomposition.pdf")
fig.savefig("../figures/cusp_decomposition.png", dpi=160)
print("written figures/cusp_decomposition.pdf/.png;", len(rows), "pairs")
