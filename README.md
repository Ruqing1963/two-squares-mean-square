# Sums of two squares with few prime factors, II: the mean square of the cuspidal error of θ²(mz)θ²(m′z), almost all N, and the size of the Petersson norm

Code, data and manuscript (draft) for

> R. Chen, *Sums of two squares with few prime factors, II: the mean square of the cuspidal error of
> θ²(mz)θ²(m′z), almost all N, and the size of the Petersson norm*, preprint, 2026.

Part I is [two-squares-theta-level](https://github.com/Ruqing1963/two-squares-theta-level)
(DOI [10.5281/zenodo.23618710](https://doi.org/10.5281/zenodo.23618710)).

## What the paper does

Let m, m′ be coprime odd integers, Θ_{m,m′}(z) = θ²(mz)θ²(m′z), G = G_{m,m′} its cuspidal part and L = 4mm′, so that
a_G(N) = T(m,m′;N) − M(m,m′;N) is the error term in the number of representations of N by m(x²+y²)+m′(z²+w²).

* **Theorem 1.1 (recalled from Part I, version 3).** ‖G_{m,m′}‖² ≪ (mm′)^ε and T − M ≪ d(N) N^{1/2} (mm′)^ε
  uniformly; Part I deduces the level of distribution 1/2 for r(n) r(N − n) and Ω(n(N − n)) ≤ 8,
  max(Ω(n), Ω(N − n)) ≤ 5 for every large N ≡ 2 (mod 4). (Until the restructuring of 2026-10-24 the norm bound and
  its proof were the main result of this paper; they now form Theorem 1.2 and Section 2.5 of Part I.)
* **Theorem 1.2 (mean square).** For every G ∈ S₂(Γ₀(L)) and Y ≥ 1, Σ_{n ≤ Y} |a_G(n)|² ≤ 4πe² Y (1 + 54Y/L) ‖G‖²
  (covering argument); hence Σ_{N ≤ Y} |T − M|² ≪ Y² (mm′)^{−1+ε} for Y ≥ mm′.
* **Theorem 1.3 / Corollary 1.4.** Level of distribution 1 − ε for all but O(X^{1−δ}) integers N ≡ 2 (mod 4) up
  to X, and for those N: Ω(n(N − n)) ≤ 5 and max(Ω(n), Ω(N − n)) ≤ 4, constant 0.34 (parameters
  θ = 19/20, α = 2/25, λ = 1/3, β₁ = 69/100 certified in interval arithmetic).
* For almost all N the circle method (Montgomery–Vaughan with characters mod 4) gives N = p + q with
  p ≡ q ≡ 1 (mod 4), so Corollary 1.4 is not new as an arithmetic statement; the new content is Theorem 1.3
  (the level 1 − ε for the weighted sequence) and the every-N results. Section 4 states the individual bound
  |T − M| ≪ N^{1/2+ε}(mm′)^{−1/2+ε} (Conjecture 4.1) that would give level 1 − ε and (5,4) for every N
  (Proposition 4.2), and Section 5.4 gives the supremum data behind it (sup/rms between 5 and 10 at every level).
* Remark 3.2: the vector sieve of Part I has sieving limit 2 + 2√e = 5.30; the product n(N − n) is a genuine
  dimension-2 sieve problem (limit 4.27), which should improve the pairs (k, k′) without new analytic input.
* **Theorem 1.5 (lower bound for the norm, Section 6).** ‖G_{m,m′}‖² ≥ (1/4π) Σ_{Q ∈ {m,4m,m′,4m′}} Σ_n |a_{G_{1,mm′}}(n)|² n⁻¹ e^{−4πn/Q}
  (the regions above the four Atkin–Lehner cusps at which G_{m,m′} becomes G_{1,mm′}), hence
  ‖G_{m,m′}‖² ≥ (1/π)((log m)² + (log m′)²) − O(log(2mm′)(log log 3mm′)^C): the exponent ε in the upper bound
  (mm′)^ε of Part I cannot be removed. Numerically (Table 3, `cusp_lower_bound.py`) the regions above all cusps
  carry 93–98 % of the norm, and ‖G‖² ≈ (1/π)((log m)² + (log m′)²) + (5/π²) log m log m′ + (terms for the other
  divisors of mm′): the "0.37 (log mm′)²" of Section 5 is the average of this over the sample.
* Numerics: the Rankin–Selberg estimate of ‖G‖² from the data satisfies ‖G_{m,m′}‖² ≈ 0.37 (log mm′)² for
  10² ≤ mm′ ≤ 1.6·10⁴; an exact check at level 20 (G_{1,5} = (8/3) f_20, ‖f_20‖² = 0.08135 from the period lattice
  of 20a1 and from Rankin–Selberg) validates the whole numerical chain. `local_densities.py` makes the Siegel main
  term exact for every n (square-free m, m′), checked against the levels 4, 12 (no cusp forms: T = M exactly) and
  20 (T − M = (8/3) a_f(n) for all n); with it the unbiased estimate for (1,5) is 0.5785, the exact value.

## Repository layout

```
paper/      two_squares_mean_square.tex, two_squares_mean_square.pdf
code/       mean_square_check.py      mean square of T - M over N <= Y and the estimate of ||G||^2 (Table 1)
            norm_level20.py           exact check at level 20: ||f_20||^2 from the eta product (Rankin-Selberg) and from the period lattice
            certify_modular.py        Part I's interval-arithmetic certification, extended by the arguments 2/3 and 19/20
            sup_error_scan.py         Part I: T(m,m';N) by FFT and the exact Siegel main term (used by mean_square_check.py)
            modular_sieve_scan.py     Part I: floating-point parameter scan (uses sieve_core.py)
            sieve_core.py             linear sieve functions
            plot_mean_square.py       Figure 1: running norm estimate against y/L, and ||G||^2 against mm'
            plot_sieve_thresholds.py  Figure 2: best sieve thresholds against the level theta
            local_densities.py        exact Siegel main term for every n (and for the odd vectors); self-test at levels 4, 12, 20
            cusp_lower_bound.py       Table 3: the lower bound from the regions above the cusps, cusp by cusp
figures/    mean_square, sieve_thresholds  (.pdf and .png)
results/    mean_square_check_x21.txt              Table 1
            cusp_lower_bound.txt, cusp_lower_bound_large.txt   Table 3 (X = 2^21; X = 2^24 for (1,15015), (1009,1013))
            norm_level20.txt                       the exact check at level 20
            certify_modular_output_twothirds.txt   certified parameters at level 2/3 (bounds 6 and 4)
            certify_modular_output_nineteen.txt    certified parameters at level 19/20 (Corollary 1.4: bounds 5 and 4)
            sieve_thresholds.txt                   data of Figure 2
```

## Reproducing the numbers

Requirements: Python ≥ 3.10 with `numpy`, `scipy`, `mpmath` (versions used: `requirements.txt`).

```
cd code
python certify_modular.py 19/20        # Corollary 1.4: level 19/20, bounds 5 and 4 (certified, directed rounding)
python certify_modular.py 2/3          # level 2/3, bounds 6 and 4 (certified)
python certify_modular.py 1/2          # level 1/2, bounds 8 and 5 (certified; the parameters of Part I, Table 2, first line)
python mean_square_check.py            # Table 1: X = 2^21, ten pairs (m, m') (~1 min, floating point)
python norm_level20.py                 # exact check at level 20 (~1 min)
python modular_sieve_scan.py 0.95      # floating-point scan at level 19/20 (not certified)
python plot_mean_square.py             # Figure 1 (~2 min)
python plot_sieve_thresholds.py        # Figure 2 from results/sieve_thresholds.txt (--recompute: several minutes)
python local_densities.py              # self-test of the exact main term (levels 4, 12, 20; ~10 s)
python cusp_lower_bound.py --detail    # Table 3 (~2 min); python cusp_lower_bound.py 16777216 1,15015 1009,1013 for the two large pairs
```

What is rigorous: `certify_modular.py` (interval arithmetic and directed rounding) enters the proofs of
Corollary 1.4 (and, at the level 1/2, Theorem 1.1 of Part I); `mean_square_check.py`, `norm_level20.py` and the
plots are floating point and are evidence only (Section 5 of the paper).

## License

Code: MIT (see `LICENSE`). Manuscript (`paper/`): CC BY 4.0.
