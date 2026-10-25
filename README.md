# Sums of two squares with few prime factors, II: the mean square of the cuspidal error of Î¸Â²(mz)Î¸Â²(mâ€²z), almost all N, and the size of the Petersson norm

Code, data and manuscript (draft) for

> R. Chen, *Sums of two squares with few prime factors, II: the mean square of the cuspidal error of
> Î¸Â²(mz)Î¸Â²(mâ€²z), almost all N, and the size of the Petersson norm*, preprint, 2026.

Part I is [two-squares-theta-level](https://github.com/Ruqing1963/two-squares-theta-level)
(DOI [10.5281/zenodo.23618710](https://doi.org/10.5281/zenodo.23618710)).

## What the paper does

Let m, mâ€² be coprime odd integers, Î˜_{m,mâ€²}(z) = Î¸Â²(mz)Î¸Â²(mâ€²z), G = G_{m,mâ€²} its cuspidal part and L = 4mmâ€², so that
a_G(N) = T(m,mâ€²;N) âˆ’ M(m,mâ€²;N) is the error term in the number of representations of N by m(xÂ²+yÂ²)+mâ€²(zÂ²+wÂ²).

* **Theorem 1.1 (recalled from Part I, version 3).** â€–G_{m,mâ€²}â€–Â² â‰ª (mmâ€²)^Îµ and T âˆ’ M â‰ª d(N) N^{1/2} (mmâ€²)^Îµ
  uniformly; Part I deduces the level of distribution 1/2 for r(n) r(N âˆ’ n) and Î©(n(N âˆ’ n)) â‰¤ 8,
  max(Î©(n), Î©(N âˆ’ n)) â‰¤ 5 for every large N â‰¡ 2 (mod 4). (Until the restructuring of 2026-10-24 the norm bound and
  its proof were the main result of this paper; they now form Theorem 1.2 and Section 2.5 of Part I.)
* **Theorem 1.2 (mean square).** For every G âˆˆ Sâ‚‚(Î“â‚€(L)) and Y â‰¥ 1, Î£_{n â‰¤ Y} |a_G(n)|Â² â‰¤ 4Ï€eÂ² Y (1 + 54Y/L) â€–Gâ€–Â²
  (covering argument); hence Î£_{N â‰¤ Y} |T âˆ’ M|Â² â‰ª YÂ² (mmâ€²)^{âˆ’1+Îµ} for Y â‰¥ mmâ€².
* **Theorem 1.3 / Corollary 1.4.** Level of distribution 1 âˆ’ Îµ for all but O(X^{1âˆ’Î´}) integers N â‰¡ 2 (mod 4) up
  to X, and for those N: Î©(n(N âˆ’ n)) â‰¤ 5 and max(Î©(n), Î©(N âˆ’ n)) â‰¤ 4, constant 0.34 (parameters
  Î¸ = 19/20, Î± = 2/25, Î» = 1/3, Î²â‚ = 69/100 certified in interval arithmetic).
* For almost all N the circle method (Montgomeryâ€“Vaughan with characters mod 4) gives N = p + q with
  p â‰¡ q â‰¡ 1 (mod 4), so Corollary 1.4 is not new as an arithmetic statement; the new content is Theorem 1.3
  (the level 1 âˆ’ Îµ for the weighted sequence) and the every-N results. Section 4 states the individual bound
  |T âˆ’ M| â‰ª N^{1/2+Îµ}(mmâ€²)^{âˆ’1/2+Îµ} (Conjecture 4.1) that would give level 1 âˆ’ Îµ and (5,4) for every N
  (Proposition 4.2), and Section 5.4 gives the supremum data behind it (sup/rms between 5 and 10 at every level).
* Remark 3.2: the vector sieve of Part I has sieving limit 2 + 2âˆše = 5.30; the product n(N âˆ’ n) is a genuine
  dimension-2 sieve problem (limit 4.27), which should improve the pairs (k, kâ€²) without new analytic input.
* **Theorem 1.5 (lower bound for the norm, Section 6).** â€–G_{m,mâ€²}â€–Â² â‰¥ (1/4Ï€) Î£_{Q âˆˆ {m,4m,mâ€²,4mâ€²}} Î£_n |a_{G_{1,mmâ€²}}(n)|Â² nâ»Â¹ e^{âˆ’4Ï€n/Q}
  (the regions above the four Atkinâ€“Lehner cusps at which G_{m,mâ€²} becomes G_{1,mmâ€²}), hence
  â€–G_{m,mâ€²}â€–Â² â‰¥ (1/Ï€)((log m)Â² + (log mâ€²)Â²) âˆ’ O(log(2mmâ€²)(log log 3mmâ€²)^C): the exponent Îµ in the upper bound
  (mmâ€²)^Îµ of Part I cannot be removed. Numerically (Table 3, `cusp_lower_bound.py`) the regions above all cusps
  carry 93â€“98 % of the norm, and â€–Gâ€–Â² â‰ˆ (1/Ï€)((log m)Â² + (log mâ€²)Â²) + (5/Ï€Â²) log m log mâ€² + (terms for the other
  divisors of mmâ€²): the "0.37 (log mmâ€²)Â²" of Section 5 is the average of this over the sample.
* **Theorem 6.5 (the exact order).** â€–G_{m,mâ€²}â€–Â² â‰ª (log 3mmâ€²)Â² (log log 3mmâ€²)^C, and â‰ª_k (log 3mmâ€²)Â² when
  Ï‰(mmâ€²) â‰¤ k; with Theorem 1.5, â€–G_{m,mâ€²}â€–Â² â‰ (log mmâ€²)Â² for mmâ€² with boundedly many prime factors, and the
  uniform bound T âˆ’ M â‰ª d(N) N^{1/2} (log 3mmâ€²)(log log 3mmâ€²)^C replaces the (mmâ€²)^Îµ of Part I. Proof: the two
  Îµ-losses of Part I (d(Î±)d(Î²) and Waibel's genus bound) are replaced by the mean value of r(Î±)Â² and by the local
  densities. Corollary 6.6: the spectral second moment Î£_Ï† |âŸ¨Î˜,Ï†âŸ©|Â²/â€–Ï†â€–Â² is â‰ (log mmâ€²)Â². For (1,p), p â‰¤ 10â¶,
  Î›(1,p) = 0.319 (log p)Â² + 0.42 log p âˆ’ 1.5 (fit), against the leading term 1/Ï€ = 0.318.
* **Theorem 6.7 (Rankinâ€“Selberg formula).** For every newform Ï† of weight 2 and level 4mmâ€² (mmâ€² odd square-free),
  âŸ¨Î˜_{m,mâ€²}, Ï†âŸ© = (4/Ï€Â²) Ï‡â‚‹â‚„(mâ€²) Î·_{mâ€²}(Ï†) L(1,Ï†) L(1,Ï†âŠ—Ï‡â‚‹â‚„) â€” all local factors cancel, and |âŸ¨Î˜_{m,mâ€²},Ï†âŸ©| is
  independent of the factorisation. Proof: Hecke's weight-1 Eisenstein series Î¸Â² unfolded over the cusps of
  Î“â‚€(4mmâ€²) above âˆž, which are the Atkinâ€“Lehner cusps W_Qâˆž. Verified numerically at levels 20, 52, 116
  (`test_rs_formula.py`: ratios 1.00000, 0.99975, 0.99997/0.99878/0.99760). Corollary 6.9: the newform second
  moment Î£ L(1,Ï†)Â²L(1,Ï†âŠ—Ï‡)Â²/âŸ¨Ï†,Ï†âŸ© = (Ï€â´/16)â€–pr_new Gâ€–Â² â‰ª (log 3mmâ€²)Â²(log log 3mmâ€²)^C.
* **Theorem 7.1 (two-dimensional sieve, Section 7).** Replacing the vector sieve of Part I by the Diamondâ€“Halberstamâ€“Richert
  sieve of dimension 2 (sieving limit Î²â‚‚ = 4.2665 instead of 2 + 2âˆše = 5.30) gives, from the same level 1/2,
  Î©(n(N âˆ’ n)) â‰¤ 7 and max(Î©(n), Î©(N âˆ’ n)) â‰¤ 5 for â‰¥ 0.25 Î£_r(N) N/(log N)Â² integers n, for every large N â‰¡ 2 (mod 4):
  (7,5) instead of (8,5). `dhr_sieve.py` computes Fâ‚‚, fâ‚‚ from the differential-delay system and determines
  Î±â‚‚ = 5.357727, Î²â‚‚ = 4.266450 (agreeing to 1e-8 with the ball-arithmetic values of Bookerâ€“Browning); the two
  inequalities the theorem needs, fâ‚‚(9) â‰¥ 0.9995 and the Fâ‚‚-integral â‰¤ 2.625, are certified in interval arithmetic
  (`certify_dhr.py`: fâ‚‚(9) âˆˆ [0.999711, 0.999738], integral â‰¤ 2.623903, from the Bookerâ€“Browning enclosures of
  Î±â‚‚, Î²â‚‚ and the monotonicity of Fâ‚‚, fâ‚‚), so the theorem is unconditional. (7,4) is missed by 0.03; at level 2/3 the
  pair (6,4) gets margin 0.196, at 19/20 (5,3) gets 0.24.
* Section 8: the remainder beyond level 1/2 is Î£_{de â‰¤ N^Î¸} Î»_d Î»â€²_e R(d,e;N), i.e. the distribution of the products
  dÎ± (weighted by Î»_d r(Î±)) in the class N mod e with the extra weight r((N âˆ’ dÎ±)/e): a dispersion problem of
  Bombieriâ€“Friedlanderâ€“Iwaniec / Titchmarsh-divisor type. `modulus_sum_scan.py` measures the sums of a_{G_{d,e}}(N)
  over all coprime odd square-free (d,e) with de â‰¤ N^Î¸ at fixed N â‰ˆ 2Â·10â¶: the signed sums (all +, Î¼(de), random
  signs) are all of size (Î£ aÂ²)^{1/2}, i.e. square-root cancellation, no sign correlation.
* Numerics: the Rankinâ€“Selberg estimate of â€–Gâ€–Â² from the data satisfies â€–G_{m,mâ€²}â€–Â² â‰ˆ 0.37 (log mmâ€²)Â² for
  10Â² â‰¤ mmâ€² â‰¤ 1.6Â·10â´; an exact check at level 20 (G_{1,5} = (8/3) f_20, â€–f_20â€–Â² = 0.08135 from the period lattice
  of 20a1 and from Rankinâ€“Selberg) validates the whole numerical chain. `local_densities.py` makes the Siegel main
  term exact for every n (square-free m, mâ€²), checked against the levels 4, 12 (no cusp forms: T = M exactly) and
  20 (T âˆ’ M = (8/3) a_f(n) for all n); with it the unbiased estimate for (1,5) is 0.5785, the exact value.

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
            dhr_sieve.py              DHR sieve functions F_2, f_2 (alpha_2, beta_2 by the adjoint method) and the weighted sieve (Table 4)
            certify_dhr.py            interval-arithmetic enclosures of F_2, f_2 and the certification of Lemma 7.2 (Theorem 7.1)
            modulus_sum_scan.py       Table 5: sums of a_{G_{d,e}}(N) over the moduli at fixed N (cancellation test)
            test_rs_formula.py        numerical test of Theorem 6.7 (newform components of G_{1,p} via T_3, norms, central values)
            plot_cusp_decomposition.py  Figure 3: the norm cusp by cusp (from results/cusp_lower_bound.txt)
            plot_dhr_sieve.py         Figure 4: F_2, f_2 against the vector sieve; thresholds against theta (cache results/dhr_thresholds.txt)
            plot_modulus_sums.py      Figure 5: cancellation in the sum over the moduli, N = 3000002
figures/    mean_square, sieve_thresholds, cusp_decomposition, dhr_sieve, modulus_sums  (.pdf and .png)
results/    mean_square_check_x21.txt              Table 1
            cusp_lower_bound.txt, cusp_lower_bound_large.txt   Table 3 (X = 2^21; X = 2^24 for (1,15015), (1009,1013))
            cusp_lower_bound_1p.txt                (1,p) for p up to 10^6 (Section 6.6(f))
            dhr_sieve.txt                          Table 4 and the parameters of Theorem 7.1
            certify_dhr_output.txt                 certified enclosures (Lemma 7.2)
            modulus_sum_scan.txt                   Table 5
            dhr_thresholds.txt                     data of Figure 4(b)
            norm_level20.txt                       the exact check at level 20
            certify_modular_output_twothirds.txt   certified parameters at level 2/3 (bounds 6 and 4)
            certify_modular_output_nineteen.txt    certified parameters at level 19/20 (Corollary 1.4: bounds 5 and 4)
            sieve_thresholds.txt                   data of Figure 2
```

## Reproducing the numbers

Requirements: Python â‰¥ 3.10 with `numpy`, `scipy`, `mpmath` (versions used: `requirements.txt`).

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
python dhr_sieve.py 0.5 2/3 0.95       # F_2, f_2, alpha_2, beta_2 and Table 4 (~3 min)
python certify_dhr.py 20000            # Lemma 7.2: certified f_2(9) and the F_2-integral (directed rounding, ~1 min)
python modulus_sum_scan.py             # Table 5: four values of N ~ 2e6 (~1 min)
python plot_cusp_decomposition.py      # Figure 3 (seconds)
python plot_dhr_sieve.py               # Figure 4 (~5 min the first time; --recompute redoes the threshold scan)
python plot_modulus_sums.py            # Figure 5 (~20 s)
```

What is rigorous: `certify_modular.py` (interval arithmetic and directed rounding) enters the proofs of
Corollary 1.4 (and, at the level 1/2, Theorem 1.1 of Part I), and `certify_dhr.py` enters the proof of Theorem 7.1;
`mean_square_check.py`, `norm_level20.py`, `cusp_lower_bound.py`, `dhr_sieve.py`, `modulus_sum_scan.py` and the
plots are floating point and are evidence only.

## License

Code: MIT (see `LICENSE`). Manuscript (`paper/`): CC BY 4.0.
