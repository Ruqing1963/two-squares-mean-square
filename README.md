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
* **Theorem 1.5 (lower bound for the norm, Section 6).** (i) ‖G_{m,m′}‖² ≥ (1/4π) Σ_{Q ∈ {m,4m,m′,4m′}} Σ_n |a_{G_{1,mm′}}(n)|² n⁻¹ e^{−4πn/Q}
  (the regions above the four Atkin–Lehner cusps at which G_{m,m′} becomes G_{1,mm′}); (ii) ‖G_{m,m′}‖ = ‖G_{1,mm′}‖:
  the norm depends only on mm′ (W_m is an isometry); (iii) ‖G_{m,m′}‖² ≥ (1/π)(log mm′)² − O(log(2mm′)(log log 3mm′)^C):
  the exponent ε in the upper bound (mm′)^ε of Part I cannot be removed. Numerically (Table 3, `cusp_lower_bound.py`)
  the regions above all cusps carry 93–98 % of the norm; the cusp sums Λ(m,m′) of the (m,m′) fundamental domain are
  ≈ (1/π)((log m)² + (log m′)²) + (5/π²) log m log m′ + (terms for the other divisors): they depend on the factorisation,
  the norm does not (checked: (3,5)/(1,15), (3,7)/(1,21), (5,7)/(1,35), (13,29)/(1,377) agree to 0.3 %,
  `results/norm_split_check.txt`). The "0.37 (log mm′)²" of Section 5 is (1/π)(log mm′)² plus lower-order terms.
* **Theorem 6.5 (the exact order).** ‖G_{m,m′}‖² ≪ (log 3mm′)² (log log 3mm′)^C, and ≪_k (log 3mm′)² when
  ω(mm′) ≤ k; with Theorem 1.5, ‖G_{m,m′}‖² ≍ (log mm′)² for mm′ with boundedly many prime factors, and the
  uniform bound T − M ≪ d(N) N^{1/2} (log 3mm′)(log log 3mm′)^C replaces the (mm′)^ε of Part I. Proof: the two
  ε-losses of Part I (d(α)d(β) and Waibel's genus bound) are replaced by the mean value of r(α)² and by the local
  densities. Corollary 6.6: the spectral second moment Σ_φ |⟨Θ,φ⟩|²/‖φ‖² is ≍ (log mm′)². For (1,p), p ≤ 10⁶,
  Λ(1,p) = 0.319 (log p)² + 0.42 log p − 1.5 (fit), against the leading term 1/π = 0.318.
* **Theorem 6.7 (Rankin–Selberg formula).** For every newform φ of weight 2 and level 4mm′ (mm′ odd square-free),
  ⟨Θ_{m,m′}, φ⟩ = (4/π²) χ₋₄(m′) η_{m′}(φ) L(1,φ) L(1,φ⊗χ₋₄) — all local factors cancel; and for every newform f of
  level 2^e M₁ | 4n and every d | 4n/M, ⟨Θ_{1,n}, f(dz)⟩ = (4/π²) χ₋₄(M₁) η_{M₁}(f) ∏_{p | n/M₁}(1+χ₋₄(p)) L(1,f) L(1,f⊗χ₋₄)/d
  (zero unless every prime of n/M₁ is 1 mod 4; the a_f(p)-dependence of the three local factors cancels). Proof:
  Hecke's weight-1 Eisenstein series θ² unfolded over the cusps of Γ₀(4n) above ∞, the Atkin–Lehner cusps W_Q∞,
  plus the action of W_Q on oldforms (Lemma 6.8). Verified numerically for newforms at levels 20, 52, 116, 148
  (`test_rs_formula.py`: ratios 1.00000; 0.99975; 0.99997, 0.99878, 0.99760; 1.00054, 0.99932) and for oldforms of
  levels 11, 26, 58, 74 inside levels 44, 52, 116, 148 (`test_rs_old.py`: 1.00071, 1.00184, 0.99873, 1.00068, 1.00035).
  Corollary 6.10: ‖G_{1,n}‖² = (16/π⁴) Σ_f γ_f 4^{ω(n/M₁)} L(1,f)² L(1,f⊗χ)² over the newforms f of all levels dividing
  4n (γ_f from the Gram matrix of {f(dz)}), so the norm bounds are bounds for this second moment; the newform part
  Σ L(1,φ)²L(1,φ⊗χ)²/⟨φ,φ⟩ ≪ (log 3n)²(log log 3n)^C.
* **Theorem 7.1 (two-dimensional sieve, Section 7).** Replacing the vector sieve of Part I by the Diamond–Halberstam–Richert
  sieve of dimension 2 (sieving limit β₂ = 4.2665 instead of 2 + 2√e = 5.30) gives, from the same level 1/2,
  Ω(n(N − n)) ≤ 7 and max(Ω(n), Ω(N − n)) ≤ 5 for ≥ 0.25 Σ_r(N) N/(log N)² integers n, for every large N ≡ 2 (mod 4):
  (7,5) instead of (8,5). `dhr_sieve.py` computes F₂, f₂ from the differential-delay system and determines
  α₂ = 5.357727, β₂ = 4.266450 (agreeing to 1e-8 with the ball-arithmetic values of Booker–Browning); the two
  inequalities the theorem needs, f₂(9) ≥ 0.9995 and the F₂-integral ≤ 2.625, are certified in interval arithmetic
  (`certify_dhr.py`: f₂(9) ∈ [0.999711, 0.999738], integral ≤ 2.623903, from the Booker–Browning enclosures of
  α₂, β₂ and the monotonicity of F₂, f₂), so the theorem is unconditional. (7,4) is missed by 0.03; at level 2/3 the
  pair (6,4) gets margin 0.196, at 19/20 (5,3) gets 0.24.
* Section 8: the remainder beyond level 1/2 is Σ_{de ≤ N^θ} λ_d λ′_e R(d,e;N), i.e. the distribution of the products
  dα (weighted by λ_d r(α)) in the class N mod e with the extra weight r((N − dα)/e): a dispersion problem of
  Bombieri–Friedlander–Iwaniec / Titchmarsh-divisor type. `modulus_sum_scan.py` measures the sums of a_{G_{d,e}}(N)
  over all coprime odd square-free (d,e) with de ≤ N^θ at fixed N ≈ 2·10⁶: the signed sums (all +, μ(de), random
  signs) are all of size (Σ a²)^{1/2}, i.e. square-root cancellation, no sign correlation.
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
            dhr_sieve.py              DHR sieve functions F_2, f_2 (alpha_2, beta_2 by the adjoint method) and the weighted sieve (Table 4)
            certify_dhr.py            interval-arithmetic enclosures of F_2, f_2 and the certification of Lemma 7.2 (Theorem 7.1)
            modulus_sum_scan.py       Table 5: sums of a_{G_{d,e}}(N) over the moduli at fixed N (cancellation test)
            test_rs_formula.py        numerical test of Theorem 6.7(i) (newform components of G_{1,p} via T_3, norms, central values)
            test_rs_old.py            numerical test of Theorem 6.7(ii) (oldform components)
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
