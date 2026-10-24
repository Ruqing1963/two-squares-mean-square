# Sums of two squares with few prime factors, II: the cuspidal error in mean square and almost all N

Code, data and manuscript for

> R. Chen, *Sums of two squares with few prime factors, II: the cuspidal error in mean square and almost all N*, preprint, 2026.

Part I is [two-squares-theta-level](https://github.com/Ruqing1963/two-squares-theta-level)
(DOI [10.5281/zenodo.23618710](https://doi.org/10.5281/zenodo.23618710)).

## What the paper does

Let m, m′ be coprime odd integers, Θ_{m,m′}(z) = θ²(mz)θ²(m′z), G = G_{m,m′} its cuspidal part and L = 4mm′, so that
a_G(N) = T(m,m′;N) − M(m,m′;N) is the error term in the number of representations of N by m(x²+y²)+m′(z²+w²).

* **Theorem A (mean square).** For every G ∈ S₂(Γ₀(L)) and Y ≥ 1,
  Σ_{n ≤ Y} |a_G(n)|² ≤ e^{8π} Y (3 + 17Y/L) ‖G‖². With Waibel's norm bound this gives
  Σ_{N ≤ Y} |T − M|² ≪ Y² (mm′)^{−1/2+ε} for Y ≥ mm′: on average over N the error is N^{1/2}(mm′)^{−1/4},
  which is what the experiments of Part I show. The proof is a covering argument (the function y²|G|² is
  Γ₀(L)-invariant and the rectangle [0,1] × [1/Y, 2/Y] covers the fundamental domain at most 3 + 17Y/L times).
* **Theorem B.** The weighted sequence r(n) r(N − n) has level of distribution 2/3 for all but O(X^{1−δ})
  integers N ≡ 2 (mod 4) up to X (Part I proves 2/5 for every N; Conjecture 5.1 of Part I would give 1/2).
* **Corollary C.** For these N there are ≥ 0.13 𝔖(N) N/(log N)² integers n < N with n and N − n sums of two
  squares, Ω(n(N − n)) ≤ 6 and max(Ω(n), Ω(N − n)) ≤ 4 (Part I: 9 and 6 for every N). The sieve parameters
  θ = 2/3, α = 3/50, λ = 1/3, β₁ = 51/100 are certified in interval arithmetic.
* The uniform conjecture of Part I is rewritten as |a_G(N)|² = 4πN Δ_L(N,N) ‖G‖² cos²∠(G, P_N) with P_N the
  N-th Poincaré series: Theorem A says the average of cos² over N is O(1/L), the random-vector level, and the
  conjecture asks that no single N exceed the average by more than L^{1/2}.
* Numerically, ‖G_{m,m′}‖² ≈ 3 (mm′)^{1/4} for 10² ≤ mm′ ≤ 1.6·10⁴ (Rankin–Selberg estimate from the data),
  smaller than the general bound (mm′)^{1/2+ε}; a norm bound ‖G‖ ≪ (mm′)^{1/8+ε} would give the levels 4/9
  (every N) and 4/5 (almost all N).

## Repository layout

```
paper/      two_squares_mean_square.tex, two_squares_mean_square.pdf
code/       mean_square_check.py      mean square of T - M over N <= Y and the estimate of ||G||^2 (Table 1)
            certify_modular.py        Part I's interval-arithmetic certification, extended by the argument 2/3
            sup_error_scan.py         Part I: T(m,m';N) by FFT and the exact Siegel main term (used by mean_square_check.py)
            modular_sieve_scan.py     Part I: floating-point parameter scan (uses sieve_core.py)
            sieve_core.py             linear sieve functions
results/    mean_square_check_x21.txt              Table 1
            certify_modular_output_twothirds.txt   the certified parameters of Corollary C
```

## Reproducing the numbers

Requirements: Python ≥ 3.10 with `numpy`, `scipy`, `mpmath` (versions used: `requirements.txt`).

```
cd code
python certify_modular.py 2/3          # Corollary C: level 2/3, two parameter sets (certified, directed rounding)
python mean_square_check.py            # Table 1: X = 2^21, ten pairs (m, m') (~1 min, floating point)
python modular_sieve_scan.py 0.666667  # floating-point scan at level 2/3 (not certified)
python modular_sieve_scan.py 0.444444  # the same at level 4/9 (Proposition 5.3)
```

What is rigorous: `certify_modular.py` (interval arithmetic and directed rounding) enters the proof of
Corollary C; `mean_square_check.py` is floating point and is evidence only (Section 5 of the paper).

## License

Code: MIT (see `LICENSE`). Manuscript (`paper/`): CC BY 4.0.
