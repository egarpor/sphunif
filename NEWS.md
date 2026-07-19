# sphunif 1.0.0

* Initial release.

# sphunif 1.0.1

* Switch to `doFuture` for parallel backend and `progressr` for progress monitoring in Monte Carlo simulations.

# sphunif 1.1.0

* Update references.
* Drop C++11 requirement to adhere to new CRAN policies.
* Drop `personList()` and `citEntry()`.
* Fix broken URLs.
* Update `comets` dataset.
* More unit tests.

# sphunif 1.3.0

* Add `"Poisson"`, `"Softmax"`, and `"Stereo"` tests.
* Add `"Sobolev"` tests.
* Vectorization of test-specific parameters in `unif_stat()`, `unif_test(p_value = "MC")`, and `unif_stat_MC()`.
* Use `doRNG::%dorng%` in `unif_stat_MC()` and `int_sph_MC()` to fix a bug when `seeds` was not set to `NULL`.

# sphunif 1.4.0

* New uniform spherical cap distribution in `unif_cap`.
* New non-uniform data generating-processes in `r_alt()`: `"MC"` and `"AUD"`.
* Add replication codes for several publications.
* Add support for `method` argument in `unif_stat_distr()` and `unif_test()`.
* Add warnings on absence of argument `n` in `unif_stat_distr()`.
* Rename argument `axial_MvMF` to `axial_mix` in `r_alt()`.
* Normalize Poisson kernel to improve numerical stability.
* Fix bugs in the coherence between asymptotic distributions and statistics of `"Poisson"` and `"Softmax"` tests.
* Add `"Stereo"` test null distribution.
* Fix bug when supplying `crit_val` to `unif_stat_MC()` in vectorized statistics.

# sphunif 1.4.1

* Fix bug on `g_i_k(k = 0)` for `p = 2` and make more robust the arguments for `g_i_k`.
* Add `"Stein"` test.
* Add `rot_ab` and `H_ab` convenience functions to rotate data.
* Fix NOTE on Rd cross-references.

# sphunif 1.4.2

* Update the reference Fernández-de-Marcos and García-Portugués (2024).

# sphunif 1.4.3

* Fix bug in `weights_dfs_Sobolev(type = "Watson")` (did not affect `unif_test(type = "Watson")`).
* Add further unit testing for `p_Sobolev()`.
* Update reference García-Portugués et al. (2025).
* Update vignette.
* Fix typos in documentation.

# sphunif 1.4.4

* Add `cmb` dataset

# sphunif 1.5.0

* Efficiency improvements in `unif_stat_MC()`, `unif_test()`, and
  `sph_stat_Sobolev()`.
* Fix the rejection direction in `unif_test(p_value = "crit_val")` (rejections
  were reported inverted).
* Fix `d_proj_unif_cap()`, whose density was identically zero due to a wrong
  support bound.
* Fix the Stephens (1970) modification in `p_cir_stat_Watson()` (it used the
  Kuiper constants, breaking coherence with `d_cir_stat_Watson()`).
* Fix the Stephens (1970) small-sample constant of the Kolmogorov-Smirnov
  statistic in `cir_stat_Kuiper(KS = TRUE)` (was `0.21 / n`, should be
  `0.11 / n`, per Stephens (1970) and Fernández-de-Marcos and
  García-Portugués (2024, Table 8)).
* Fix `unif_stat()`, `unif_test()`, and `unif_stat_distr()` to accept numeric
  `type` vectors of length greater than one.
* Fix `unif_stat_distr(approx = "MC")` with unsorted evaluation points `x`.
* Fix `unif_stat()` so that `Pycke` for `p >= 4` and vectorized columns for
  `type = "all"` are handled correctly.
* Fix `unif_test(p_value = "MC")` with a single significance level.
* `unif_stat_MC()` now simulates `ceiling(M / chunks)` replications per chunk.
* Forward `Stein_cf` in `d_Sobolev()`, `p_Sobolev()`, and `q_Sobolev()`.
* More unit tests.
