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

# sphunif 1.3.1

* Normalize Poisson kernel to improve numerical stability.
* Add uniform spherical cap distribution in `unif_cap`.
* New non-uniform data generating-processes in `r_alt()`: `"MC"` and `"AUD"`.
* Rename argument `axial_MvMF` to `axial_mix` in `r_alt()`.
