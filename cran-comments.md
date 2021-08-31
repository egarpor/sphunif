## Test environments

* local R installation, R 4.1.0
* ubuntu 16.04 (on travis-ci), R 3.5.0
* win-builder (release, devel)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
* Windows Server 2008 R2 SP1, R-patched, 32/64 bit
* macOS 10.13.6 High Sierra, R-release, brew
* macOS 10.13.6 High Sierra, R-release, CRAN's setup
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC
* Debian Linux, R-release, GCC
* Debian Linux, R-devel, GCC

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a resubmission where I have
    - added missing @return in r_alt;
    - removed the uses of ::: in documentation examples.
    - (I have also changed the parallel backend from doParallel to doFuture.)
* "Possibly mis-spelled words in DESCRIPTION" -> Double-checked, they are correctly spelled.
* "Found the following (possibly) invalid URLs" -> Double-checked, they are correct. They come from the usage of \doi{} with correct doi's. The redirections might be causing timeouts.
* Please observe that the use of `set.seed` in `unif_stat_MC` and `int_sph_MC` serves the very important purpose of ensuring full reproducibility and facilitating debugging when running simulations in parallel. Since the default is **not** to set any seeds (`seeds = NULL`), unless the user explicitly demands otherwise by overriding `seeds`, I believe the package is compliant with CRAN's policy.