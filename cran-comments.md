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

* This is a new release.
* "Possibly mis-spelled words in DESCRIPTION" -> Double-checked, they are correctly spelled.
* "Found the following (possibly) invalid URLs" -> Double-checked, they are correctly spelled.
* Observe that the use of `set.seed` in `unif_stat_MC` and `int_sph_MC` serves the purpose of reproducibility and debugging when running simulations in parallel. Since the default is *not* to set any seeds (`seeds = NULL`), unless the user explicitly demands otherwise, the package is compliant with the policy "Packages should not modify the global environment (user’s workspace)."