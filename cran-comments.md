## Test environments

* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel, release)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (on R-hub)

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release.
* "Uses the superseded package: 'doSNOW'". Use of the 'doSNOW' package as opposed to the 'doParallel' package is required due to the support of the printed txtProgressBar in the 'doSNOW' package.