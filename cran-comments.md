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
* "Possibly mis-spelled words in DESCRIPTION" ->  Double-checked, they are correctly spelled.
* "Uses the superseded package: 'doSNOW'" -> Use of the 'doSNOW' package as opposed to the 'doParallel' package is required due to the support of the printed txtProgressBar in the 'doSNOW' package.
* "Found the following (possibly) invalid URLs" -> Double-checked, they are correctly spelled.
* "Found the following URLs which should use \doi (with the DOI name only): File 'rhea.Rd': https://dx.doi.org/10.1002/2015JE004940" -> But the URL can not be casted as \doi{10.1002/2015JE004940} since it is used inside an \href{}{} and an "! Undefined control sequence" will appear if doing so.
