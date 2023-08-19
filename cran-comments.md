## Test environments

* local R installation, R 4.2.2
* win-builder (release, devel)
* Windows Server 2022, R-release, 32/64 bit
* Windows Server 2022, R-devel, 64 bit
* Windows Server 2022, R-oldrel, 32/64 bit
* Windows Server 2022, R-patched, 32/64 bit
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC
* Debian Linux, R-release, GCC
* Debian Linux, R-devel, GCC

## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

* In this resubmission I **fix the previous WARNING**: sph_stat.cpp:357:7: note: cast one or both operands to int to silence this warning
* There is a **NOTE that is a false positive**:
  - Possibly misspelled words in DESCRIPTION:
      replicability (26:5)
      **This has been double-checked. 'replicability' is a word in common use in science.**
  - Found the following (possibly) invalid URLs:
      **All the URLs have been double-checked. They work fine.**