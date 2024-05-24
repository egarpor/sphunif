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

False positive NOTE in check_win_release() and check_win_devel() about the spelling of the words "al", "de", "et", and "ndez" in the DESCRIPTION file. These words are correctly spelled.