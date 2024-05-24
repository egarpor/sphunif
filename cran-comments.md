## Test environments

* local R installation, R 4.2.2
* win-builder (release, devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

Added arXiv doi's in DESCRIPTION and documentation.

False positive NOTE in check_win_release() and check_win_devel() about the spelling of the words "al", "de", "et", and "ndez" in the DESCRIPTION file. These words are correctly spelled.

Other testing environments from devtools::check_rhub() and rhub::check_for_cran() are not available due to the following error message:
Error in curl::curl_fetch_memory(url, handle = handle) : 
  SSL peer certificate or SSH remote key was not OK: [builder.r-hub.io] SSL certificate problem: self signed certificate