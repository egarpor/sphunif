# sphunif

[![License:
GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R build
status](https://github.com/egarpor/sphunif/workflows/R-CMD-check/badge.svg)](https://github.com/egarpor/sphunif/actions)
[![R build
status](https://github.com/egarpor/sphunif/workflows/test-coverage/badge.svg)](https://github.com/egarpor/sphunif/actions)
[![](https://codecov.io/gh/egarpor/sphunif/branch/master/graph/badge.svg)](https://app.codecov.io/gh/egarpor/sphunif)
[![](https://www.r-pkg.org/badges/version/sphunif?color=green)](https://cran.r-project.org/package=sphunif)
[![](http://cranlogs.r-pkg.org/badges/grand-total/sphunif)](https://cran.r-project.org/package=sphunif)
[![](http://cranlogs.r-pkg.org/badges/last-month/sphunif)](https://cran.r-project.org/package=sphunif)

<!-- <img src="" alt="sphunif hexlogo" align="right" width="200" style="padding: 0 15px; float: right;"/> -->

## Overview

Implementation of more than 35 tests of uniformity on the circle,
sphere, and hypersphere. Software companion for the (evolving) review
*An overview of uniformity tests on the hypersphere* (García-Portugués
and Verdebout, 2018) and the paper *On a projection-based class of
uniformity tests on the hypersphere* (García-Portugués, Navarro-Esteban,
and Cuesta-Albertos, 2023). The package also provides several novel
datasets and gives the replicability for the data
applications/simulations in [different
publications](https://github.com/egarpor/sphunif?tab=readme-ov-file#replicability).

## Installation

Get the latest version from GitHub:

``` r
# Install the package and the vignettes
library(devtools)
install_github("egarpor/sphunif", build_vignettes = TRUE)

# Load package
library(sphunif)

# See main vignette
vignette("sphunif")
```

## Usage

### Circular data

You want to test if a sample of *circular* data is uniformly
distributed. For example, the following circular uniform sample in
*radians*:

``` r
set.seed(987202226)
cir_data <- runif(n = 10, min = 0, max = 2 * pi)
```

Call the main function in the `sphunif` package, `unif_test`, specifying
the `type` of test to be performed. For example, the `"Watson"` test:

``` r
library(sphunif)
unif_test(data = cir_data, type = "Watson") # An htest object
#> 
#>  Watson test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.036003, p-value = 0.8694
#> alternative hypothesis: any alternative to circular uniformity
```

By default, the calibration of the test statistic is done by Monte
Carlo. This can be changed with `p_value = "asymp"` to employ asymptotic
distributions (faster, but not available for all tests):

``` r
unif_test(data = cir_data, type = "Watson", p_value = "MC") # Monte Carlo
#> 
#>  Watson test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.036003, p-value = 0.8831
#> alternative hypothesis: any alternative to circular uniformity
unif_test(data = cir_data, type = "Watson", p_value = "asymp") # Asymp. distr.
#> 
#>  Watson test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.036003, p-value = 0.8694
#> alternative hypothesis: any alternative to circular uniformity
```

You can perform *several* tests with a *single* call to `unif_test`.
Choose the available circular uniformity tests from

``` r
avail_cir_tests
#>  [1] "Ajne"           "Bakshaev"       "Bingham"        "Cressie"       
#>  [5] "CCF09"          "FG01"           "Gine_Fn"        "Gine_Gn"       
#>  [9] "Gini"           "Gini_squared"   "Greenwood"      "Hermans_Rasson"
#> [13] "Hodges_Ajne"    "Kuiper"         "Log_gaps"       "Max_uncover"   
#> [17] "Num_uncover"    "PAD"            "PCvM"           "Poisson"       
#> [21] "PRt"            "Pycke"          "Pycke_q"        "Range"         
#> [25] "Rao"            "Rayleigh"       "Riesz"          "Rothman"       
#> [29] "Sobolev"        "Softmax"        "Stein"          "Vacancy"       
#> [33] "Watson"         "Watson_1976"
```

For example:

``` r
# A *list* of htest objects
unif_test(data = cir_data, type = c("Watson", "PAD", "Ajne"))
#> $Watson
#> 
#>  Watson test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.036003, p-value = 0.8694
#> alternative hypothesis: any alternative to circular uniformity
#> 
#> 
#> $PAD
#> 
#>  Projected Anderson-Darling test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.47247, p-value = 0.9051
#> alternative hypothesis: any alternative to circular uniformity
#> 
#> 
#> $Ajne
#> 
#>  Ajne test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.11016, p-value = 0.7361
#> alternative hypothesis: any non-axial alternative to circular uniformity
```

### Spherical data

You want to test if a sample of *spherical* data is uniformly
distributed. Consider a *non*-uniformly-generated sample of directions
in **Cartesian coordinates**:

``` r
# Sample data on S^2
set.seed(987202226)
theta <- runif(n = 50, min = 0, max = 2 * pi)
phi <- runif(n = 50, min = 0, max = pi)
sph_data <- cbind(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi))
```

The available spherical uniformity tests:

``` r
avail_sph_tests
#>  [1] "Ajne"        "Bakshaev"    "Bingham"     "CCF09"       "CJ12"       
#>  [6] "Gine_Fn"     "Gine_Gn"     "PAD"         "PCvM"        "Poisson"    
#> [11] "PRt"         "Pycke"       "Sobolev"     "Softmax"     "Stein"      
#> [16] "Stereo"      "Rayleigh"    "Rayleigh_HD" "Riesz"
```

The default `type = "all"` equals `type = avail_sph_tests`:

``` r
head(unif_test(data = sph_data, type = "all", p_value = "MC"))
#> $Ajne
#> 
#>  Ajne test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 0.079876, p-value = 0.9578
#> alternative hypothesis: any non-axial alternative to spherical uniformity
#> 
#> 
#> $Bakshaev
#> 
#>  Bakshaev (2010) test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.2727, p-value = 0.4418
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $Bingham
#> 
#>  Bingham test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 22.455, p-value < 2.2e-16
#> alternative hypothesis: scatter matrix different from constant
#> 
#> 
#> $CCF09
#> 
#>  Cuesta-Albertos et al. (2009) test of spherical uniformity with k = 50
#> 
#> data:  sph_data
#> statistic = 1.4619, p-value = 0.2775
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $CJ12
#> 
#>  Cai and Jiang (2012) test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 27.401, p-value = 0.262
#> alternative hypothesis: unclear consistency
#> 
#> 
#> $Gine_Fn
#> 
#>  Gine's Fn test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.8889, p-value = 0.2216
#> alternative hypothesis: any alternative to spherical uniformity
unif_test(data = sph_data, type = "Rayleigh", p_value = "asymp")
#> 
#>  Rayleigh test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 0.52692, p-value = 0.9129
#> alternative hypothesis: mean direction different from zero
```

### Higher dimensions

The *hyperspherical* setting is treated analogously to the spherical
setting, and the available tests are exactly the same
(`avail_sph_tests`). An example of testing uniformity with a sample of
size `100` on the $`9`$-sphere:

``` r
# Sample data on S^9
set.seed(987202226)
hyp_data <- r_unif_sph(n = 50, p = 10)

# Test
unif_test(data = hyp_data, type = "Rayleigh", p_value = "asymp")
#> 
#>  Rayleigh test of spherical uniformity
#> 
#> data:  hyp_data
#> statistic = 11.784, p-value = 0.2997
#> alternative hypothesis: mean direction different from zero
```

## Replicability

### *On a projection-based class of uniformity tests on the hypersphere*

The data application in García-Portugués, Navarro-Esteban, and
Cuesta-Albertos (2023) regarding the testing of uniformity of locations
for craters in Rhea can be reproduced through the script
[data-application-ecdf.R](https://github.com/egarpor/sphunif/blob/master/applications/data-application-ecdf.R).
The code snippet below is a simplified version.

``` r
# Load data
data(rhea)

# Add Cartesian coordinates
rhea$X <- cbind(cos(rhea$theta) * cos(rhea$phi),
                sin(rhea$theta) * cos(rhea$phi),
                sin(rhea$phi))

# Distribution of diameter
quantile(rhea$diameter)
#>       0%      25%      50%      75%     100% 
#>  10.0000  13.2475  17.0500  24.5600 449.8200

# Subsets of craters, according to diameter
ind_15_20 <- rhea$diameter > 15 & rhea$diameter < 20
ind_20 <- rhea$diameter > 20

# Sample sizes
nrow(rhea)
#> [1] 3596
sum(ind_15_20)
#> [1] 867
sum(ind_20)
#> [1] 1373

# Tests to be performed
type_tests <- c("PCvM", "PAD", "PRt")

# Tests
tests_rhea_15_20 <- unif_test(data = rhea$X[ind_15_20, ], type = type_tests,
                              p_value = "asymp", K_max = 5e4)
tests_rhea_20 <- unif_test(data = rhea$X[ind_20, ], type = type_tests,
                           p_value = "asymp", K_max = 5e4)
tests_rhea_15_20
#> $PCvM
#> 
#>  Projected Cramer-von Mises test of spherical uniformity
#> 
#> data:  rhea$X[ind_15_20, ]
#> statistic = 0.26452, p-value = 0.1176
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PAD
#> 
#>  Projected Anderson-Darling test of spherical uniformity
#> 
#> data:  rhea$X[ind_15_20, ]
#> statistic = 1.6854, p-value = 0.07209
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PRt
#> 
#>  Projected Rothman test of spherical uniformity with t = 0.333
#> 
#> data:  rhea$X[ind_15_20, ]
#> statistic = 0.31352, p-value = 0.1856
#> alternative hypothesis: any alternative to spherical uniformity if t is irrational
tests_rhea_20
#> $PCvM
#> 
#>  Projected Cramer-von Mises test of spherical uniformity
#> 
#> data:  rhea$X[ind_20, ]
#> statistic = 3.6494, p-value = 2.202e-09
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PAD
#> 
#>  Projected Anderson-Darling test of spherical uniformity
#> 
#> data:  rhea$X[ind_20, ]
#> statistic = 18.482, p-value < 2.2e-16
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PRt
#> 
#>  Projected Rothman test of spherical uniformity with t = 0.333
#> 
#> data:  rhea$X[ind_20, ]
#> statistic = 5.3485, p-value = 3.481e-09
#> alternative hypothesis: any alternative to spherical uniformity if t is irrational
```

### *A Cramér–von Mises test of uniformity on the hypersphere*

The data application in García-Portugués, Navarro-Esteban, and
Cuesta-Albertos (2021) regarding the testing of the uniformity of
locations for craters in Venus can be reproduced through the script
[data-application-cvm.R](https://github.com/egarpor/sphunif/blob/master/applications/data-application-cvm.R).

### *A stereographic test of spherical uniformity*

The script
[simulations-stereo.R](https://github.com/egarpor/sphunif/blob/master/simulations/simulations-stereo.R)
contains some of the numerical experiments in Fernández-de-Marcos and
García-Portugués (2024).

### *On a class of Sobolev tests for symmetry of directions, their detection thresholds, and asymptotic powers*

The data application in García-Portugués, Paindaveine, and Verdebout
(2024) regarding the symmetry of comet orbits can be reproduced through
the script
[data-application-sobolev.R](https://github.com/egarpor/sphunif/blob/master/applications/data-application-sobolev.R).

## References

Fernández-de-Marcos, A. and García-Portugués, E. (2024). A stereographic
test of spherical uniformity. *Statistics and Probability Letters*,
215:110218.
[doi:10.1016/j.spl.2024.110218](https://doi.org/10.1016/j.spl.2024.110218).

García-Portugués, E., Navarro-Esteban, P., and Cuesta-Albertos, J. A.
(2023). On a projection-based class of uniformity tests on the
hypersphere. *Bernoulli*, 29(1):181–204.
[doi:10.1007/978-3-030-69944-4_12](https://doi.org/10.3150/21-BEJ1454).

García-Portugués, E., Navarro-Esteban, P., and Cuesta-Albertos, J. A.
(2021). A Cramér–von Mises test of uniformity on the hypersphere. In
Balzano, S., Porzio, G. C., Salvatore, R., Vistocco, D., and Vichi, M.
(Eds.), *Statistical Learning and Modeling in Data Analysis*, Studies in
Classification, Data Analysis and Knowledge Organization, pp. 107–116.
Springer, Cham.
[doi:10.1007/978-3-030-69944-4_12](https://doi.org/10.1007/978-3-030-69944-4_12).

García-Portugués, E., Paindaveine, D., and Verdebout, T. (2024). On a
class of Sobolev tests for symmetry of directions, their detection
thresholds, and asymptotic powers. *arXiv:2108.09874v2*.
[doi:10.48550/arXiv.2108.09874](https://doi.org/10.48550/arXiv.2108.09874).

García-Portugués, E. and Verdebout, T. (2018). An overview of uniformity
tests on the hypersphere. *arXiv:1804.00286*.
[doi:10.48550/arXiv.1804.00286](https://doi.org/10.48550/arXiv.1804.00286).
