# sphunif

[![License:
GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://app.travis-ci.com/egarpor/sphunif.svg?branch=master)](https://app.travis-ci.com/egarpor/sphunif)
[![](https://codecov.io/gh/egarpor/sphunif/branch/master/graph/badge.svg)](https://codecov.io/gh/egarpor/sphunif)
[![](https://www.r-pkg.org/badges/version/sphunif?color=green)](https://cran.r-project.org/package=sphunif)
[![](http://cranlogs.r-pkg.org/badges/grand-total/sphunif?color=green)](https://cran.r-project.org/package=sphunif)
[![](http://cranlogs.r-pkg.org/badges/last-month/sphunif?color=green)](https://cran.r-project.org/package=sphunif)

<!-- <img src="" alt="sphunif  hexlogo" align="right" width="200" style="padding: 0 15px; float: right;"/> -->

## Overview

Implementation of more than 30 tests of uniformity on the circle,
sphere, and hypersphere. Software companion for the (evolving) review
“*An overview of uniformity tests on the hypersphere*” (García-Portugués
and Verdebout, 2018) and the paper “*On a projection-based class of
uniformity tests on the hypersphere*” (García-Portugués, Navarro-Esteban
and Cuesta-Albertos, 2020).

## Installation

Get the latest version from GitHub:

``` r
# Install the package
library(devtools)
install_github("egarpor/sphunif")

# Load package
library(sphunif)
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
unif_test(data = cir_data, type = "Watson", verbose = FALSE) # An htest object
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
unif_test(data = cir_data, type = "Watson", p_value = "MC",
          verbose = FALSE) # Monte Carlo
#> 
#>  Watson test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.036003, p-value = 0.8817
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
#> [17] "Num_uncover"    "PAD"            "PCvM"           "PRt"           
#> [21] "Pycke"          "Pycke_q"        "Range"          "Rao"           
#> [25] "Rayleigh"       "Riesz"          "Rothman"        "Vacancy"       
#> [29] "Watson"         "Watson_1976"
```

For example:

``` r
unif_test(data = cir_data, type = c("Watson", "PAD", "Ajne"),
          verbose = FALSE) # A *list* of htest objects
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
#>  [1] "Ajne"        "Bakshaev"    "Bingham"     "CJ12"        "CCF09"      
#>  [6] "Gine_Fn"     "Gine_Gn"     "PAD"         "PCvM"        "PRt"        
#> [11] "Pycke"       "Rayleigh"    "Rayleigh_HD" "Riesz"
```

The default `type = "all"` equals `type = avail_sph_tests`:

``` r
head(unif_test(data = sph_data, type = "all", p_value = "MC", verbose = FALSE))
#> $Ajne
#> 
#>  Ajne test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 0.079876, p-value = 0.956
#> alternative hypothesis: any non-axial alternative to spherical uniformity
#> 
#> 
#> $Bakshaev
#> 
#>  Bakshaev (2010) test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.2727, p-value = 0.4419
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $Bingham
#> 
#>  Bingham test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 22.455, p-value = 3e-04
#> alternative hypothesis: scatter matrix different from constant
#> 
#> 
#> $CJ12
#> 
#>  Cai and Jiang (2012) test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 27.401, p-value = 0.2746
#> alternative hypothesis: unclear consistency
#> 
#> 
#> $CCF09
#> 
#>  Cuesta-Albertos et al. (2009) test of spherical uniformity with k = 50
#> 
#> data:  sph_data
#> statistic = 1.4619, p-value = 0.272
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $Gine_Fn
#> 
#>  Gine's Fn test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.8889, p-value = 0.2247
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
size `100` on the 9-sphere:

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

## Data application in astronomy

The data application in García-Portugués, Navarro-Esteban and
Cuesta-Albertos (2020) can be reproduced through the script
[data-application-ecdf.R](https://github.com/egarpor/sphunif/blob/master/application/data-application-ecdf.R).
The code snippet below illustrates the testing of the uniformity of
craters in Rhea.

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

## References

García-Portugués, E., Navarro-Esteban, P., and Cuesta-Albertos, J. A.
(2021). A Cramér–von Mises test of uniformity on the hypersphere. In
Balzano, S., Porzio, G. C., Salvatore, R., Vistocco, D., and Vichi, M.
(Eds.), *Statistical Learning and Modeling in Data Analysis*, Studies in
Classification, Data Analysis and Knowledge Organization, pp. 107–-116.
Springer, Cham.
[doi:10.1007/978-3-030-69944-4_12](https://doi.org/10.1007/978-3-030-69944-4_12).

García-Portugués, E., Navarro-Esteban, P., and Cuesta-Albertos, J. A.
(2020). On a projection-based class of uniformity tests on the
hypersphere. *arXiv:2008.09897*. <https://arxiv.org/abs/2008.09897>.

García-Portugués, E. and Verdebout, T. (2018). An overview of uniformity
tests on the hypersphere. *arXiv:1804.00286*.
<https://arxiv.org/abs/1804.00286>.

García-Portugués, E., Paindaveine, D., and Verdebout, T. (2021). On the
power of Sobolev tests for isotropy under local rotationally symmetric
alternatives. *arXiv:2108.XXXXX*. <https://arxiv.org/abs/2108.XXXXX>.
