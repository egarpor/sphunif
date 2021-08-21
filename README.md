# sphunif

[![License:
GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://travis-ci.com/egarpor/sphunif.svg?branch=master)](https://travis-ci.com/egarpor/sphunif)
[![](https://codecov.io/gh/egarpor/sphunif/branch/master/graph/badge.svg)](https://codecov.io/gh/egarpor/sphunif)
[![](https://www.r-pkg.org/badges/version/sphunif?color=green)](https://cran.r-project.org/package=sphunif)
[![](http://cranlogs.r-pkg.org/badges/grand-total/sphunif?color=green)](https://cran.r-project.org/package=sphunif)
[![](http://cranlogs.r-pkg.org/badges/last-month/sphunif?color=green)](https://cran.r-project.org/package=sphunif)

<!-- <img src="" alt="sphunif  hexlogo" align="right" width="200" style="padding: 0 15px; float: right;"/> -->

## Overview

Implementation of circa 40 tests of uniformity on the circle, sphere,
and hypersphere. Software companion for the (evolving) review “*An
overview of uniformity tests on the hypersphere*” (García-Portugués and
Verdebout, 2018) and the paper “*On a projection-based class of
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
the `type` of test to be performed. For example, the `"Kuiper"` test:

``` r
library(sphunif)
unif_test(data = cir_data, type = "Kuiper", verbose = FALSE) # An htest object
#> 
#>  Kuiper test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.89659, p-value = 0.8452
#> alternative hypothesis: any alternative to circular uniformity
```

By default, the calibration of the test statistic is done by Monte
Carlo. This can be changed with `p_value = "asymp"` to employ asymptotic
distributions (faster, but not available for all tests):

``` r
unif_test(data = cir_data, type = "Kuiper", p_value = "MC",
          verbose = FALSE) # Monte Carlo
#> 
#>  Kuiper test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.89659, p-value = 0.8469
#> alternative hypothesis: any alternative to circular uniformity
unif_test(data = cir_data, type = "Kuiper", p_value = "asymp") # Asymp. distr.
#> 
#>  Kuiper test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.89659, p-value = 0.8452
#> alternative hypothesis: any alternative to circular uniformity
```

You can perform *several* tests with a *single* call to `unif_test`.
Choose the available circular uniformity tests from

``` r
avail_cir_tests
#>  [1] "Ajne"            "Bakshaev"        "Bingham"         "Cressie"        
#>  [5] "Cuesta_Albertos" "Feltz_Goldin"    "Gine_Fn"         "Gine_Gn"        
#>  [9] "Gini"            "Gini_squared"    "Greenwood"       "Hermans_Rasson" 
#> [13] "Hodges_Ajne"     "Kuiper"          "Log_gaps"        "Max_uncover"    
#> [17] "Num_uncover"     "PAD"             "PCvM"            "PRt"            
#> [21] "Pycke"           "Pycke_q"         "Range"           "Rao"            
#> [25] "Rayleigh"        "Riesz"           "Rothman"         "Vacancy"        
#> [29] "Watson"          "Watson_1976"
```

For example:

``` r
unif_test(data = cir_data, type = c("Kuiper", "Watson", "Ajne"),
          verbose = FALSE) # A *list* of htest objects
#> $Kuiper
#> 
#>  Kuiper test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.89659, p-value = 0.8452
#> alternative hypothesis: any alternative to circular uniformity
#> 
#> 
#> $Watson
#> 
#>  Watson test of circular uniformity
#> 
#> data:  cir_data
#> statistic = 0.036003, p-value = 0.8694
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
#>  [1] "Ajne"            "Bakshaev"        "Bingham"         "Cai"            
#>  [5] "Cuesta_Albertos" "Gine_Fn"         "Gine_Gn"         "PAD"            
#>  [9] "PCvM"            "PRt"             "Pycke"           "Rayleigh"       
#> [13] "Rayleigh_HD"     "Riesz"
```

The default `type = "all"` equals `type = avail_sph_tests`:

``` r
unif_test(data = sph_data, type = "all", p_value = "MC", verbose = FALSE)
#> $Ajne
#> 
#>  Ajne test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 0.079876, p-value = 0.9576
#> alternative hypothesis: any non-axial alternative to spherical uniformity
#> 
#> 
#> $Bakshaev
#> 
#>  Bakshaev test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.2727, p-value = 0.4388
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
#> $Cai
#> 
#>  Cai test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 27.401, p-value = 0.2655
#> alternative hypothesis: unclear consistency
#> 
#> 
#> $Cuesta_Albertos
#> 
#>  Cuesta-Albertos et al. (2009) test of spherical uniformity with k = 50
#> 
#> data:  sph_data
#> statistic = 1.4619, p-value = 0.2688
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $Gine_Fn
#> 
#>  Gine's Fn test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.8889, p-value = 0.2229
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $Gine_Gn
#> 
#>  Gine's Gn test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.5694, p-value = 4e-04
#> alternative hypothesis: any axial alternative to spherical uniformity
#> 
#> 
#> $PAD
#> 
#>  Projected Anderson-Darling test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.1068, p-value = 0.3196
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PCvM
#> 
#>  Projected Cramer-von Mises test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 0.15909, p-value = 0.4388
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PRt
#> 
#>  Projected Rothman test of spherical uniformity with t = 0.333
#> 
#> data:  sph_data
#> statistic = 0.19003, p-value = 0.5131
#> alternative hypothesis: any alternative to spherical uniformity if t is irrational
#> 
#> 
#> $Pycke
#> 
#>  Pycke test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 0.091169, p-value = 0.1853
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $Rayleigh
#> 
#>  Rayleigh test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 0.52692, p-value = 0.9184
#> alternative hypothesis: mean direction different from zero
#> 
#> 
#> $Rayleigh_HD
#> 
#>  HD-standardized Rayleigh test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = -1.0096, p-value = 0.9184
#> alternative hypothesis: mean direction different from zero
#> 
#> 
#> $Riesz
#> 
#>  Riesz test of spherical uniformity
#> 
#> data:  sph_data
#> statistic = 1.2727, p-value = 0.4388
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

## Data applications in astronomy

The data applications in García-Portugués, Navarro-Esteban and
Cuesta-Albertos (2020) can be reproduced through the script
[data-application-ecdf.R](https://github.com/egarpor/sphunif/blob/master/application/data-application-ecdf.R).
The code snippet below illustrates the testing of the uniformity of
orbits of planets and long-period comets.

``` r
## Planets

# Load data
data("planets")

# Add normal vectors
planets$normal <- cbind(sin(planets$i) * sin(planets$om),
                       -sin(planets$i) * cos(planets$om),
                       cos(planets$i))

# Tests without Pluto
unif_test(data = planets$normal[-9, ], type = c("PCvM", "PAD", "PRt"),
          p_value = "MC")
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> $PCvM
#> 
#>  Projected Cramer-von Mises test of spherical uniformity
#> 
#> data:  planets$normal[-9, ]
#> statistic = 1.2894, p-value < 2.2e-16
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PAD
#> 
#>  Projected Anderson-Darling test of spherical uniformity
#> 
#> data:  planets$normal[-9, ]
#> statistic = 7.6553, p-value < 2.2e-16
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PRt
#> 
#>  Projected Rothman test of spherical uniformity with t = 0.333
#> 
#> data:  planets$normal[-9, ]
#> statistic = 1.7251, p-value < 2.2e-16
#> alternative hypothesis: any alternative to spherical uniformity if t is irrational

## Comets

# Load data
data("comets")

# Add normal vectors
comets$normal <- cbind(sin(comets$i) * sin(comets$om),
                       -sin(comets$i) * cos(comets$om),
                       cos(comets$i))

# Exclude the C/1882 R1-X (Great September comet) records with X = B, C, D
comets_ccf2009 <- comets[comets$ccf2009, ][-c(13:15), ]

# Tests for the data in Cuesta-Albertos et al. (2009)
tests_ccf2009 <- unif_test(data = comets_ccf2009$normal,
                           type = c("PCvM", "PAD", "PRt"),
                           p_value = "asymp", K_max = 5e4)
tests_ccf2009
#> $PCvM
#> 
#>  Projected Cramer-von Mises test of spherical uniformity
#> 
#> data:  comets_ccf2009$normal
#> statistic = 0.27607, p-value = 0.1011
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PAD
#> 
#>  Projected Anderson-Darling test of spherical uniformity
#> 
#> data:  comets_ccf2009$normal
#> statistic = 1.6733, p-value = 0.07445
#> alternative hypothesis: any alternative to spherical uniformity
#> 
#> 
#> $PRt
#> 
#>  Projected Rothman test of spherical uniformity with t = 0.333
#> 
#> data:  comets_ccf2009$normal
#> statistic = 0.36307, p-value = 0.1207
#> alternative hypothesis: any alternative to spherical uniformity if t is irrational
```

## References

García-Portugués, E., Navarro-Esteban, P., and Cuesta-Albertos, J. A.
(2020). On a projection-based class of uniformity tests on the
hypersphere. *arXiv:2008.09897*. <https://arxiv.org/abs/2008.09897>

García-Portugués, E. and Verdebout, T. (2018). An overview of uniformity
tests on the hypersphere. *arXiv:1804.00286*.
<https://arxiv.org/abs/1804.00286>
