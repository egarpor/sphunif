
# Clean workspace
rm(list = ls())

# Load packages
library(rotasym)
library(sphunif)

## Sunspots

# Load data used in García-Portugués, Paindaveine, and Verdebout (2020)
data("sunspots_births")

# Add Cartesian coordinates
sunspots_births$X <-
  cbind(cos(sunspots_births$phi) * cos(sunspots_births$theta),
        cos(sunspots_births$phi) * sin(sunspots_births$theta),
        sin(sunspots_births$phi))

# Data associated to the 23rd and 22nd solar cycles
sunspots_23 <- subset(sunspots_births, cycle == 23)
sunspots_22 <- subset(sunspots_births, cycle == 22)

# Obtain the signs with respect to the north pole c(0, 0, 1) -- these are the
# positions in S^1 of the longitudes associated to the data. The latitudes are
# disregarded to investigate if the longitudes are uniformly distributed, which
# is implied if rotational symmetry about (0, 0, 1) holds
signs_23 <- signs(sunspots_23$X, theta = c(0, 0, 1))
signs_22 <- signs(sunspots_22$X, theta = c(0, 0, 1))

# Sample sizes
nrow(signs_23)
nrow(signs_22)

# Tests to be performed
type_tests <- c("PCvM", "PAD", "PRt")

# Tests
tests_sunspots_23 <- unif_test(data = signs_23, type = type_tests,
                               p_value = "asymp", K_max = 5e4)
tests_sunspots_22 <- unif_test(data = signs_22, type = type_tests,
                               p_value = "asymp", K_max = 5e4)
tests_sunspots_23
tests_sunspots_22

## Planets

# Load data
data("planets")

# Add normal vectors
planets$normal <- cbind(sin(planets$i) * sin(planets$om),
                       -sin(planets$i) * cos(planets$om),
                       cos(planets$i))

# Tests to be performed
type_tests <- c("PCvM", "PAD", "PRt")

# Tests with Pluto
unif_test(data = planets$normal, type = type_tests, M = 1e6,
          p_value = "MC", cores = 4, chunks = 100, seed = 1:100)

# Tests without Pluto
unif_test(data = planets$normal[-9, ], type = type_tests, M = 1e6,
          p_value = "MC", cores = 4, chunks = 100, seed = 1:100)

## Comets

# Load data
data("comets")

# Add normal vectors
comets$normal <- cbind(sin(comets$i) * sin(comets$om),
                       -sin(comets$i) * cos(comets$om),
                       cos(comets$i))

# Tests to be performed
type_tests <- c("PCvM", "PAD", "PRt")

# Excluding the C/1882 R1-X (Great September comet) records with X = B, C, D
comets_ccf2009 <- comets[comets$ccf2009, ][-c(13:15), ]

# Sample size
nrow(comets_ccf2009)

# Tests for the data in Cuesta-Albertos et al. (2009)
tests_ccf2009 <- c(unif_test(data = comets_ccf2009$normal, type = type_tests,
                             p_value = "asymp", K_max = 5e4),
                   list("Cuesta_Albertos" =
                          unif_test(data = comets_ccf2009$normal,
                                    type = "Cuesta_Albertos", p_value = "MC",
                                    M = 1e4, cores = 4, chunks = 100,
                                    seed = 1:100))
                   )
tests_ccf2009

# Numbered comets
comets$numbered <- substr(comets$id, 1, 1) == "c"

# Filter comets as in Cuesta-Albertos et al. (2009)
comets_ellip <- subset(x = comets, subset = !(class %in% c("HYP", "PAR")) &
                         per_y >= 200 & !numbered)

# Remove "duplicated" records (as done in Cuesta-Albertos et al. (2009))
dup <- duplicated(round(comets_ellip$normal, 2))
comets_ellip <- comets_ellip[!dup, ]

# Sample size
nrow(comets_ellip)

# Tests for the extended dataset
tests_comets <- c(unif_test(data = comets_ellip$normal, type = type_tests,
                            p_value = "asymp", K_max = 5e4),
                   list("Cuesta_Albertos" =
                          unif_test(data = comets_ellip$normal,
                                    type = "Cuesta_Albertos", p_value = "MC",
                                    M = 1e4, cores = 4, chunks = 100,
                                    seed = 1:100))
)
tests_comets

## Craters named by the IUA

# Read data of named craters
data("craters")

# Add Cartesian coordinates
craters$X <- cbind(cos(craters$theta) * cos(craters$phi),
                   sin(craters$theta) * cos(craters$phi),
                   sin(craters$phi))

# Most crater-populated bodies and its type
num_craters <- sort(table(craters$target), decreasing = TRUE)
body_type <- craters$target_type[match(names(num_craters), craters$target)]
tab_craters <- data.frame("body" = names(num_craters), "type" = body_type,
                          "craters" = c(num_craters))
rownames(tab_craters) <- NULL
tab_craters

# Tests to be performed
type_tests <- c("PCvM", "PAD", "PRt")

# Select the bodies with more than 30 craters and that are either planets
# or moons
sel_bodies <- which(tab_craters$craters >= 30 & tab_craters$type != "Asteroid")
names_sel_bodies <- tab_craters$body[sel_bodies]
tab_craters[sel_bodies, ]

# How many filtered bodies?
length(sel_bodies)

# How many craters?
sum(tab_craters$craters[sel_bodies])

# Run the tests
tests_craters <- lapply(names_sel_bodies, function(body) {
  unif_test(data = craters$X[craters$target == body, , drop = FALSE],
            type = type_tests, p_value = "asymp", K_max = 5e4)
})
names(tests_craters) <- names_sel_bodies

# Matrix with the p-values for the bodies
bodies_pvalues <- t(sapply(tests_craters, function(body) {
  sapply(body, function(test) test$p.value)
}))
bodies_pvalues <- as.data.frame(cbind("n" = num_craters[sel_bodies],
                                      bodies_pvalues))
bodies_pvalues

## Rhea craters

# Load data
data(rhea)

# Add Cartesian coordinates
rhea$X <- cbind(cos(rhea$theta) * cos(rhea$phi),
                sin(rhea$theta) * cos(rhea$phi),
                sin(rhea$phi))

# Distribution of diameter
quantile(rhea$diameter)

# Subsets of craters, according to diameter
ind_15_20 <- rhea$diameter > 15 & rhea$diameter < 20
ind_20 <- rhea$diameter > 20
ind_15 <- rhea$diameter > 15

# Sample sizes
nrow(rhea)
sum(ind_15_20)
sum(ind_20)
sum(ind_15)

# Tests to be performed
type_tests <- c("PCvM", "PAD", "PRt")

# Tests
tests_rhea_15_20 <- unif_test(data = rhea$X[ind_15_20, ], type = type_tests,
                              p_value = "asymp", K_max = 5e4)
tests_rhea_20 <- unif_test(data = rhea$X[ind_20, ], type = type_tests,
                           p_value = "asymp", K_max = 5e4)
tests_rhea_15 <- unif_test(data = rhea$X[ind_15, ], type = type_tests,
                           p_value = "asymp", K_max = 5e4)
tests_rhea_15_20
tests_rhea_20
tests_rhea_15

## Venus craters

# Load data
data(venus)

# Add Cartesian coordinates
venus$X <- cbind(cos(venus$theta) * cos(venus$phi),
                 sin(venus$theta) * cos(venus$phi),
                 sin(venus$phi))

# Sample size
nrow(venus)

# Tests
type_tests <- c("PCvM", "PAD", "PRt")
tests_venus <- unif_test(data = venus$X, type = type_tests,
                         p_value = "asymp", K_max = 5e4)
tests_venus
