
# Clean workspace
rm(list = ls())

# Load packages
library(sphunif)

## Craters named by the IAU

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
