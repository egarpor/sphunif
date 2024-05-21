
# Clean workspace
rm(list = ls())

# Load packages
library(sphunif)

# Load data
data("venus")

# Add Cartesian coordinates
venus$X <- cbind(cos(venus$theta) * cos(venus$phi),
                 sin(venus$theta) * cos(venus$phi),
                 sin(venus$phi))

# Tests
unif_test(data = venus$X, type = c("PCvM", "PAD", "PRt"), p_value = "asymp")
