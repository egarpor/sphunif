
## Retrieve data from https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf

# Extract I (inclination in degs) and Omega (longitude of the ascending
# node in degs) from Table 2a
planets <- matrix(c(7.00559432, 48.33961819,
                    3.39777545, 76.67261496,
                    0, 0,
                    1.85181869, 49.71320984,
                    1.29861416, 100.29282654,
                    2.49424102, 113.63998702,
                    0.77298127, 73.96250215,
                    1.77005520, 131.78635853,
                    17.14104260, 110.30167986),
                  ncol = 2, byrow = TRUE)

# As data.frame
names <- c("Mercury", "Venus", "Earth", "Mars", "Jupiter",
           "Saturn", "Uranus", "Neptune", "Pluto")
planets <- data.frame("planet" = names, "i" = planets[, 1], "om" = planets[, 2])

# Transform degrees to radians
planets$i <- planets$i / 180 * pi
planets$om <- planets$om / 180 * pi

# Save object
save(list = "planets", file = "planets.rda", compress = "bzip2")
