
## Retrieve data from https://planetarynames.wr.usgs.gov/AdvancedSearch

# The source is the Gazetteer of Planetary Nomenclature of the
# International Astronomical Union (IAU)

# The database was accessed on the 2020-05-31 using the following query:
# System: All.
# Target: All.
# Coordinate System: +East. 0 - 360. Planetocentric.
# Latitude Longitude: Do not enter anything.
# Feature Type: Crater, craters.
# Feature Name: Do not enter anything.
# Approval Status: Adopted by IAU.
# Feature Diameter (km): Do not select anything.
# Approval Date: Do not enter anything.
# Continent: All.
# Ethnic/Cultural Group or Country: All.
# Reference: All.
# Columns to include:
# [x] Feature ID
# [x] Feature Name
# [ ] Clean Feature Name
# [x] Target
# [x] Diameter
# [x] Center Lat/Lon
# [x] Lat/Lon Boundaries
# [ ] Coordinate System
# [ ] Continent/Ethnicity
# [ ] Feature Type
# [ ] Feature Type Code
# [ ] Quad
# [ ] Approval Status
# [ ] Approval Date
# [ ] Reference
# [ ] Origin
# [ ] Additional Info
# [ ] Last Updated
# Sorted by: Feature ID. Ascending.
# Output Format: CSV - (Comma Separated Values) for importing into Excel

# Read data
craters <- read.csv(file = "craters-accessed-2020-05-31.csv", header = TRUE)
craters$X <- NULL

## Preprocess data

# Transform degrees to radians
craters$phi <- craters$Center_Latitude / 180 * pi
craters$theta <- (craters$Center_Longitude / 180 * pi) %% (2 * pi)
craters$Center_Latitude <- NULL
craters$Center_Longitude <- NULL

# Rename
names(craters)[1:4] <- c("ID", "name", "target", "diameter")

# Type of body (planet, moon, asteroid, ...)
body_type <- read.csv(file = "bodies-type.csv", sep = ",", dec = ".")[, 1:2]
craters$target_type <- body_type$Type[match(craters$target, body_type$Body)]

# Reorder
craters <- craters[, c(1:3, 7, 4:6)]

# Save object
save(list = "craters", file = "craters.rda", compress = "bzip2")
