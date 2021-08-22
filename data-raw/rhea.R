
## Retrieve data from https://astrogeology.usgs.gov/search/map/Venus/venuscraters

# The source is https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JE004940

# Download
download.file(url = paste0("https://agupubs.onlinelibrary.wiley.com/action/",
                           "downloadSupplement?doi=10.1002%2F2015JE004940",
                           "&file=jgre20485-sup-0002-TableS1.txt"),
              destfile = "jgre20485-sup-0002-tables1.txt")

## Import and preprocess data

# Read data
rhea <- read.table(file = "jgre20485-sup-0002-TableS1.txt", skip = 2,
                   sep = "\t", dec = ".")
names(rhea) <- c("Center_Longitude", "Center_Latitude", "diameter", "name")

# Transform degrees to radians
rhea$theta <- (rhea$Center_Longitude / 180 * pi) %% (2 * pi)
rhea$phi <- rhea$Center_Latitude / 180 * pi
rhea$Center_Latitude <- NULL
rhea$Center_Longitude <- NULL

# Reorder
rhea <- rhea[, c(2, 1, 3:4)]

# Save object
save(list = "rhea", file = "rhea.rda", compress = "bzip2")
