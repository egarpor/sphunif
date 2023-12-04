
## Retrieve data from https://astrogeology.usgs.gov/search/map/Venus/venuscraters

# The source is the USGS Astrogeology Science Center

# Download
download.file(url = paste0("https://astropedia.astrogeology.usgs.gov/",
                           "download/Venus/venuscraters.csv"),
              destfile = "venuscraters.csv")

## Import and preprocess data

# Read data
venus <- read.csv(file = "venuscraters.csv", header = FALSE, sep = ",",
                  dec = ".")[, 1:4]
names(venus)[1:4] <- c("name", "Center_Latitude", "Center_Longitude",
                       "diameter")

# Remove NAs (there are empty rows)
noNA <- complete.cases(cbind(venus$Center_Latitude, venus$Center_Longitude))
venus <- venus[noNA, ]

# Remove initial spaces in names
venus$name <- gsub(pattern = " (.*)", replacement = "\\1", x = venus$name)

# Transform degrees to radians
venus$theta <- (venus$Center_Longitude / 180 * pi) %% (2 * pi)
venus$phi <- venus$Center_Latitude / 180 * pi
venus$Center_Latitude <- NULL
venus$Center_Longitude <- NULL

# Save object
save(list = "venus", file = "venus.rda", compress = "bzip2")
