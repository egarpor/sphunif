
library(rjson)

## Retrieve data from https://ssd.jpl.nasa.gov/tools/sbdb_query.html

# The source is NASA's Jet Propulsion Laboratory Small-Body Database
# Search Engine.

# The database was accessed on the 2022-05-28 using the following query:
#
# - Limit by Object Kind/Group:
#   * object group:    [x] all objects [ ] NEOs      [ ] PHAs
#   * object kind:     [ ] all objects [ ] asteroids [x] comets
#   * numbered state:  [x] all objects [ ] numbered  [ ] unnumbered
#   * comet fragments: mark or leave unmarked depending if comet fragments
#     are to be retrieved also (see below).
#   * satellite(s): leave unmarked
# - Limit by Orbit Class:
#   * Asteroid Orbit Classes: mark none
#   * Comet Orbit Classes: mark all
# - Custom Object/Orbit Constraints: add no constraints
# - Output Selection Controls:
#   Select the following fields from Available Fields:
#   * "SPK-ID": object primary SPK-ID.
#   * "object fullname": object full name/designation.
#   * "prim. desig.": object primary designation.
#   * "diameter": object diameter (from equivalent sphere) (km).
#   * "e": eccentricity.
#   * "a": semi-major axis (au).
#   * "i": inclination; angle with respect to x-y ecliptic plane (deg).
#   * "node": longitude of the ascending node (deg).
#   * "peri": argument of perihelion (deg).
#   * "period": sidereal orbital period (years).
#   * "orbit class": classification.
#   * "date of first obs.": date of first observation used in the orbit fit (UT).
#   * "date of last obs.": date of last observation used in the orbit fit (UT).
# - Mark Full Precision in Output Fields
#
# The classes of comet orbits are:
#
# - "COM" = Comet. Comet orbit not matching any defined orbit class.
# - "CTc" = Chiron-type Comet. Chiron-type comet, as defined by Levison and
#           Duncan (T_Jupiter > 3; a > a_Jupiter).
# - "ETc" = Encke-type Comet. Encke-type comet, as defined by Levison and
#           Duncan (T_Jupiter > 3; a < a_Jupiter).
# - "HTC" = Halley-type Comet. Halley-type comet, classical definition
#           (20y < P < 200y).
# - "HYP" = Hyperbolic Comet. Comets on hyperbolic orbits (e > 1.0).
# - "JFc" = Jupiter-family Comet. Jupiter-family comet, as defined by
#           Levison and Duncan (2 < T_Jupiter < 3).
# - "JFC" = Jupiter-family Comet. Jupiter-family comet, classical
#           definition (P < 20y).
# - "PAR" = Parabolic Comet. Comets on parabolic orbits (e = 1.0).

# The urls with the database query from the above options, as done in
# 2022-05-28 through https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html
# (version 1.0), are (for all comets and excluding comet fragments):
url_all <- paste0("https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=spkid,",
                  "full_name,pdes,diameter,e,a,i,om,w,per_y,class,first_obs,",
                  "last_obs&full-prec=true&sb-class=HYP,PAR,COM,JFC,HTC,ETc,",
                  "CTc,JFc&sb-kind=c&www=1")
url_nofr <- paste0("https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=spkid,",
                   "full_name,pdes,diameter,e,a,i,om,w,per_y,class,first_obs,",
                   "last_obs&full-prec=true&sb-class=HYP,PAR,COM,JFC,HTC,ETc,",
                   "CTc,JFc&sb-kind=c&sb-xfrag=1&www=1")

# Download the JSON including comet fragments
json_all <- rjson::fromJSON(file = url_all, simplify = FALSE)

# Put in a tidy data frame
df <- sapply(seq_along(json_all$fields), function(i) {
  sapply(json_all$data, function(x) {
    y <- x[[i]]
    y[which(is.null(y))] <- NA
    return(unlist(y))
    })
  })
colnames(df) <- unlist(json_all$fields)
df <- as.data.frame(df)
for (fac in c("id", "full_name", "pdes", "class")) {
  df[[fac]] <- as.factor(df[[fac]])
}
for (num in c("spkid", "diameter", "e", "a", "i", "om", "w", "per_y")) {
  df[[num]] <- as.numeric(df[[num]])
}

# Download JSON with comet fragments
json_nofr <- rjson::fromJSON(file = url_nofr, simplify = FALSE)

# Add flag with comet fragments
spkid_nofr <- sapply(json_nofr$data, function(x) {
  y <- x[[2]]
  y[which(is.null(y))] <- NA
  return(unlist(y))
})
df$frag <- is.na(match(x = df$spkid, table = spkid_nofr))

# Put the flag after pdes
df <- df[, c(1:4, 15, 5:14)]

## Retrieve discovery date from https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html

# The discovery date is not accessible through the previous database. One needs
# to query the NASA's Jet Propulsion Laboratory Small-Body Database Lookup via
# its API using the comet's primary identifiers previously obtained.

# Access discovery dates retrieving the corresponding JSON
disc <- rep(NA, nrow(df))
pb <- txtProgressBar(style = 3)
for (i in seq_len(nrow(df))) {
  pdes <- gsub(x = df$pdes[i], pattern = " ", replacement = "%20")
  url_disc <- paste0("https://ssd-api.jpl.nasa.gov/sbdb.api?des=",
                     pdes, "&discovery=1")
  date <- rjson::fromJSON(file = url_disc)$discovery$date
  disc[i] <- ifelse(is.null(date), NA, date)
  setTxtProgressBar(pb, i / nrow(df))
}

# Add discovery date to the data frame
df$discovery <- disc

# Put the discovery before first_obs
df <- df[, c(1:13, 16, 14:15)]

## Save data frame as a csv

write.csv(x = df, file = "comets-accessed-2022-05-28.csv", row.names = FALSE)

## Import and preprocess data

# Read data
comets <- read.csv(file = "comets-accessed-2022-05-28.csv", header = TRUE,
                   stringsAsFactors = TRUE)

# Transform degrees to radians
comets$i <- comets$i / 180 * pi
comets$om <- comets$om / 180 * pi
comets$w <- comets$w / 180 * pi

# Process observation dates as Date
comets$first_obs <- as.Date(comets$first_obs)
comets$last_obs <- as.Date(comets$last_obs)
# Two BC dates are lost in this process:
# 1000590 C/-43 K1
# 1000589 C/-146 P1

# Process discovery dates as Date
comets$discovery <- as.Date(comets$discovery, tryFormats = "%Y-%B-%d")

# Most of the discovery dates are missing... However, we know from
# https://www.iau.org/public/themes/naming/#comets that:
# "When a periodic comet is observed after its second apparition, the IAU’s
# Minor Planet Center (MPC) gives it a sequential number indicating the order
# of the discovery." Therefore, the pdes and the id fields encode the
# sequential order of the discovery of the comets. Notice that the dataset is
# sorted according to id. We remove the discovery field since it is not very
# relevant.
comets$discovery <- NULL

# Compute the normal vectors to the ecliptic planes of the comets
comets_normal <- cbind(sin(comets$i) * sin(comets$om),
                       -sin(comets$i) * cos(comets$om),
                       cos(comets$i))

## Compare with the data in Cuesta-Albertos et al. (2009)

# Read data directly retrieved from the authors' file
comets_ccf09 <- readxl::read_xls("LP-unumb.xls")
comets_ccf09_normal <- cbind(comets_ccf09$...5,
                             comets_ccf09$...6,
                             comets_ccf09$...7)
ccf09 <- match(comets_ccf09$full_name,
               substr(comets$full_name, start = 6, stop = 100))

# Matches are correct and the same if done with the id's
stopifnot(!anyNA(ccf09))
stopifnot(all(ccf09 == match(comets_ccf09$id, comets$id))) 

# Check cosine equality of observations
summary(1 - rowSums(comets_normal[ccf09, ] * comets_ccf09_normal))

# The biggest difference (0.035608) happens for "C/1905 F1 (Giacobini)", which
# seems to have been updated in 2021-May-05 09:27:44 (solution date):
# https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?des=1905%20F1
# Excluding this comet, the maximum difference in the normal vectors is
max(1 - rowSums(comets_normal[ccf09, ] * comets_ccf09_normal)[-20])

# Add flag to larger database
comets$ccf09 <- FALSE
comets$ccf09[ccf09] <- TRUE

# Trim whitespaces in names
comets$full_name <- stringr::str_trim(comets$full_name)

# Save object
save(list = "comets", file = "comets.rda", compress = "bzip2")
