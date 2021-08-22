
## Retrieve data from https://ssd.jpl.nasa.gov/sbdb_query.cgi

# The source is NASA's Jet Propusion Laboratory Small-Body Database
# Search Engine

# The database was accessed on the 2020-05-07 using the following query:
#
# Step 1.
# - Limit by object type/group:
#   * object group:    [x] All objects [ ] NEOs      [ ] PHAs
#   * object kind:     [ ] All objects [ ] Asteroids [x] Comets
#   * numbered state:  [x] All objects [ ] Numbered  [ ] Unnumbered
# - Limit to selected orbit class(es):
#   * Asteroid Orbit Classes: mark none
#   * Comet Orbit Classes: mark all
# - Limit by object characteristics: add no constraints
#
# Step 2.
# - Object Fields: "object internal database ID", "object full name/designation"
# - Orbital and Model Parameter Fields:
#   "[i] inclination; angle with respect to x-y ecliptic plane (deg)",
#   "longitude of the ascending node (deg)", "sidereal orbital period (years)",
#   "orbit classification"
#
# Step 3, 4, 5. Not relevant.
#
# Step 6. [ ] HTML (default) [x] CSV (download)
#
# The returned dataset was saved as "comets-accessed-2020-05-07.csv". The urls
# with the database query (done in 2020-05-07) are
# CSV: https://ssd.jpl.nasa.gov/sbdb_query.cgi?obj_group=all;obj_kind=com;obj_numbered=all;com_orbit_class=HYP;com_orbit_class=ETc;com_orbit_class=PAR;com_orbit_class=CTc;com_orbit_class=JFC;com_orbit_class=JFc;com_orbit_class=HTC;com_orbit_class=COM;OBJ_field=0;ORB_field=0;table_format=CSV;format_option=full;query=Generate%20Table;c_fields=AaAcBjBkBsCiAp;c_sort=;.cgifields=format_option;.cgifields=ast_orbit_class;.cgifields=table_format;.cgifields=obj_kind;.cgifields=obj_group;.cgifields=obj_numbered;.cgifields=com_orbit_class
# HTML: https://ssd.jpl.nasa.gov/sbdb_query.cgi?obj_group=all;obj_kind=com;obj_numbered=all;com_orbit_class=HYP;com_orbit_class=ETc;com_orbit_class=PAR;com_orbit_class=CTc;com_orbit_class=JFC;com_orbit_class=JFc;com_orbit_class=HTC;com_orbit_class=COM;OBJ_field=0;ORB_field=0;table_format=HTML;max_rows=500;format_option=full;query=Generate%20Table;c_fields=AaAcBjBkBsCiAp;c_sort=;.cgifields=format_option;.cgifields=ast_orbit_class;.cgifields=table_format;.cgifields=obj_kind;.cgifields=obj_group;.cgifields=obj_numbered;.cgifields=com_orbit_class

# Alternatively
download.file(url = paste0("https://ssd.jpl.nasa.gov/sbdb_query.cgi?",
                           "obj_group=all;obj_kind=com;obj_numbered=all;",
                           "com_orbit_class=HYP;com_orbit_class=ETc;",
                           "com_orbit_class=PAR;com_orbit_class=CTc;",
                           "com_orbit_class=JFC;com_orbit_class=JFc;",
                           "com_orbit_class=HTC;com_orbit_class=COM;",
                           "OBJ_field=0;ORB_field=0;table_format=CSV;",
                           "format_option=full;query=Generate%20Table;",
                           "c_fields=AaAcBjBkBsCiAp;c_sort=;",
                           ".cgifields=format_option;",
                           ".cgifields=ast_orbit_class;",
                           ".cgifields=table_format;.cgifields=obj_kind;",
                           ".cgifields=obj_group;.cgifields=obj_numbered;",
                           ".cgifields=com_orbit_class"),
              destfile = paste0("comets-accessed-", Sys.Date(), ".csv"))

# The classes of comet orbits are:
#
# - "COM" = Comet. Comet orbit not matching any defined orbit class.
# - "CTc" = Chiron-type Comet. Chiron-type comet, as defined by Levison and
#           Duncan (T_Jupiter > 3; a > a_Jupiter).
# - "ETc" = Encke-type Comet. Encke-type comet, as defined by Levison and
#           Duncan (T_Jupiter > 3; a < a_Jupiter).
# - "HTC" = Halley-type Comet*. Halley-type comet, classical definition
#           (20y < P < 200y).
# - "HYP" = Hyperbolic Comet. Comets on hyperbolic orbits (e > 1.0).
# - "JFc" = Jupiter-family Comet. Jupiter-family comet, as defined by
#           Levison and Duncan (2 < T_Jupiter < 3).
# - "JFC" = Jupiter-family Comet*. Jupiter-family comet, classical
#           definition (P < 20y).
# - "PAR" = Parabolic Comet. Comets on parabolic orbits (e = 1.0).

## Import and preprocess data

# Read data
comets <- read.csv(file = "comets-accessed-2020-05-07.csv", header = TRUE)

# Transform degrees to radians
comets$i <- comets$i / 180 * pi
comets$om <- comets$om / 180 * pi

# Normal vectors to the ecliptic planes of the comets
comets_normal <- cbind(sin(comets$i) * sin(comets$om),
                       -sin(comets$i) * cos(comets$om),
                       cos(comets$i))

## Compare with the data in Cuesta-Albertos et al. (2009)

# Read data directly retrieved from the authors' file
comets_ccf2009 <- readxl::read_xls("LP-unumb.xls")
comets_ccf2009_normal <- cbind(comets_ccf2009$...5,
                               comets_ccf2009$...6,
                               comets_ccf2009$...7)
ccf2009 <- match(comets_ccf2009$id, comets$id)

# Check equality of observations
max(rowSums((comets_normal[ccf2009, ] - comets_ccf2009_normal)^2))

# Add flag to larger database
comets$ccf2009 <- FALSE
comets$ccf2009[ccf2009] <- TRUE

# Save object
save(list = "comets", file = "comets.rda", compress = "bzip2")
