
## Standardized spherical-harmonic coefficients of the Planck SMICA CMB map

# The map is the foreground-cleaned SMICA temperature map of the Planck
# mission, read with rcosmo on the HEALPix pixelization of the sphere.
# rcosmo is archived on CRAN; install with
#   install.packages(c("geoR", "FITSio", "mmap"))
#   install.packages(paste0("https://cran.r-project.org/src/contrib/Archive/",
#                           "rcosmo/rcosmo_1.1.3.tar.gz"), repos = NULL)

library(rcosmo)

## Import the raw data

# Download the Planck SMICA temperature map (168 MB)
map_file <- "smica_cmb_1024.fits"
if (!file.exists(map_file)) {

  downloadCMBMap(foreground = "smica", nside = 1024, destfile = map_file)

}
cmb_map <- CMBDataFrame(map_file)
ordering(cmb_map) <- "nested" # Keep each coarse cell's children contiguous

# Equal-area grid for the harmonic transform (HEALPix Nside): supports degrees
# up to ~2 * NSIDE = 128, safely above l_max = 100, divides the native 1024 grid
NSIDE <- 64

# Coarsen to the NSIDE grid by averaging HEALPix pixels
Tt <- cmb_map$I
r <- (nside(cmb_map) / NSIDE)^2
Tg <- colMeans(matrix(Tt, nrow = r), na.rm = TRUE)
xyz <- as.data.frame(
  pix2coords(nside = NSIDE, ordering = "nested"))[, c("x", "y", "z")]

## Computation of spherical harmonic coefficients

# Spherical harmonic degrees l = 2, ..., l_max; the monopole (l = 0) and dipole
# (l = 1) are non-cosmological and excluded, as standard in CMB isotropy
# (Planck 2015 XVI, l >= 2): https://doi.org/10.1051/0004-6361/201526681
l_max <- 100

# Standardized coefficients by HEALPix quadrature: divide each degree by its
# root-mean-square to pool into a single N(0, 1) sample under H0.
area <- 4 * pi / length(Tg)
cmb <- unlist(lapply(2:l_max, function(l) {
  a <- sapply(-l:l,
              function(m) area * sum(Tg * sphericalHarmonics(l, m, xyz)))
  a / sqrt(mean(a^2))
}))

# Save object
save(list = "cmb", file = "cmb.rda", compress = "bzip2")
