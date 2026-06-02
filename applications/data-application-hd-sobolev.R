
# ## Import the raw data
#
# # rcosmo is archived on CRAN; install with
# #   install.packages(c("geoR", "FITSio", "mmap"))
# #   install.packages(paste0("https://cran.r-project.org/src/contrib/Archive/",
# #                           "rcosmo/rcosmo_1.1.3.tar.gz"), repos = NULL)
# library(rcosmo)
# library(polykde)
# library(FNN)
#
# # Download the Planck SMICA temperature map (168 MB)
# map_file <- "smica_cmb_1024.fits"
# if (!file.exists(map_file)) {
#
#   downloadCMBMap(foreground = "smica", nside = 1024, destfile = map_file)
#
# }
# cmb <- CMBDataFrame(map_file)
# ordering(cmb) <- "nested" # Keep each coarse cell's children contiguous
#
# # Equal-area grid for the harmonic transform (HEALPix Nside): supports degrees
# # up to ~2 * NSIDE = 128, safely above l_max = 100, divides the native 1024 grid
# NSIDE <- 64L
#
# # Coarsen to the NSIDE grid by averaging HEALPix pixels
# Tt <- cmb$I
# r <- (nside(cmb) / NSIDE)^2
# Tg <- colMeans(matrix(Tt, nrow = r), na.rm = TRUE)
# xyz <- as.data.frame(
#   pix2coords(nside = NSIDE, ordering = "nested"))[, c("x", "y", "z")]
#
# ## Hammer-projection map of the CMB on the NSIDE grid
#
# # 2D coordinates of the grid cells
# H <- sph_to_hammer(as.matrix(xyz))
#
# # Temperature palette in microkelvin
# muK <- Tg * 1e6
# lim <- as.numeric(quantile(abs(muK), 0.99))
# brk <- seq(-lim, lim, length.out = 257)
# pal <- colorRampPalette(
#   c("navy", "steelblue", "white", "tomato", "darkred"))(256)
#
# # Rasterize the map by inverse Hammer projection to avoid the white gaps and
# # overplotting of a projected HEALPix point cloud: build a fine 2D grid, recover
# # each cell's direction on the sphere, and colour it with its nearest HEALPix
# # pixel. The ellipse spanned by the projection has semi-axes 2 * sqrt(2)
# # (horizontal) and sqrt(2) (vertical).
# ax <- 2 * sqrt(2)
# ay <- sqrt(2)
# nx <- 2400L
# ny <- 1200L
# gx <- seq(-ax, ax, length.out = nx)
# gy <- seq(-ay, ay, length.out = ny)
# G <- expand.grid(x = gx, y = gy)
#
# # Keep only the 2D grid points lying inside the projection ellipse
# inside <- (G$x^2 / 8 + G$y^2 / 2) <= 1
#
# # Inverse Hammer projection: plane (x, y) -> longitude lambda, latitude phi
# z <- sqrt(pmax(1 - (G$x / 4)^2 - (G$y / 2)^2, 0))
# lambda <- 2 * atan2(z * G$x, 2 * (2 * z^2 - 1))
# phi <- asin(pmax(pmin(z * G$y, 1), -1))
#
# # Directions on the sphere in the data convention (lon = lambda + pi)
# lon_g <- lambda + pi
# sph <- cbind(cos(phi) * sin(lon_g), cos(phi) * cos(lon_g), sin(phi))
#
# # Nearest HEALPix pixel for each inside grid point, with its temperature
# nn <- get.knnx(as.matrix(xyz), sph[inside, , drop = FALSE], k = 1)$nn.index
# muK_grid <- rep(NA_real_, nx * ny)
# muK_grid[inside] <- pmin(pmax(muK[nn], -lim), lim)
# Z <- matrix(muK_grid, nrow = nx, ncol = ny)
#
# pdf("cmb_hammer.pdf", width = 7, height = 4.1)
# # png("cmb_hammer.png", width = 7, height = 4.1, res = 300, units = "in",
# #     bg = "transparent")
#
# # Plot the rasterized CMB map (cells outside the ellipse stay transparent)
# par(mar = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))
# plot(NA, xlim = c(-ax, ax), ylim = c(-ay, ay), asp = 1, axes = FALSE,
#      xlab = "", ylab = "")
# image(gx, gy, Z, col = pal, breaks = brk, add = TRUE, useRaster = TRUE)
#
# # Meridians and parallels every 30 deg
# lat <- seq(-pi / 2, pi / 2, length.out = 200)
# lon <- seq(0, 2 * pi, length.out = 400)
# grid_xy <- function(lon, lat)
#   sph_to_hammer(cbind(cos(lat) * sin(lon), cos(lat) * cos(lon), sin(lat)))
# for (lon0 in seq(0, 2 * pi, by = pi / 6)) {
#   lines(grid_xy(rep(lon0, 200), lat), col = "gray", lwd = 0.75)
# }
# for (lat0 in seq(-pi / 3, pi / 3, by = pi / 6)) {
#   lines(grid_xy(lon, rep(lat0, 400)), col = "gray", lwd = 0.75)
# }
#
# dev.off()
#
# ## Computation of spherical harmonic coefficients
#
# # Spherical harmonic degrees l = 2, ..., l_max; the monopole (l = 0) and dipole
# # (l = 1) are non-cosmological and excluded, as standard in CMB isotropy
# # (Planck 2015 XVI, l >= 2): https://doi.org/10.1051/0004-6361/201526681
# l_max <- 100
#
# # Standardized coefficients by HEALPix quadrature: divide each degree by its
# # root-mean-square to pool into a single N(0, 1) sample under H0.
# area <- 4 * pi / length(Tg)
# cmb <- unlist(lapply(2:l_max, function(l) {
#   a <- sapply(-l:l,
#               function(m) area * sum(Tg * sphericalHarmonics(l, m, xyz)))
#   a / sqrt(mean(a^2))
# }))
#
# # Save as an rda
# save(list = "cmb", file = "cmb.rda", compress = "bzip2")

## Create the high-dimensional sample

library(sphunif)
library(testthat)
library(goftest)

# Load dataset
stopifnot(packageVersion("sphunif") >= "1.4.4")
data(cmb)

# Balanced high-dimensional design n = d
n <- 100
p <- 101
stopifnot(n * p <= length(cmb))
X <- matrix(cmb[seq_len(n * p)], nrow = n, ncol = p)

# Projections and radii
projs <- X / sqrt(rowSums(X^2))
radii_2 <- rowSums(X^2)

## Tests of uniformity and normality

# Asymptotic standardization of Sobolev statistics
hd_std_asymp <- function(stat, n, p, vk2) {

  # Only accept scalar n and p
  stopifnot(length(n) == 1, length(p) == 1)

  # Remove biases to make it a U-statistic
  bk <- vk2_to_bk(vk2, p = p)
  nonzero_vk2 <- which(vk2 != 0)
  bias <- drop(bk[nonzero_vk2] %*%
                 Gegen_polyn(theta = 0, k = nonzero_vk2, p = p))

  # # \sigma_n^{-1} for a single k0 s.t. v_{k0} != 0
  # # \sqrt{(2k_0!) / d_n^{k_0}} / 2 = \sqrt{k_0! / (2 * d_n^{k_0})}
  # # The division by 2 is done because the computed statistic is 2/n \sum_{i<j}
  # k_0 <- which(vk2 != 0)
  # inv_sigma <- sqrt(factorial(k_0) / (2 * (p - 1)^k_0))

  # General \sigma_n^{-1}
  k <- seq_along(vk2)
  inv_sigma <- 1 / sqrt(2 * sum(exp(2 * log(vk2) +
                                      d_p_k(p = p, k = k, log = TRUE))))

  # Standardize
  return(inv_sigma * (stat - bias))

}

# Verification
test_that("HD-std Sobolev stats coincide with equations (9-10) in the paper", {

  # Weights
  K0 <- 6
  vk2 <- diag(K0)

  # k0-Sobolev uniformity tests with asymptotic standardization
  stats_sphunif <- unif_stat(array(projs, c(n, p, 1)), type = "Sobolev",
                             Sobolev_vk2 = vk2)
  std_stats_sphunif <- sapply(seq_along(stats_sphunif), function(i)
    hd_std_asymp(stats_sphunif[[i]], n = n, p = p, vk2 = vk2[i, ]))

  # Exact form (9) in the paper
  theta_dist <- c(Psi_mat(array(projs, c(n, p, 1))))
  stats_9_exact_factor <- sapply(seq_len(K0), function(k0) {
    bk <- 1 + 2 * k0 / (p - 2)
    bk * sqrt(2 / d_p_k(p = p, k = k0, log = FALSE)) / n *
      sum(Gegen_polyn(theta = theta_dist, k = k0, p = p))
  })

  # Asymptotic form (10) in the paper
  stats_9_asymp_factor <- sapply(seq_len(K0), function(k0) {
    sqrt(2 * factorial(k0) / (p - 1)^k0) / n *
      sum(Gegen_polyn(theta = theta_dist, k = k0, p = p))
  })

  expect_equal(std_stats_sphunif, stats_9_exact_factor)
  expect_equal(std_stats_sphunif, stats_9_asymp_factor, tolerance = 0.05)

})

# k0-Sobolev uniformity tests
K0 <- 6
vk2 <- diag(K0)
proj_test <- unif_stat(array(projs, c(n, p, 1)), type = "Sobolev",
                       Sobolev_vk2 = vk2)
proj_stats <- sapply(seq_along(proj_test),
  function(i) hd_std_asymp(proj_test[[i]], n = n, p = p, vk2 = vk2[i, ]))
p_unif_1 <- pnorm(proj_stats, lower.tail = FALSE)
p_unif_2 <- pchisq(proj_stats^2, df = 1, lower.tail = FALSE)

# Radial Anderson-Darling GoF of ||X_i||^2 ~ chi^2_p (shared by all k0)
ad_test <- ad.test(radii_2, null = "pchisq", df = p)
p_radii <- ad_test$p.value

# Associated tests of Gaussianity: Fisher's aggregation, ~ chi^2_4
p_gauss_1 <- pchisq(-2 * (log(p_radii) + log(p_unif_1)), df = 4,
                    lower.tail = FALSE)
p_gauss_2 <- pchisq(-2 * (log(p_radii) + log(p_unif_2)), df = 4,
                    lower.tail = FALSE)

# Arrange the p-values in a table
tab <- rbind(Uniformity_1 = p_unif_1, Uniformity_2 = p_unif_2,
             Chisquaredness = p_radii,
             Gaussianity_1 = p_gauss_1, Gaussianity_2 = p_gauss_2)
colnames(tab) <- paste("k0 =", seq_len(K0))
print(round(tab, 3))
