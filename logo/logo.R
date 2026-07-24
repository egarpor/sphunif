
# Hex sticker for the 'sphunif' package: a globe (S^2) with a uniform random
# sample (drawn with sphunif's own r_unif_sph()) joined by great-circle arcs --
# the minimum spanning tree of geodesic distances -- over a Lambert-shaded
# sphere that draws only the front hemisphere. Built as an orthographic ggplot,
# then framed by hexSticker.

library(hexSticker)
library(ggplot2)
library(sphunif)

# Shared logo standards
font <- "Aller_Rg"
name_size <- 31.2
url_size <- 9.0
url_x <- 1.00
url_y <- 0.08
url_angle <- 30
hex_border <- 1.5
dpi <- 600

# Tunable parameters
seed <- 11 # fixed for reproducibility
n_pts <- 18 # uniform sample size on S^2
n_grid <- 320 # sphere-shading raster resolution
n_arc <- 80 # samples per great-circle arc
n_circ <- 240 # samples per wireframe circle

# Viewing rotation (radians): tilt about x, yaw about y, spin about z.
tilt <- -0.42
yaw <- 0.55
spin <- 0.30

# Light direction for the Lambert-shaded sphere (view space, unit-ish).
light <- c(-0.42, 0.52, 0.74)

# Palette (deep-blue globe, viridis points, light arcs)
sphere_ramp <- c("#07101F", "#0C1F3E", "#173E73", "#2E6BA6", "#6BA6D6",
                 "#A9D2EE")
col_wire <- "#8FC0E6"
col_arc <- "#F2FAFF"
col_arc_case <- "#0A1B33"
col_rim <- "#7FB4DD"
col_pt_out <- "#0A1730"
h_fill <- "#0C1626"
h_color <- "#3E86C4"
p_color <- "#EAF2FB"
u_color <- "#8FB6DB"

# Geometry helpers
rot_x <- function(a) matrix(c(1, 0, 0, 0, cos(a), -sin(a), 0, sin(a), cos(a)),
                            3, 3, byrow = TRUE)
rot_y <- function(a) matrix(c(cos(a), 0, sin(a), 0, 1, 0, -sin(a), 0, cos(a)),
                            3, 3, byrow = TRUE)
rot_z <- function(a) matrix(c(cos(a), -sin(a), 0, sin(a), cos(a), 0, 0, 0, 1),
                            3, 3, byrow = TRUE)

# Composite viewing rotation.
rot <- rot_x(tilt) %*% rot_y(yaw) %*% rot_z(spin)

# Orthographic projection: rotate rows of pts by rot, keep (x, y); depth = z
# (depth > 0 is the front, viewer-facing hemisphere).
project <- function(pts) {
  rotated <- pts %*% t(rot)
  data.frame(x = rotated[, 1], y = rotated[, 2], depth = rotated[, 3])
}

# Take a 3D curve (rows of pts), project it, and return only the front-facing
# portions, split into contiguous runs so ggplot does not connect across the
# parts hidden behind the globe.
front_path <- function(pts, id) {
  q <- project(pts)
  front <- q$depth > 1e-9
  if (!any(front)) return(NULL)
  runs <- rle(front)
  q$grp <- paste(id, rep(seq_along(runs$lengths), runs$lengths), sep = "_")
  q[front, c("x", "y", "grp")]
}

# Great-circle (slerp) arc between two unit vectors a and b.
slerp_arc <- function(a, b, n = n_arc) {
  omega <- acos(max(-1, min(1, sum(a * b))))
  t <- seq(0, 1, length.out = n)
  if (omega < 1e-8) return(matrix(a, n, 3, byrow = TRUE))
  (outer(sin((1 - t) * omega), a) + outer(sin(t * omega), b)) / sin(omega)
}

# Minimum spanning tree (Prim) over geodesic distances: the sparse graph of
# arcs joining observations.
prim_mst <- function(x) {
  dist_geo <- acos(pmin(pmax(tcrossprod(x), -1), 1))
  n <- nrow(x)
  intree <- c(TRUE, rep(FALSE, n - 1))
  edges <- matrix(NA_integer_, n - 1, 2)
  for (k in seq_len(n - 1)) {
    sub <- dist_geo[intree, !intree, drop = FALSE]
    idx <- which(sub == min(sub), arr.ind = TRUE)[1, ]
    i <- which(intree)[idx[1]]
    j <- which(!intree)[idx[2]]
    edges[k, ] <- c(i, j)
    intree[j] <- TRUE
  }
  edges
}

# Sphere body: Lambert-shaded raster over the unit disk
g <- seq(-1, 1, length.out = n_grid)
grid <- expand.grid(x = g, y = g)
rho2 <- grid$x^2 + grid$y^2
inside <- rho2 <= 1
z <- sqrt(pmax(0, 1 - rho2)) # front-hemisphere height (screen space)
light_unit <- light / sqrt(sum(light^2))
ndotl <- grid$x * light_unit[1] + grid$y * light_unit[2] + z * light_unit[3]
ambient <- 0.30
shade <- ambient + (1 - ambient) * pmax(0, ndotl)

# Map the shading to sphere_ramp colours
ramp <- grDevices::colorRamp(sphere_ramp)
fill <- rep(NA_character_, nrow(grid))
rgbm <- ramp(pmin(1, pmax(0, shade[inside])))
fill[inside] <- grDevices::rgb(rgbm[, 1], rgbm[, 2], rgbm[, 3],
                               maxColorValue = 255)
sphere_df <- data.frame(x = grid$x, y = grid$y, fill = fill)

# Wireframe: equator, a few meridians and parallels (front-facing only)
u <- seq(0, 2 * pi, length.out = n_circ)
wire_list <- list()

# Equator + parallels at latitudes +/- 40 deg.
for (lat in c(0, 40, -40) * pi / 180) {
  parallel_xyz <- cbind(cos(lat) * cos(u), cos(lat) * sin(u),
                        rep(sin(lat), length(u)))
  wire_list[[length(wire_list) + 1]] <-
    front_path(parallel_xyz, paste0("par", round(lat, 3)))
}
# Meridians at longitudes 0, 60, 120 deg.
for (lon in c(0, 60, 120) * pi / 180) {
  meridian_xyz <- cbind(sin(u) * cos(lon), sin(u) * sin(lon), cos(u))
  wire_list[[length(wire_list) + 1]] <-
    front_path(meridian_xyz, paste0("mer", round(lon, 3)))
}
wire_df <- do.call(rbind, wire_list)

# Observations: uniform sample on S^2 (sphunif) + MST great-circle arcs
set.seed(seed)
sample_sph <- r_unif_sph(n = n_pts, p = 3)[, , 1]

# Build the MST and its great-circle arcs
edges <- prim_mst(sample_sph)
arc_list <- lapply(seq_len(nrow(edges)), function(k) {
  front_path(slerp_arc(sample_sph[edges[k, 1], ], sample_sph[edges[k, 2], ]),
             paste0("arc", k))
})
arc_df <- do.call(rbind, arc_list)

# Points: keep front hemisphere, colour by latitude via viridis.
qp <- project(sample_sph)
qp$lat <- sample_sph[, 3]
front_pts <- qp[qp$depth > 0, ]
# begin = 0.22 lifts the dark end so low-latitude points stay visible on blue.
pal <- viridisLite::viridis(256, begin = 0.22)
front_pts$fill <-
  pal[pmax(1, pmin(256, round(1 + 255 * (front_pts$lat + 1) / 2)))]

# Silhouette rim.
rim <- data.frame(x = cos(u), y = sin(u))

# Assemble the subplot
subplot <- ggplot() +
  geom_raster(data = sphere_df, aes(x, y, fill = fill)) +
  geom_path(data = wire_df, aes(x, y, group = grp),
            colour = col_wire, linewidth = 0.30, alpha = 0.28,
            lineend = "round") +
  geom_path(data = arc_df, aes(x, y, group = grp),
            colour = col_arc_case, linewidth = 1.15, alpha = 0.55,
            lineend = "round") +
  geom_path(data = arc_df, aes(x, y, group = grp),
            colour = col_arc, linewidth = 0.55, alpha = 0.95,
            lineend = "round") +
  geom_path(data = rim, aes(x, y), colour = col_rim,
            linewidth = 0.5, alpha = 0.55) +
  geom_point(data = front_pts, aes(x, y, fill = fill),
             shape = 21, size = 2.4, stroke = 0.35, colour = col_pt_out) +
  scale_fill_identity() +
  coord_fixed(xlim = c(-1.08, 1.08), ylim = c(-1.08, 1.08), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

# Hex sticker
dir.create("logo", showWarnings = FALSE)
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

# Draw the hexagon, wordmark and URL, then mirror to man/figures
sticker(
  subplot = subplot, s_x = 1, s_y = 1.16, s_width = 1.26, s_height = 1.26,
  package = "sphunif", p_x = 1, p_y = 0.40, p_size = name_size,
  p_color = p_color, p_family = font,
  h_fill = h_fill, h_color = h_color, h_size = hex_border,
  spotlight = FALSE, white_around_sticker = FALSE,
  url = "github.com/egarpor/sphunif",
  u_x = url_x, u_y = url_y, u_angle = url_angle, u_size = url_size,
  u_color = u_color, u_family = font,
  dpi = dpi, filename = "logo/logo.png"
)
file.copy("logo/logo.png", "man/figures/logo.png", overwrite = TRUE)
