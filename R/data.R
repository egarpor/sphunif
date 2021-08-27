

#' @title Comet orbits
#'
#' @description Comet orbits data from the
#' \href{https://ssd.jpl.nasa.gov/sbdb_query.cgi}{
#' JPL Small-Body Database Search Engine}. The normal vector of a comet orbit
#' represents is a vector on \eqn{S^2}.
#'
#' @docType data
#' @format A data frame with 3633 rows and 8 variables:
#' \describe{
#'   \item{id}{database ID.}
#'   \item{full_name}{full name/designation following the
#'   \href{https://www.iau.org/public/themes/naming/#comets}{
#'   IUA naming convention}.}
#'   \item{i}{inclination; the orbit angle with respect to the ecliptic plane,
#'   in radians in \eqn{[0, \pi]}.}
#'   \item{om}{longitude of the ascending node; the angle between the normal
#'   vector of the orbit and the normal vector of the ecliptic plane, in
#'   radians in \eqn{[0, 2\pi)}.}
#'   \item{per_y}{sidereal orbital period (in years).}
#'   \item{class}{orbit classification. A factor with levels given below.}
#'   \item{diameter}{diameter from equivalent sphere (in km).}
#'   \item{ccf2009}{flag indicating if the comet was considered in the data
#'   application in Cuesta-Albertos et al. (2009); see details below.}
#' }
#' @details
#' The normal vector to the ecliptic plane of the comet with inclination
#' \eqn{i} and longitude of the ascending node \eqn{\omega} is
#' \deqn{(\sin(i) \sin(\omega), -\sin(i) \cos(\omega), \cos(i))'.}{
#' (sin(i) sin(\omega), -sin(i) cos(\omega), cos(i))'.}
#'
#' A prograde comet has positive \eqn{\cos(i)}{cos(i)}, negative
#' \eqn{\cos(i)}{cos(i)} represents a retrograde comet.
#'
#' \code{class} has the following levels:
#' \itemize{
#'  \item \code{COM}: comet orbit not matching any defined orbit class.
#'  \item \code{CTc}: Chiron-type comet, as defined by Levison and Duncan
#'  (T_Jupiter > 3; a > a_Jupiter).
#'  \item \code{ETc}: Encke-type comet, as defined by Levison and Duncan
#'  (T_Jupiter > 3; a < a_Jupiter).
#'  \item \code{HTC}: Halley-type comet, classical definition (20y < P < 200y).
#'  \item \code{HYP}: comets on hyperbolic orbits.
#'  \item \code{JFc}: Jupiter-family comet, as defined by Levison and Duncan
#'  (2 < T_Jupiter < 3).
#'  \item \code{JFC}: Jupiter-family comet, classical definition (P < 20y).
#'  \item \code{PAR}: comets on parabolic orbits.
#'}
#' Hyperbolic and parabolic comets are not periodic; only elliptical comets
#' are periodic.
#'
#' The \code{ccf2009} variable gives the observations considered in
#' Cuesta-Albertos et al. (2009) after fetching in the database in 2007-12-14
#' for the comets such that
#' \code{!(class \%in\% c("HYP", "PAR")) & per_y >= 200 & !numbered}. A
#' periodic comet is \code{numbered} by the IUA only after its second
#' perihelion passage, and then its \code{id} starts with \code{c}. Due to the
#' dynamic nature of the data, more comets were added to the database since
#' 2007 and also some observations were updated.
#'
#' The script performing the data preprocessing is available at
#' \href{https://github.com/egarpor/sphunif/blob/master/data-raw/comets.R}{
#' \code{comets.R}}. The data was retrieved on 2020-05-07.
#' @source \url{https://ssd.jpl.nasa.gov/sbdb_query.cgi}
#' @references
#' Cuesta-Albertos, J. A., Cuevas, A., Fraiman, R. (2009) On projection-based
#' tests for directional and compositional data. \emph{Statistics and
#' Computing}, 19:367--380. \doi{10.1007/s11222-008-9098-3}
#' @examples
#' # Load data
#' data("comets")
#'
#' # Add normal vectors
#' comets$normal <- cbind(sin(comets$i) * sin(comets$om),
#'                        -sin(comets$i) * cos(comets$om),
#'                        cos(comets$i))
#'
#' # Add numbered information
#' comets$numbered <- substr(comets$id, 1, 1) == "c"
#'
#' # Tests to be performed
#' type_tests <- c("PCvM", "PAD", "PRt")
#'
#' # Excluding the C/1882 R1-X (Great September comet) records with X = B, C, D
#' comets_ccf2009 <- comets[comets$ccf2009, ][-c(13:15), ]
#'
#' # Sample size
#' nrow(comets_ccf2009)
#'
#' # Tests for the data in Cuesta-Albertos et al. (2009)
#' tests_ccf2009 <- unif_test(data = comets_ccf2009$normal, type = type_tests,
#'                            p_value = "asymp")
#' tests_ccf2009
"comets"


#' @title Planet orbits
#'
#' @description Planet orbits data from the
#' \href{https://ssd.jpl.nasa.gov/?planet_pos}{
#' JPL Keplerian Elements for Approximate Positions of the Major Planets}.
#' The normal vector of a planet orbit represents is a vector on \eqn{S^2}.
#'
#' @docType data
#' @format A data frame with 9 rows and 3 variables:
#' \describe{
#'   \item{planet}{names of the planets and Pluto.}
#'   \item{i}{inclination; the orbit angle with respect to the ecliptic plane,
#'   in radians in \eqn{[0, \pi]}.}
#'   \item{om}{longitude of the ascending node; the angle between the normal
#'   vector of the orbit and the normal vector of the ecliptic plane, in
#'   radians in \eqn{[0, 2\pi)}.}
#' }
#' @details
#' The normal vector to the ecliptic plane of the planet with inclination
#' \eqn{i} and longitude of the ascending node \eqn{\omega} is
#' \deqn{(\sin(i) \sin(\omega), -\sin(i) \cos(\omega), \cos(i))'.}{
#' (sin(i) sin(\omega), -sin(i) cos(\omega), cos(i))'.}
#'
#' The script performing the data preprocessing is available at
#' \href{https://github.com/egarpor/sphunif/blob/master/data-raw/planets.R}{
#' \code{planets.R}}. The data was retrieved on 2020-05-16.
#' @source Table 2b in \url{https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf}
#' @examples
#' # Load data
#' data("planets")
#'
#' # Add normal vectors
#' planets$normal <- cbind(sin(planets$i) * sin(planets$om),
#'                        -sin(planets$i) * cos(planets$om),
#'                        cos(planets$i))
#'
#' # Tests to be performed
#' type_tests <- c("PCvM", "PAD", "PRt")
#'
#' # Tests with Pluto
#' unif_test(data = planets$normal, type = type_tests, p_value = "MC")
#'
#' # Tests without Pluto
#' unif_test(data = planets$normal[-9, ], type = type_tests, p_value = "MC")
"planets"


#' @title Craters named by the IUA
#'
#' @description \emph{Named} craters of the Solar System by the
#' \href{https://planetarynames.wr.usgs.gov}{Gazetteer of Planetary
#' Nomenclature} of the International Astronomical Union (IUA).
#'
#' @docType data
#' @format A data frame with 5235 rows and 7 variables:
#' \describe{
#'   \item{ID}{database ID.}
#'   \item{name}{name of the crater.}
#'   \item{target}{name of the celestial body. A factor with 43 levels,
#'   such as \code{"Moon"}, \code{"Venus"}, or \code{"Europa"}.}
#'   \item{target_type}{type of celestial body. A factor with 3 levels:
#'   \code{"Planet"}, \code{"Moon"}, \code{"Dwarf planet"}, or
#'   \code{"Asteroid"}.}
#'   \item{diameter}{diameter of the crater (in km).}
#'   \item{theta}{longitude angle \eqn{\theta \in [0, 2\pi)} of the
#'   crater center.}
#'   \item{phi}{latitude angle \eqn{\phi \in [-\pi/2, \pi/2]} of the
#'   crater center.}
#' }
#' @details
#' "Craters" are understood in the Gazetteer of Planetary Nomenclature as
#' roughly circular depressions resulting from impact or volcanic activity
#' (the geological origin is
#' \href{https://planetarynames.wr.usgs.gov/DescriptorTerms}{unspecified}).
#'
#' Be aware that the dataset only contains \emph{named} craters by the IUA.
#' Therefore, there is likely a \bold{high uniform bias} on the distribution
#' of craters. Presumably the naming process attempts to cover the planet in
#' a somehow uniform fashion (distant craters are more likely to be named than
#' neighboring craters). Also, there are substantially more craters in the
#' listed bodies than those named by the IUA. See \code{\link{venus}} and
#' \code{\link{rhea}} for more detailed and specific crater datasets.
#'
#' The \eqn{(\theta, \phi)} angles are such their associated planetocentric
#' coordinates are:
#' \deqn{(\cos(\phi) \cos(\theta), \cos(\phi) \sin(\theta), \sin(\phi))',}{
#' (cos(\phi) cos(\theta), cos(\phi) sin(\theta), sin(\phi))',}
#' with \eqn{(0, 0, 1)'} denoting the north pole.
#'
#' The script performing the data preprocessing is available at
#' \href{https://github.com/egarpor/sphunif/blob/master/data-raw/craters.R}{
#' \code{craters.R}}. The data was retrieved on 2020-05-31.
#' @source \url{https://planetarynames.wr.usgs.gov/AdvancedSearch}
#' @examples
#' # Load data
#' data("craters")
#'
#' # Add Cartesian coordinates
#' craters$X <- cbind(cos(craters$theta) * cos(craters$phi),
#'                    sin(craters$theta) * cos(craters$phi),
#'                    sin(craters$phi))
#'
#' # Tests to be performed
#' type_tests <- c("PCvM", "PAD", "PRt")
#'
#' # Tests for Venus and Rhea
#' unif_test(data = craters$X[craters$target == "Venus", ], type = type_tests,
#'           p_value = "asymp")
#' unif_test(data = craters$X[craters$target == "Rhea", ], type = type_tests,
#'           p_value = "asymp")
"craters"


#' @title Venus craters
#'
#' @description Craters on Venus from the
#' \href{https://astrogeology.usgs.gov/search/map/Venus/venuscraters}{
#' USGS Astrogeology Science Center}.
#'
#' @docType data
#' @format A data frame with 967 rows and 4 variables:
#' \describe{
#'   \item{name}{name of the crater (if named).}
#'   \item{diameter}{diameter of the crater (in km).}
#'   \item{theta}{longitude angle \eqn{\theta \in [0, 2\pi)} of the
#'   crater center.}
#'   \item{phi}{latitude angle \eqn{\phi \in [-\pi/2, \pi/2]} of the
#'   crater center.}
#' }
#' @details
#' The \eqn{(\theta, \phi)} angles are such their associated planetocentric
#' coordinates are:
#' \deqn{(\cos(\phi) \cos(\theta), \cos(\phi) \sin(\theta), \sin(\phi))',}{
#' (cos(\phi) cos(\theta), cos(\phi) sin(\theta), sin(\phi))',}
#' with \eqn{(0, 0, 1)'} denoting the north pole.
#'
#' The script performing the data preprocessing is available at
#' \href{https://github.com/egarpor/sphunif/blob/master/data-raw/venus.R}{
#' \code{venus.R}}.
#' @source \url{https://astrogeology.usgs.gov/search/map/Venus/venuscraters}
#' @examples
#' # Load data
#' data("venus")
#'
#' # Add Cartesian coordinates
#' venus$X <- cbind(cos(venus$theta) * cos(venus$phi),
#'                  sin(venus$theta) * cos(venus$phi),
#'                  sin(venus$phi))
#'
#' # Tests
#' unif_test(data = venus$X, type = c("PCvM", "PAD", "PRt"), p_value = "asymp")
"venus"


#' @title Rhea craters from Hirata (2016)
#'
#' @description Craters on Rhea from
#' \href{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JE004940}{
#' Hirata (2016)}.
#'
#' @docType data
#' @format A data frame with 3596 rows and 4 variables:
#' \describe{
#'   \item{name}{name of the crater (if named).}
#'   \item{diameter}{diameter of the crater (in km).}
#'   \item{theta}{longitude angle \eqn{\theta \in [0, 2\pi)} of the
#'   crater center.}
#'   \item{phi}{latitude angle \eqn{\phi \in [-\pi/2, \pi/2]} of the
#'   crater center.}
#' }
#' @details
#' The \eqn{(\theta, \phi)} angles are such their associated planetocentric
#' coordinates are:
#' \deqn{(\cos(\phi) \cos(\theta), \cos(\phi) \sin(\theta), \sin(\phi))',}{
#' (cos(\phi) cos(\theta), cos(\phi) sin(\theta), sin(\phi))',}
#' with \eqn{(0, 0, 1)'} denoting the north pole.
#'
#' The script performing the data preprocessing is available at
#' \href{https://github.com/egarpor/sphunif/blob/master/data-raw/rhea.R}{
#' \code{rhea.R}}.
#' @source \url{https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002\%2F2015JE004940&file=jgre20485-sup-0002-TableS1.txt}
#' @references
#' Hirata, N. (2016) Differential impact cratering of Saturn's satellites by
#' heliocentric impactors. \emph{Journal of Geophysical Research: Planets},
#' 121:111--117. \doi{10.1002/2015JE004940}
#' @examples
#' # Load data
#' data("rhea")
#'
#' # Add Cartesian coordinates
#' rhea$X <- cbind(cos(rhea$theta) * cos(rhea$phi),
#'                 sin(rhea$theta) * cos(rhea$phi),
#'                 sin(rhea$phi))
#'
#' # Tests
#' unif_test(data = rhea$X[rhea$diam > 15 & rhea$diam < 20, ],
#'           type = c("PCvM", "PAD", "PRt"), p_value = "asymp")
"rhea"


#' #' @title TODO
#' #'
#' #' @description TODO.
#' #'
#' #' @docType data
#' #' @format TODO
#' #' \describe{
#' #'   TODO
#' #' }
#' #' @details
#' #' TODO
#' #' @source TODO
#' #' @references
#' #' TODO
#' #' @examples
#' #' # TODO
#' "ambrosia"
#'
#'
#' #' @rdname ambrosia
#' "angustifolia"
#'
#'
#' #' @rdname ambrosia
#' "chenopodium"
