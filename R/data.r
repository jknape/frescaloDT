#' Bryophyte data
#'
#' A subset of bryophyte observations analyzed in Hill 2012, and
#' distributed with the fortran frescalo code.
#'
#' @name bryophyte
#' @format ## `bryophyte`
#' A data frame with 78332 rows and 3 columns:
#' \describe{
#'   \item{hectad}{Id of hectad (site) where species was observed.}
#'   \item{species}{Id of species.}
#'   \item{year}{Year of observation.}
#' }
#' @source <https://www.brc.ac.uk/biblio/frescalo-computer-program-analyse-your-biological-records>
#' @references Hill, M. O. 2012. Local frequency as a key to interpreting species occurrence data when recording effort is not known.
"bryophyte"


#' @rdname bryophyte
#' @format ## `hectad_locations`
#' A data frame with 404 rows and 3 columns:
#' \describe{
#'   \item{hectad}{Id of target hectad.}
#'   \item{x}{x coordinate of hectad.}
#'   \item{y}{y coordinate of hectad.}
#' }
"hectad_locations"

#' @rdname bryophyte
#' @format ## `vascular_plants`
#' A data frame with 404 rows and 1192 columns:
#' \describe{
#'   \item{hectad}{Id of target hectad.}
#'   \item{other columns}{Presence of absence of vascular plant species in hectad. One column for each vascular plant species.}
#' }
"vascular_plants"


