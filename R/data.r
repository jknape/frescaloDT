#' Bryophyte data
#'
#' A subset of bryophyte observations analyzed in Hill 2012, and
#' distributed with its fortran code.
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
#' @format ## `bryophyte_weights`
#' A data frame with 35549 rows and 3 columns:
#' \describe{
#'   \item{hectad}{Id of target hectad.}
#'   \item{neighbour}{Id of neighbour hectad.}
#'   \item{weight}{Similarity weight.}
#' }
"bryophyte_weights"

