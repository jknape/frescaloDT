#' Title
#'
#' @param data A data frame with sitenames and x and y coordinates.
#' @param site Name of column with site/location name or id.
#' @param lat Name of column containig y coordinate.
#' @param lon Name of column containig x coordinate.
#' @param trunc Only distances smaller than this will be kept in the output. Useful for keeping the size of the output down when there are many sites.
#' @param rank If TRUE, the ranks of the distances from each site are computed.
#'
#' @returns A data frame with distances.
#' @export
#'
#' The function computes simple euclidian distances. For more accurate results,
#' use a dedicated function, such as st_dist from the sf package.
#'
#' @examples
euclid_dist = function(data, site, lat, lon, trunc = Inf, rank = TRUE)  {
  dat = data[c(site, lon, lat)]
  setDT(dat)
  n = nrow(dat)
  out = dat[, .(to = dat[[site]], dist = sqrt((get(lon) - dat[[lon]])^2 +(get(lat) - dat[[lat]])^2)), by = site][dist < trunc]
  setnames(out, c("from", "to", "dist"))
  setorder(out, from, dist, to)
  if (rank) {
    out[, rank:= 1:.N, by = from]
  }
  setDF(out)
  out
}

