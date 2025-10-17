#' Compute euclidead distances and their ranks coordinates.
#'
#' @param data A data frame with site names and x and y coordinates.
#' @param site Name of column with site/location name or id.
#' @param x Name of column containing x coordinate.
#' @param y Name of column containing y coordinate.
#' @param nneigh The largest number of neighbours to keep. Defaults to keeping the 200 closest neighbours of each site.
#' @param max_dist Distances larger than this are ignored in the output. Not applied by default.
#'
#' @returns A data frame with distances.
#'
#' The function computes simple euclidian distances. For more accurate results,
#' use a dedicated function, such as st_dist from the sf package.
#'
#' @examples
#' @export
euclid_dist = function(data, site, x, y, nneigh = 200, max_dist = Inf)  {
  nneigh = min(as.integer(nneigh), nrow(data))
  dists = apply(data, 1, function(row) {ds =  sqrt((row[[x]] - data[[x]]) ^ 2 + (row[[y]] - data[[y]]) ^ 2)
                                        if (max_dist < Inf) {ds = ds[ds < max_dist]}
                                        n = min(nneigh, length(ds));
                                        ds = sort(ds, index.return = TRUE, method = "quick"); # nneighate by rank (slower)
  list(ds$ix[1:nneigh], ds$x[1:nneigh], 1:nneigh)})
  n = sapply(dists, function(l) {length(l[[1]])})
  ntot = sum(n)
  out = data.frame(from = rep(data[[site]], times = n))
  out$to =   do.call(c, lapply(dists, function(l) {data[[site]][l[[1]]]}))
  out$dist =   do.call(c, lapply(dists, function(l) {l[[2]]}))
  out$dist_rank =   do.call(c, lapply(dists, function(l) {l[[3]]}))
  out
}
# euclid_dist = function(data, site, lat, lon, nneigh = Inf, rank = TRUE)  {
#   dat = data[c(site, lon, lat)]
#   setDT(dat)
#   n = nrow(dat)
#   out = dat[, .(to = dat[[site]], dist = sqrt((get(lon) - dat[[lon]])^2 +(get(lat) - dat[[lat]])^2)), by = site][dist < nneigh]
#   setnames(out, c("from", "to", "dist"))
#   setorder(out, from, dist, to)
#   if (rank) {
#     out[, rank:= 1:.N, by = from]
#   }
#   setDF(out)
#   out
# }

