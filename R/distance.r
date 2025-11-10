#' Compute euclidean distances and their ranks coordinates.
#'
#' @param data A data frame with site names and x and y coordinates.
#' @param site Name of column with site/location name or id.
#' @param x Name of column containing x coordinate.
#' @param y Name of column containing y coordinate.
#' @param max_neigh The largest number of neighbours to keep. Defaults to keeping the 200 closest neighbours of each site.
#' @param max_dist Distances larger than this are ignored in the output. Not applied by default.
#'
#' @returns A data frame with distances.
#'
#' The function computes simple euclidian distances. For more accurate results,
#' use a dedicated function, such as st_dist from the sf package.
#'
#' @examples
#' @export
euclid_dist = function(data, site, x, y, max_neigh = 200, max_dist = Inf)  {
  max_neigh = min(as.integer(max_neigh), nrow(data))
  dists = apply(data, 1, function(row) {ds =  sqrt((row[[x]] - data[[x]]) ^ 2 + (row[[y]] - data[[y]]) ^ 2)
                                        if (max_dist < Inf) {ds = ds[ds < max_dist]}
                                        n = min(max_neigh, length(ds));
                                        ds = sort(ds, index.return = TRUE, method = "quick");
  list(ds$ix[1:max_neigh], ds$x[1:max_neigh], 1:max_neigh)})
  n = sapply(dists, function(l) {length(l[[1]])})
  ntot = sum(n)
  out = data.frame(site = rep(data[[site]], times = n))
  out$neigh =   do.call(c, lapply(dists, function(l) {data[[site]][l[[1]]]}))
  out$dist =   do.call(c, lapply(dists, `[[`, i = 2))
  out$dist_rank =   do.call(c, lapply(dists, `[[`, i = 3))
  out
}


similarity = function(data, site, species) {
  isDT = is.data.table(data)
  if (!isDT) {
    setDT(data)
  }
  dat = data[, list(splist = list(unique(species))), by = site]
  if (!isDT) {
    setDF(data)
  }
  sapply(dat$splist[1:2000], l = dat$splist, FUN = neigh_sorensen)
  #dists = apply(data, 1, function(row) {neigh_sorensen});
}

neigh_sorensen = function(x, l) {
  #inters = lapply(l, intersect, x = x)
  #lapply(l, match, table = x, nomatch = 0)
  inters = sapply(lapply(l, match, table = x, nomatch = 0), function(x) {sum(x>0)})
  1 - 2 * inters / (length(x) + sapply(l, length))
}


#' Convenience function for converting a distance matrix, generated e.g. by dist or vegan::dist
#' to a data.frame suitable for the frescalo function.
#'
#' @param D Distance matrix
#' @param ids Labels or ids for the sites.
#'
#' @returns description A data frame with distances between pairs of sites.
# @export
#'
#' @examples
dist2df.mat <- function(D, ids = seq_len(nrow(D))) {
  n <- nrow(D)
  a <- rep(1L:n, times = n)
  b <- rep(1L:n, each = n)
  out = data.table(site = ids[a], neigh = ids[b], dist = as.numeric(D[cbind(a, b)]))
  if (!is.null(pairs)) { # This requires first creating full n x n data.frame out. Could be avoided, but indexing gets complicated.
    setDT(pairs)
    setkeyv(pairs, colnames(pairs))
    setkey(out, "site", "neigh")
    out = out[pairs]
  }
  setorderv(out, c("site", "dist"))
  out[, dist_rank := 1:.N, by = "site"]
  if (max_neigh < n) {
    out = out[dist_rank <= max_neigh]
  }
  setDF(out)
  out
}


dist2df.dist <- function(D, labels = NULL, max_neigh = 200, pairs = NULL) {
  n <- nrow(D)
  if (is.null(labels)) {
    labels = attr(D, "Labels")
  }
  if (is.null(labels)) {
    labels = seq_len(n)
  }
  a = rep(1L:n, times = n:1)
  b = sequence(n:1L) + rep(0L:(n-1L), n:1L)
  self = seq(1, by = n+1, l = n) - cumsum(c(0,1:(n-1)))
  # Include reverse distances, larger data frame needed, but necessary to compute ranks.
  a_temp = c(a, b[-self])
  b = c(b, a[-self])
  a = a_temp
  out = data.table(site = labels[a], neigh = labels[b], dist = numeric(length(a)))
  out$dist[-self] = D # D recycled, could also use = c(D, D)
  if (!is.null(pairs)) { # This requires first creating full n x n data.frame out. Could be avoided, but indexing gets complicated.
    setDT(pairs)
    setkeyv(pairs, colnames(pairs))
    setkeyv(out, c("site", "neigh"))
    out = out[pairs]
  }
  setorderv(out, c("site", "dist"))
  out[, dist_rank := 1:.N, by = "site"]
  if (max_neigh < n) {
    out = out[dist_rank <= max_neigh]
  }
  setDF(out)
  out
}




neig_weights.data.frame = function(data) {}

neigh_weights.matrix = function(s_dist) {}


sorensen_dist = function(data, site, vars, max_neigh = 200)  {
  max_neigh = min(as.integer(max_neigh), nrow(data))
  browser()
  dists = apply(data, 1, function(row) {ds =  sorensen_dist1(unlist(row[vars]), data[vars])
                                        n = min(max_neigh, length(ds));
                                        ds = sort(ds, index.return = TRUE, method = "quick");
                                        list(ds$ix[1:max_neigh], ds$x[1:max_neigh], 1:max_neigh)})
  n = sapply(dists, function(l) {length(l[[1]])})
  ntot = sum(n)
  out = data.frame(from = rep(data[[site]], times = n))
  out$to =   do.call(c, lapply(dists, function(l) {data[[site]][l[[1]]]}))
  out$dist =   do.call(c, lapply(dists, function(l) {l[[2]]}))
  out$dist_rank =   do.call(c, lapply(dists, function(l) {l[[3]]}))
  out
}



sorensen_dist1 <- function(v, mat) {
  n <- nrow(mat); p <- ncol(mat)
  stopifnot(length(v) == p)
  if (n == 0) return(matrix(0, 0, 0))
  # row-wise sums required for denominator
  row_sums <- rowSums(mat)
  num <- matrix(0,1,n)
  for (k in seq_len(p)) {
    colk <- mat[, k]
    num = num + outer(v[k], colk, pmin)
  }
  denom <- sum(v) + rowSums(mat)
  as.numeric(1 - 2 * num / denom)
}


#sdist = sorensen_dist(grid, "gid", c("c_lat", "c_lon"), max_neigh = 100)

