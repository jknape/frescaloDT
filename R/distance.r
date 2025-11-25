#' Compute euclidean distances and their ranks coordinates.
#'
#' @param data A data frame with site names and x and y coordinates.
#' @param site Name of column with site/location name or id.
#' @param x Name of column containing x coordinate.
#' @param y Name of column containing y coordinate.
#' @param max_neigh The largest number of neighbours to keep. Defaults to keeping the 200 closest neighbours of each site.
#' @param max_dist Distances larger than this are omitted in the output. Not applied by default.
#'
#' @returns A data frame with distances.
#'
#' The function computes simple euclidian distances. For more accurate results,
#' use a dedicated function, such as st_dist from the sf package.
#'
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


#' Convert a distance matrix, generated e.g. by dist or vegan::dist
#' to a data.frame suitable for the frescalo function. EXPERIMENTAL.
#'
#' @param D Distance matrix
#' @param ids Labels or ids for the sites.
#' @param pairs A data.frame with pairs of site for which distance should be computed.
#'
#' @returns description A data frame with distances between pairs of sites.
#' @export
dist2df.mat <- function(D, ids = seq_len(nrow(D)), pairs = NULL) {
  max_neigh = dist_rank = NULL
  n <- nrow(D)
  a <- rep(1L:n, times = n)
  b <- rep(1L:n, each = n)
  out = data.table(site = ids[a], neigh = ids[b], dist = as.numeric(D[cbind(a, b)]))
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


eco_dist <- function(x, labels, max_neigh = 200, pairs = NULL, method = 4) {
  n <- nrow(x)
  if (is.null(labels)) {
    labels = attr(D, "Labels")
  }
  if (is.null(labels)) {
    labels = seq_len(n)
  }
  stopifnot(identical(nrow(x), length(labels)))
  if (is.null(pairs)) {
    out = data.frame(site = rep(labels, each = length(labels)), neigh = rep(labels, length(labels)))
  } else {
    out = data.frame(site = pairs[,1], neigh = pairs[,2])
  }
  ind_mat = cbind(match(out$site, labels), match(out$neigh, labels))
  if (any(is.na(ind_mat))) {
    stop("All elements in pairs not found in labels.")
  }
  out$dist = vegd_test_p(as.matrix(x), ind_mat, method)
  setDT(out)
  setorderv(out, c("site", "dist", "neigh"))
  out[, "dist_rank" := 1:.N, by = "site"]
  if (max_neigh < n) {
    out = out[dist_rank <= max_neigh]
  }
  setDF(out)
  out
}




neig_weights.data.frame = function(data) {}

neigh_weights.matrix = function(s_dist) {}

#' @useDynLib frescaloDT, .registration = TRUE
NULL


#'@export
vegd_test_p = function(x, pairs, method = 4) {
  nr = nrow(x)
  stopifnot(ncol(pairs) == 2)
  pairs = as.integer(pairs)
  stopifnot(all(pairs > 0) & all(pairs <= nr))
  d = .Call("do_vegdist_p", x, pairs, as.integer(method))
  d[d < 1e-15] = 0
  d
}

#sdist = sorensen_dist(grid, "gid", c("c_lat", "c_lon"), max_neigh = 100)

