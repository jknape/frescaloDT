#' Compute euclidead distances and their ranks coordinates.
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
  out = data.frame(from = rep(data[[site]], times = n))
  out$to =   do.call(c, lapply(dists, function(l) {data[[site]][l[[1]]]}))
  out$dist =   do.call(c, lapply(dists, function(l) {l[[2]]}))
  out$dist_rank =   do.call(c, lapply(dists, function(l) {l[[3]]}))
  out
}


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

jaccard_dist0 = function(v1, v2) {
  browser()
  1 - sum(pmin(v1,v2)) / sum(pmax(v1,v2))
}


# GPT-5
jaccard_dist_tmp = function(mat) {
  n <- nrow(v2); p <- ncol(v2)
  row_sums <- rowSums(mat)
  # compute pairwise denominator = sum(pmax) = row_i + row_j - sum(pmin)
  # compute pairwise numerator = sum(pmin) using decomposition by columns
  # accumulate sum(pmin) for all pairs using crossprod of indicator-weighted values per column
  out <- numeric(0, n)
  for (k in seq_len(p)) {
    colk <- mat[, k]
    # outer(pmin) for one column can be computed as pmin(colk, colk')
    # efficient trick: sort unique values? but simplest is outer with pmin which is vectorized C
    out <- out + outer(colk, colk, pmin)
  }
  den_mat <- outer(row_sums, row_sums, "+") - num_mat
  # avoid division by zero
  D <- matrix(0, n, n)
  nz <- den_mat != 0
  D[nz] <- 1 - num_mat[nz] / den_mat[nz]
  diag(D) <- 0
  rownames(D) <- colnames(D) <- if (!is.null(rownames(mat))) rownames(mat) else paste0("r", seq_len(n))
  D
}


sorensen_dist_mat <- function(mat) {
  n <- nrow(mat); p <- ncol(mat)
  if (n == 0) return(matrix(0, 0, 0))
  # row-wise sums required for denominator
  row_sums <- rowSums(mat)
  # numerator: sum of pairwise minima across columns
  num_mat <- matrix(0, n, n)
  for (k in seq_len(p)) {
    colk <- mat[, k]
    num_mat <- num_mat + outer(colk, colk, pmin)
  }
  # Sorensen distance (also known as Bray-Curtis similarity variant):
  # distance = 1 - (2 * sum(min)) / (sum(row_i) + sum(row_j))
  denom_mat <- outer(row_sums, row_sums, "+")
  D <- matrix(0, n, n)
  nz <- denom_mat != 0
  D[nz] <- 1 - (2 * num_mat[nz]) / denom_mat[nz]
  diag(D) <- 0
  D
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


sdist = sorensen_dist(grid, "gid", c("c_lat", "c_lon"), max_neigh = 100)

