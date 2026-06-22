# #' Compute euclidean distances and their ranks coordinates.
# #'
# #' @param data A data frame with site names and x and y coordinates.
# #' @param site Name of column with site/location name or id.
# #' @param x Name of column containing x coordinate.
# #' @param y Name of column containing y coordinate.
# #' @param max_neigh The largest number of neighbours to keep. Defaults to keeping the 200 closest neighbours of each site.
# #' @param max_dist Distances larger than this are omitted in the output. Not applied by default.
# #'
# #' @returns A data frame with distances.
# #'
# #' The function computes simple euclidian distances. For more accurate results,
# #' use a dedicated function, such as st_dist from the sf package.
# #'
# #' @export
# euclid_dist = function(data, site, x, y, max_neigh = 200, max_dist = Inf)  {
#   max_neigh = min(as.integer(max_neigh), nrow(data))
#   dists = apply(data[,c(x,y)], 1, function(row) {ds =  sqrt((row[[x]] - data[[x]]) ^ 2 + (row[[y]] - data[[y]]) ^ 2)
#                                         if (max_dist < Inf) {ds = ds[ds < max_dist]}
#                                         n = min(max_neigh, length(ds));
#                                         ds = sort(ds, index.return = TRUE, method = "quick");
#                                         list(ds$ix[1:max_neigh], ds$x[1:max_neigh], 1:max_neigh)})
#   n = sapply(dists, function(l) {length(l[[1]])})
#   ntot = sum(n)
#   out = data.frame(site = rep(data[[site]], times = n))
#   out$neigh =   do.call(c, lapply(dists, function(l) {data[[site]][l[[1]]]}))
#   out$dist =   do.call(c, lapply(dists, `[[`, i = 2))
#   out$dist_rank =   do.call(c, lapply(dists, `[[`, i = 3))
#   out
# }

#
# similarity = function(data, site, species) {
#   isDT = is.data.table(data)
#   if (!isDT) {
#     setDT(data)
#   }
#   dat = data[, list(splist = list(unique(species))), by = site]
#   if (!isDT) {
#     setDF(data)
#   }
#   sapply(dat$splist[1:2000], l = dat$splist, FUN = neigh_sorensen)
#   #dists = apply(data, 1, function(row) {neigh_sorensen});
# }
#
# neigh_sorensen = function(x, l) {
#   #inters = lapply(l, intersect, x = x)
#   #lapply(l, match, table = x, nomatch = 0)
#   inters = sapply(lapply(l, match, table = x, nomatch = 0), function(x) {sum(x>0)})
#   1 - 2 * inters / (length(x) + sapply(l, length))
# }


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


#' Compute ecological or habitat similarity between pairs of sites.
#'
#' @param x A numeric matrix with positive entries representing for example abundances of different species (columns) or amount or proportions of habitat classes.
#' @param labels A vector of names of the sites. Should have length equal to the number of rows in x. If NULL, the row names of x are taken as labels.
#' @param max_neigh The maximum number of of neighbours to keep from each site. The closest neighbours are kept.
#' @param include A two column matrix or data frame, or a list with two vectors, that specify a subset of pairs of sites for which to compute distances.
#'                Elements need to match the labels.
#'                Useful when there are many sites which would lead to very large complete distance matrices.
#' @param method Methods used to compute distances. A subset of the methods available for vegdist in package vegan can be used.
#'
#' @returns A data frame with pairwise distances and their within site ranks.
#'
#' @note
#' Underlying C code for distance computations is based on code from the vegdist function in the vegan package.
#' The original code is modified mainly to avoid the need to compute full distance matrices corresponding to all pairs of sites.
#' @export
#' @examples
#' # Simulate distances
#' X = exp(matrix(rnorm(30), ncol = 5))
#' # Keep the two neighbours closest to each site:
#' eco_dist(X, max_neigh = 2)
#' # Compute distances from site 1 to sites 2, 3 and 4:
#' eco_dist(X, include = list(c(1,1,1), 2:4))
eco_dist <- function(x, labels = NULL, max_neigh = 200, include = NULL, method = "bray") {
  n <- nrow(x)
  if (is.null(labels)) {
    labels = rownames(x)
  }
  if (is.null(labels)) {
    labels = seq_len(n)
  }
  stopifnot(identical(nrow(x), length(labels)))
  if (is.null(include)) {
    out = data.frame(site = rep(labels, each = length(labels)), neigh = rep(labels, length(labels)))
  } else {
    if (is.matrix(include)) {
      out = data.frame(site = include[,1], neigh = include[,2])
    } else {
      out = data.frame(site = include[[1]], neigh = include[[2]])
    }
  }
  ind_mat = cbind(match(out$site, labels), match(out$neigh, labels))
  if (any(is.na(ind_mat))) {
    stop("All elements in include not found in labels.")
  }
  out$dist = vegd_p(as.matrix(x), ind_mat, method)
  setDT(out)
  setorderv(out, c("site", "dist", "neigh"))
  out[, "dist_rank" := 1:.N, by = "site"]
  if (max_neigh < n) {
    dist_rank = NULL
    out = out[dist_rank <= max_neigh]
  }
  setDF(out)
  out
}




neig_weights.data.frame = function(data) {}

neigh_weights.matrix = function(s_dist) {}

#' @useDynLib frescaloDT, .registration = TRUE
NULL


vegd_p = function(x, pairs, method = "bray") {
  x = as.matrix(x)
  nr = nrow(x)
  method = get_dist_method(x, method)
  stopifnot(ncol(pairs) == 2)
  pairs = as.integer(pairs)
  stopifnot(all(pairs > 0) & all(pairs <= nr))
  d = .Call("do_vegdist_p", x, pairs, as.integer(method))
  d[d < .Machine$double.eps] = 0
  if (any(is.na(d)))
    warning("missing values in results")
  d
}



# Adapted from vegan::vegdist
get_dist_method  = function (x, method = "bray", na.rm = FALSE) {
    if (!is.na(pmatch(method, "euclidian")))
      method <- "euclidean"
    ## the order of METHODS below *MUST* match the #define'd numbers
    ## in vegdist.c
    METHODS <- c("manhattan", "euclidean", "canberra", "bray", # 4
                 "kulczynski", "gower", "morisita", "horn", #8
                 "mountford", "jaccard", "raup", "binomial", "chao", #13
                 "altGower", "cao", "mahalanobis", "clark", "chisq", "chord", #19
                 "hellinger", "aitchison", "robust.aitchison") # 22
    method <- pmatch(method, METHODS)
    if (method %in% c(6, 16, 18, 21, 22)) {
      stop(paste("method", METHODS[method], "not available."))
    }
    inm <- METHODS[method]
    if (is.na(method))
      stop("invalid distance method")
    if (method == -1)
      stop("ambiguous distance method")
    if (!na.rm && anyNA(x))
      stop("missing values are not allowed with argument 'na.rm = FALSE'")
    ## all vegdist indices need numeric data (Gower included).
    if (!(is.numeric(x) || is.logical(x)))
      stop("input data must be numeric")

    if (!method %in% c(1,2,6,16,18) && any(rowSums(x, na.rm = TRUE) == 0))
      warning("you have empty rows: their dissimilarities may be
                 meaningless in method ",
              dQuote(inm))
    ## 1 manhattan, 2 euclidean, 3 canberra, 6 gower, 16 mahalanobis, 19 chord
    if (!method %in% c(1,2,3,6,16,19,20) && any(x < 0, na.rm = TRUE))
      warning("results may be meaningless because data have negative entries
                 in method ",
              dQuote(inm))
    if (method %in% c(11,18) && any(colSums(x) == 0, na.rm = TRUE))
      warning("data have empty species which influence the results in
                 method ",
              dQuote(inm))
    #if (method == 6) # gower, but no altGower
    #  x <- decostand(x, "range", 2, na.rm = TRUE, ...)
    #if (method == 16) # mahalanobis
    #  x <- veganMahatrans(scale(x, scale = FALSE), na.rm = na.rm)
    #if (method == 18) # chisq
    #  x <- decostand(x, "chi.square", na.rm = na.rm)
    #if (method == 21)  # aitchison
    #  x <- decostand(x, "clr", ...)  # dots to pass possible pseudocount
    #if (method == 22)  # robust.aitchison
    #  x <- decostand(x, "rclr", na.rm = na.rm, ...) # No pseudocount for rclr
    N <- nrow(x)
    if (method %in% c(7, 13, 15) && !isTRUE(all.equal(x, round(x))))
      warning("results may be meaningless with non-integer data in method ",
              dQuote(inm))
    if (method %in% 7) { # morisita
      if (any(x[x>0] < 1))
        warning("results may be meaningless with positive values < 1 in ",
                dQuote(inm))
      if (round(max(x)) == 1)
        warning("results may be meaningless with largest integer 1 in ",
                dQuote(inm))
      else if (any(r1 <- apply(x, 1, max) <= 1))
        warning(dQuote(inm),
                " expects some counts > 1, none in row(s) ",
                paste(which(r1), collapse=", "))
    }
    method
}
