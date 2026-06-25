#' Compute euclidean distances and their ranks coordinates.
#'
#' @param crd A data frame with x and y coordinates.
#' @param labels Vector with names of sites. Should have length equal to number of columns in crd.
#' @param max_neigh The largest number of neighbours to keep. Defaults to keeping the 200 closest neighbours of each site.
#' @param max_dist Distances larger than this are omitted in the output. Not applied by default.
#'
#' @returns A data frame with distances.
#'
#' The function computes simple euclidian distances. For more accurate results,
#' use a dedicated function, such as st_dist from the sf package.
#'
#' @export
euclid_dist = function(crd, labels = NULL, max_neigh = 200, max_dist = Inf)  {
  crd = as.data.frame(crd)
  n <- nrow(crd)
  if (is.null(labels)) {
    labels = rownames(crd)
  }
  if (is.null(labels)) {
    labels = seq_len(n)
  }
  stopifnot(identical(nrow(crd), length(labels)))
  max_neigh = min(as.integer(max_neigh), nrow(crd))
  dists = apply(crd, 1, function(row) {ds =  sqrt((row[[1]] - crd[[1]]) ^ 2 + (row[[2]] - crd[[2]]) ^ 2)
    if (max_dist < Inf) {ds = ds[ds < max_dist]}
    n = min(max_neigh, length(ds));
    ds = sort(ds, index.return = TRUE, method = "quick");
    list(ds$ix[1:n], ds$x[1:n], 1:n)})
  n = sapply(dists, function(l) {length(l[[1]])})
  ntot = sum(n)
  out = data.frame(site = rep(labels, times = n))
  out$neigh =   do.call(c, lapply(dists, function(l) {labels[l[[1]]]}))
  out$dist =   do.call(c, lapply(dists, `[[`, i = 2))
  out$dist_rank =   do.call(c, lapply(dists, `[[`, i = 3))
  out
}


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



#' Compute neighbourhood weights according to Hill 2012.
#'
#' @param coords Data frame specifying locations of sites. Should have three columns:
#'               a column with site labels, a columns with x coordinates and a columns with y coordinates.
#' @param attributes Data frame containing habitat or community data. Should have a column with site labels
#'                   and columns with habitat or comunity composition.
#' @param max_sp Maximum number of spatial neighbours to consider.
#' @param max_neigh Maximum neighbourhood size. Must be smaller than max_sp.
#' @param cols A list or vector with elements named site, x, and y and values giving the corresponding
#'            column names in coords and attributes. If NULL, the order of the columns in coords and attributes
#'            are used to identify the content.
#' @returns A data frame with columns site, neigh, and the corresponding weight.
#' @export
compute_weights = function(coords, attributes, max_sp = 200, max_neigh = 100, cols = NULL) {
  coords = as.data.frame(coords)
  attributes = as.data.frame(attributes)
  if (!identical(nrow(coords), nrow(attributes))) {
    stop("coords and attributes need to have the same number of rows")
  }
  if (max_neigh > max_sp) {
    stop("max_neigh must be less than max_sp")
  }
  if (!is.null(cols)) {
    crd = coords[, match_cols2(coords, cols, c("x", "y"))]
    site_labs = coords[, match_cols2(coords, cols, c("site"))]
    site_col_att = match_cols2(attributes, cols, c("site"))
    site_labs_a = attributes[[site_col_att]]
    attrib = attributes[, -(site_col_att)]
  } else {
    crd = coords[, 2:3]
    site_labs = coords[, 1]
    site_col_att = 1
    site_labs_a = attributes[[site_col_att]]
    attrib = attributes[, -(site_col_att)]
  }
  if (!all(sapply(crd, is.numeric))) {
    stop("Coordinates are not numeric.")
  }
  if (!all(sapply(attrib, is.numeric))) {
    stop("Some attributes are not numeric")
  }
  if (!identical(site_labs, site_labs_a)) {
    ord = match(site_labs, site_labs_a, nomatch = 0)
    if (!identical(site_labs, site_labs_a[ord])) {
      stop("Site labels in coords not matching labels in attributes.")
    }
    attrib = attrib[ord,]
  }
  calc_weights(crd, attrib, site_labs, max_sp, max_neigh)
}


# Assumes coords, attributes and labels have the same ordering.
calc_weights = function(coords, attributes, labels, max_dist, max_neigh) {
  dist = weight = dist_rank = hab_rank = NULL
  #euclid_dist = eco_dist(coords, labels = labels, max_neigh = max_dist, method = "euclidean")
  spat_dist = euclid_dist(coords, labels = labels, max_neigh = max_dist)
  hab_dist = eco_dist(attributes, labels = labels, max_neigh = max_neigh, subset = spat_dist[,c("site", "neigh")], method = "bray")
  setDT(spat_dist)
  spat_dist[, dist := NULL]
  setDT(hab_dist)
  hab_dist[, dist := NULL]
  setnames(hab_dist, "dist_rank", "hab_rank")
  weights = spat_dist[hab_dist, on = c("site", "neigh")][, weight := (1 - (dist_rank - 1)^2 / max_dist^2)^4 * (1 - (hab_rank - 1)^2 / max_neigh^2)^4]
  weights[, c("dist_rank", "hab_rank") := NULL]
  setDF(weights)
  weights
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
    pairs = data.table::as.data.table(pairs)
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
#' @param subset A two column matrix or data frame, or a list with two vectors, that specify a subset of pairs of sites for which to compute distances.
#'                Elements need to match the labels.
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
#' eco_dist(X, subset = list(c(1,1,1), 2:4))
eco_dist <- function(x, labels = NULL, max_neigh = 200, subset = NULL, method = "bray") {
  n <- nrow(x)
  if (is.null(labels)) {
    labels = rownames(x)
  }
  if (is.null(labels)) {
    labels = seq_len(n)
  }
  stopifnot(identical(nrow(x), length(labels)))
  if (is.null(subset)) {
    out = data.frame(site = rep(labels, each = length(labels)), neigh = rep(labels, length(labels)))
  } else {
    if (is.matrix(subset)) {
      out = data.frame(site = subset[,1], neigh = subset[,2])
    } else {
      out = data.frame(site = subset[[1]], neigh = subset[[2]])
    }
  }
  ind_mat = cbind(match(out$site, labels), match(out$neigh, labels))
  if (any(is.na(ind_mat))) {
    stop("All elements in subset not found in labels.")
  }
  out$dist = vegd_p(as.matrix(x), ind_mat, method)
  out = data.table::as.data.table(out) # Cannot use setDT here, may cause global modifications.
  setorderv(out, c("site", "dist", "neigh"))
  out[, "dist_rank" := 1:.N, by = "site"]
  if (max_neigh < n) {
    dist_rank = NULL
    out = out[dist_rank <= max_neigh]
  }
  setDF(out)
  out
}





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
