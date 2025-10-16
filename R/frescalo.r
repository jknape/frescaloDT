#' Analyse species occurrence data with the frescalo algorithm of Hill 2012.
#'
#' @param data Data frame with the samples. By default, first column is interpreted as the name or id of the sampled site,
#'             the second column as the observed species, and the third column the time of observation. This can be changed
#'             via tha colnames argument.
#' @param weights Data frame with neighbourhood weights where the first column is the target site, second column is the
#'                neighbour, and third column is the weight in the neigbourhood of the target site.
#' @param phi_target Target value for adjusted frequency weighted mean frequencies. The default value, 0.74, follows the default of Hill,
#'                   but is arbitrary.
#' @param Rstar Threshold for species to be considered as benchmarks when computing time factors.
#' @param bench_exclude Vector of names of species not to be used as benchmarks when computing time factors.
#' @param colnames A list with elements named location, species, time, location2 and weigth and values equal to the corresponding
#'                 column names in data and weight. Defaults to NULL in which the order of the columns is used.
#'
#' @returns An object of class frescalo
#' @export
#'
#' The implentation uses similar conventions to the original fortran program. E.g. small constants are added in strategic places
#' to avoid divisions by zero or other issues that can cause the algorithm to otherwise fail.
#'
#'
#' @examples
frescalo = function(data, weights, phi_target = .74, Rstar = 0.27, bench_exclude = NULL, colnames = NULL) {
  samp_id = samp = spec_id = time_id = samp1_id = NULL # To avoid Notes in R CMD check

  if (is.null(colnames)) {
    data = data[,1:3]
    data_names = colnames(data)
    setDT(data) # Should copy, otherwise may be reordered on return!
    setnames(data, c("samp", "spec", "time"))

    weights = weights[,1:3]
    weight_names = colnames(weights)
    setDT(weights)
    setnames(weights, c("samp", "samp1", "wgt"))
  }
  stopifnot(setequal(unique(weights$samp), unique(weights$samp1))) # How handle this??
  sites = data.table(samp = sort(unique(c(weights$samp))))[, samp_id := 1:.N]

  # Argument checking
  if (phi_target > 1) {
    stop("Argument phi_target needs to be less than 1.")
  }
  if (phi_target < 0) {
    stop("Argument phi_target needs to be positive.")
  }

  if (Rstar < 0) {
     stop("Argument Rstar should be positive.")
  }

  # Filter sites in data not present in weights
  exclude_sites = setdiff(unique(data$samp), sites$samp)
  if (length(exclude_sites) > 0) {
    message(paste("Site(s)", paste(exclude_sites, collapse = ", "), "not present in weights, removed."))
  }
  data =  data[!(samp %in% exclude_sites)]

  # Key tables for species and times
  species = data.table(spec = sort(unique(data$spec)))[, spec_id := 1:.N] # Note: species may have been removed, if only present in excluded sites.
  times = data.table(time = sort(unique(data$time)))[, time_id := 1:.N]

  # Add keys to weights
  weights[ , samp_id := sites$samp_id[pmatch(weights$samp, sites$samp ,duplicates.ok = TRUE)]]
  weights[ , samp1_id := sites$samp_id[pmatch(weights$samp1, sites$samp ,duplicates.ok = TRUE)]]

  # Add keys to data
  data[ , spec_id := species$spec_id[match(data$spec, species$spec)]]
  data[ , samp_id := sites$samp_id[match(data$samp, sites$samp)]]
  data[ , time_id := times$time_id[match(data$time, times$time)]]

  # Compute frequency weighted mean frequencies
  freqs = nfcalc(data, weights, sites, species)
  nc = nrow(freqs)
  set(freqs, j = c("Freq_1", "SD_Frq1", "rank", "rank1") ,
      value = list(numeric(nc), numeric(nc), integer(nc), numeric(nc))) # rank1 = R´, Hill P200.
  setkey(freqs, samp_id) # Not needed, minimal improvement if any?

  freqs[, c("Freq_1", "SD_Frq1", "rank", "rank1") := frescaDT2(.SD, sites, phi_target = phi_target, irepmax = 100), keyby = list(samp_id), .SDcols = c("samp_id", "freq")]

  tfs = tfcalc(data, freqs, species, sites, times, Rstar = Rstar, no_bench = bench_exclude)

  freqs$species = species$spec[match(freqs$spec_id, species$spec_id)]
  freqs$samp = sites$samp[match(freqs$samp_id, sites$samp_id)]
  freqs$spec_id = NULL
  freqs$samp_id = NULL
  setDF(freqs)

  tfs$species = species$spec[match(tfs$spec_id, species$spec_id)]
  tfs$spec_id = NULL
  tfs$time = times$time[match(tfs$time_id, times$time_id)]
  tfs$time_id = NULL
  setcolorder(tfs, c("species", "time", "tf", "se", "jtot", "sptot", "esttot"))
  setDF(tfs)
  out = list(freqs = freqs, tfs = tfs, sites = sites, species = species, times = times,
             excluded_sites = exclude_sites, phi_target = phi_target, Rstar = Rstar)
  class(out) = "frescalo"
  check_phi(out)
  out
}

frequencies = function(object) {
  # freqs = object$freqs # Results in additional copy of large table.
  #freqs$spec_id = NULL
  #freqs$samp_id = NULL
  #setDF(freqs)
  object$freqs
}

timefactors = function(object) {
  #if (is.null(object$tfs)) {
    #tfalc....
  #}
  object$tfs
}

check_phi = function(object, prob = .985, plot = FALSE) {
  obs_p = stats::quantile(object$sites$phi_orig, probs = prob)
  phi_ok = obs_p < object$phi_target
  message(paste0(round(100 * prob, 1), " percentile of input phi = ", round(obs_p,2), "\n",
                "Target phi = ", round(object$phi_target, 2)))
  if (!phi_ok) {
    warning("Target value of phi may be too small.")
  }
  if (plot) {
    graphics::hist(object$sites$phi_orig, xlim = c(0, 1), xlab = "phi", main = "")
    graphics::abline(v = object$phi_target, col = "red")
  }
  if (!all(object$sites$conv.fresca)) {
    conv_fail = !object$sites$conv.fresca
    warning(paste("phi did not converge to phi_target for all sites, increasing max.iter might help. Convergence failed for site(s):",
                  paste(object$sites$samp[conv_fail], collapse = ", ")))
  }
  phi_ok
}


# Is this meaningful?
check_r = function(object) {
  graphics::hist(object$freqs$rank1, xlab = "R", main = "")
  graphics::abline(v = object$Rstar, col = "red")
}

#check_rescaling = function(object) {
#  op = par()
#  par(mfcol = c(1,2))
#  graphics::hist(object$freqs$rank1, xlab = "R", main = "")
#  graphics::abline(v = object$Rstar, col = "red")
#}


# ~ Fig 2 & 3 in Hill
check_rescaling = function(object, max_sites = 500) {
  scaled = NULL
  if (nrow(object$sites) > max_sites) {
    keep = object$freqs[["samp"]] %in% sample(object$sites[["samp"]], max_sites)
  } else {
    keep = TRUE
  }
  pldat = object$freqs[keep, c("samp", "rank", "rank1", "freq", "Freq_1")]
  setDT(pldat)
  pldat[, rank := as.numeric(rank)]
  pldat = data.table::melt(pldat,
                           id.vars = "samp", measure.vars = list(freq = c("freq", "Freq_1"), rank = c("rank", "rank1")),
                           variable.name = "scaled", variable.factor = FALSE)
  pldat[, scaled := c("unscaled", "rescaled")[as.integer(scaled)]]
  ggplot(data = pldat, aes(x = .data[["rank"]], y = .data[["freq"]], group = .data[["samp"]])) +
    geom_line(alpha = 100/min(nrow(object$sites), max_sites)) + facet_wrap("scaled", scales = "free_x")

}

# Estimated probability of occurence under standard effort, sit = 1 meaning all benchmarks found (Bijlsma). Note that this depends on the proportion of benchmarks (see Prescott 2025).
# Note: this is rather probability of detection under a sampling effort sufficient for all benchmarks to be found(?)
# Also assumes identical trends across all sites...
# Computed for a subset of species to avoid huge output.
prob_occ = function(object, spec, s = 1) {
  setDT(fr$freqs)
  setDT(fr$tfs)
  out = merge(object$freqs[species %in% spec, list(species, samp, Freq_1)],
                          object$tfs[species %in% spec, list(species, time, tf)], allow.cartesian = TRUE)
  out[, p_occ := 1 - (1-s * Freq_1)^tf]
  out$tf = NULL
  out$Freq_1 = NULL
  setDF(fr$freqs)
  setDF(fr$tfs)
  setDF(out)
  out
}



