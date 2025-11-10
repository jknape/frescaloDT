#' Analyse species occurrence data with the frescalo algorithm of Hill 2012.
#'
#' @param data Data frame with the samples. By default, first column is interpreted as the name or id of the sampled site,
#'             the second column as the observed species, and the third column the time of observation. This can be changed
#'             via the colnames argument.
#' @param weights Data frame with neighborhood weights where the first column is the target site, second column is the
#'                neighbor, and third column is the weight in the neighborhood of the target site.
#' @param phi_target Target value for adjusted frequency weighted mean frequencies. The default value, 0.74, follows the default of Hill,
#'                   but is arbitrary. If NA, target is set to the quantile of the input mean frequencies corresponding to phi_prob.
#' @param Rstar Threshold for species to be considered as benchmarks when computing time factors.
#' @param phi_prob Used to check that the phi_target is not too low. A warning is generated if thee quantile of input
#'                 mean frequencies corresponding to phi_prob is larger than phi_target.
#'                 If phi_target is set to NA, the quantile corresponding to phi_prob is taken as the target. Defaults to 0.985.
#' @param bench_exclude Vector of names of species not to be used as benchmarks when computing time factors.
#' @param col_names TODO, NOT WORKING YET. A list with elements named location, species, time, location2 and weight and values equal to the corresponding
#'                 column names in data and weight. Defaults to NULL in which case the order of the columns is used.
#'
#' @returns A frescalo object.
#'
#' @details
#' The implementation uses similar conventions to the original fortran program. E.g. small constants are added in strategic places
#' to avoid divisions by zero or other issues that can cause the algorithm to otherwise fail numerically.
#'
#'
#' @examples
#' data(bryophyte)
#' fr = frescalo(bryophyte, bryophyte_weights)
#' summary(fr)
#' @export
frescalo = function(data, weights, phi_target = .74, Rstar = 0.27, phi_prob = .985, bench_exclude = NULL, col_names = NULL) {
  samp_id = samp = spec_id = time_id = samp1_id = NULL # To avoid Notes in R CMD check

  if (is.null(col_names)) {
    data = data[,1:3]
    data_names = colnames(data)
    setDT(data) # Should copy, otherwise may be reordered on return!
    setnames(data, c("samp", "spec", "time"))

    weights = weights[,1:3]
    weight_names = colnames(weights)
    setDT(weights)
    setnames(weights, c("samp", "samp1", "wgt"))
  } else {
    expected_cols = c("location", "species", "time", "location2", "weight")
    names(col_names) = match.arg(names(col_names), expected_cols, several.ok = TRUE)
    missing_cols = setdiff(expected_cols, names(col_names))
    if (length(missing_cols)>0) {
      stop(paste0("Name of ", paste0(missing_cols, collapse = ", "), " not found in col_names."))
    }
    data_names = match.arg(colnames(data), unlist(col_names[c("location", "species", "time")]), several.ok = TRUE)
    if (length(data_names) != 3 ) {
      missing_cols = setdiff(expected_cols[1:3], names(data_names))
      stop(paste0("Column(s) ", paste0(col_names[missing_cols], collapse = ", "), " not found data."))
    }
    weight_names = match.arg(colnames(data), unlist(col_names[c("location", "location2", "weight")]), several.ok = TRUE)
    if (length(weight_names) != 3 ) {
      missing_cols = setdiff(expected_cols[c(1,4:5)], names(weight_names))
      stop(paste0("Column(s) ", paste0(col_names[missing_cols], collapse = ", "), " not found weights."))
    }
    browser()
    data = data[, c(data_names)]
  }
  stopifnot(setequal(unique(weights$samp), unique(weights$samp1))) # How handle this??
  sites = data.table(samp = sort(unique(c(weights$samp))))[, samp_id := 1:.N]

  # Argument checking
  if (!is.na(phi_target)) {
    if (phi_target > 1) {
      stop("Argument phi_target needs to be less than 1.")
    }
    if (phi_target < 0) {
      stop("Argument phi_target needs to be positive.")
    }
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

  if (length(bench_exclude) > 0) {
    bm = match(bench_exclude, species$spec, nomatch = 0)
    if (sum(bm == 0) > 0) {
      writeLines(paste("Species", paste(bench_exclude[bm == 0], collapse = ", "), "to exclude from benchmarks not found, ignored."))
    }
    bench_exclude = bench_exclude[bm > 0]
  }

  # Add keys to weights
  weights[ , samp_id := sites$samp_id[pmatch(weights$samp, sites$samp ,duplicates.ok = TRUE)]]
  weights[ , samp1_id := sites$samp_id[pmatch(weights$samp1, sites$samp ,duplicates.ok = TRUE)]]

  # Add keys to data
  data[ , spec_id := species$spec_id[match(data$spec, species$spec)]]
  data[ , samp_id := sites$samp_id[match(data$samp, sites$samp)]]
  data[ , time_id := times$time_id[match(data$time, times$time)]]

  # Compute frequency weighted mean frequencies
  freqs = nfcalc(data, weights, sites, species)
  if (is.na(phi_target)) {
    phi_target = freqs[ , .(phi0 = sum(freq^2)/sum(freq)), by = "samp_id"][, quantile(phi0, prob = phi_prob, names = FALSE)]
  }
  nc = nrow(freqs)
  set(freqs, j = c("Freq_1", "SD_Frq1", "rank", "rank1"),
      value = list(numeric(nc), numeric(nc), integer(nc), numeric(nc))) # rank1 = R´, Hill P200.
  setkey(freqs, samp_id) # Not needed, minimal improvement if any?

  freqs[, c("Freq_1", "SD_Frq1", "rank", "rank1") := frescaDT2(.SD, sites, phi_target = phi_target, irepmax = 100), keyby = list(samp_id), .SDcols = c("samp_id", "freq")]

  # Compute effort per site and time as proportion of benchmarks found.
  sampling_effort = benchmark_proportions(data, freqs, species, Rstar = Rstar, bench_exclude = bench_exclude)

  # Compute time factors.
  tfs = tfcalc(data, freqs, species, sites, times, sampling_effort)

  freqs$species = species$spec[match(freqs$spec_id, species$spec_id)]
  freqs$samp = sites$samp[match(freqs$samp_id, sites$samp_id)]
  freqs$spec_id = NULL
  freqs$samp_id = NULL
  setDF(freqs)

  tfs$species = species$spec[match(tfs$spec_id, species$spec_id)]
  tfs$spec_id = NULL
  tfs$time = times$time[match(tfs$time_id, times$time_id)]
  tfs$time_id = NULL
  setcolorder(tfs, c("species", "time", "tf", "tf_se", "n_obs", "sptot", "esttot"))
  setDF(tfs)
  out = list(call = match.call(), freqs = freqs, tfs = tfs, sites = sites, species = species, times = times,
             phi_target = phi_target, Rstar = Rstar, phi_prob = phi_prob, excluded_sites = exclude_sites,
             bench_exclude = bench_exclude, sampling_effort = sampling_effort,
             n_obs = nrow(data), n_weights = nrow(weights))
  class(out) = "frescalo"
  check_phi(out, prob = phi_prob)
  out
}

#' @exportS3Method base::summary
summary.frescalo = function(object, ...) {
  out = list(call = object$call,
             nsp = nrow(object$species),
        nsite = nrow(object$sites),
        nt = nrow(object$times),
        n_obs = object$n_obs,
        n_weights = object$n_weights,
        phi_target = object$phi_target,
        phi_in_quant = quantile(object$sites$phi_orig, probs = c(.25,.5,.75, object$phi_prob)),
        mean_nsp = c(mean(object$sites$n_spec), mean(object$sites$spnum_orig), mean(object$sites$spnum_new)))
  class(out) = "summary.frescalo"
  out
}

#' @exportS3Method base::print
print.summary.frescalo = function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\n######################################")
  cat("\n\n  Number of sites:", x$nsite)
  cat("\n  Number of species:", x$nsp)
  cat("\n  Number of time periods:", x$nt)
  cat("\n  Number of observations:", x$n_obs)
  cat("\n  Number of weights:", x$n_weights)
  cat("\n\n######################################")
  cat("\n")
  cat("\n Target phi:", x$phi_target)
  cat("\n")
  cat("\n Quantiles of input phi:\n")
  print(x$phi_in_quant, digits = 2)
  cat("\n######################################")
  cat("\n\n  Mean number of species per site")
  cat("\n  Observed:", round(x$mean_nsp[1],1))
  cat("\n  Expected, no adjustment:", round(x$mean_nsp[2],1))
  cat("\n  Expected, after adjustment:", round(x$mean_nsp[3],1))
  cat("\n\n######################################")
}

#' Extract species frequencies from a frescalo object.
#'
#' @param object An object as returned from the frescalo function.
#'
#' @returns A data frame with species frequencies across locations.
#'
#'
#' @examples
#' data(bryophyte)
#' fr = frescalo(bryophyte, bryophyte_weights)
#' frequencies(fr)
#' @export
frequencies = function(object) {
  # freqs = object$freqs # Results in additional copy of large table.
  #freqs$spec_id = NULL
  #freqs$samp_id = NULL
  #setDF(freqs)
  object$freqs
}


#' Extract time factors from a frescalo object.
#'
#' @param object An object as returned from the frescalo function.
#'
#' @returns A data frame with time factors across species.
#'
#' @examples
#' data(bryophyte)
#' fr = frescalo(bryophyte, bryophyte_weights)
#' timefactors(fr)
#' @export
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


#' Plot unscaled and rescaled neighbourhood species frequency curves.
#'
#' @param object An object as returned from the frescalo function.
#' @param max_sites Maximum number of sites for which to plot curves.
#'                  If less than the total number of sites, a random sample is taken.
#'
#' @returns A ggplot object.
#'
#' The plots correspond to Fig. 2 and 3 in Hill 2012, a comparison of neighborhoud frequency curves
#' before and after rescaling.
#' @examples
#' data(bryophyte)
#' fr = frescalo(bryophyte, bryophyte_weights)
#' check_rescaling(fr)
#' @export
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
    geom_line(alpha = 100/min(nrow(object$sites), max_sites)) + facet_wrap("scaled", scales = "free_x") + theme_minimal()

}





#' Compute estimated occupancy probabilities for a given level of effort.
#'
#' @param object An object as returned from the frescalo function.
#' @param species A character vector of species names for which to compute probabilities.
#' @param s Level of effort.
#'
#' @returns
#'
#' @note
#' Estimated probability of occurence under standard effort, sit = 1 meaning all benchmarks found (Bijlsma).
#' This depends on the proportion of benchmarks (see Prescott 2025).
#' Note: this is rather probability of detection under a sampling effort sufficient for all benchmarks to be found(?)
#  By assumption of the frescalo method trends in occupancy probability are identical across sites.
#' @examples
#' @export
occupancy_prob = function(object, species, s = 1) {

  setDT(fr$freqs)
  setDT(fr$tfs)
  keep_f = which(object$freqs$species %in% species)
  keep_tfs = which(object$tfs$species %in% species)
  out = merge(object$freqs[keep_f, list(species, samp, Freq_1)],
                          object$tfs[keep_tfs, list(species, time, tf)], allow.cartesian = TRUE)
  out[, p_occ := 1 - (1-s * Freq_1)^tf]
  out$tf = NULL
  out$Freq_1 = NULL
  setDF(fr$freqs)
  setDF(fr$tfs)
  setDF(out)
  out
}



