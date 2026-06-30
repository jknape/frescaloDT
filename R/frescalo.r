#' Analyse species occurrence data with the frescalo algorithm of Hill 2012.
#'
#' @param data Data frame with the samples. By default, first column is interpreted as the name or id of the sampled site,
#'             the second column as the observed species, and the third column the time of observation. This can be changed
#'             via the \code{cols} argument.
#' @param weights Data frame with neighborhood weights where the first column is the target site, second column is the
#'                neighbor, and third column is the weight in the neighborhood of the target site.
#' @param phi_target Target value for adjusted frequency weighted mean frequencies. The default value, 0.74, follows the default of Hill,
#'                   but is arbitrary. If NA, target is set to the quantile of the input mean frequencies corresponding to \code{phi_prob}.
#' @param Rstar Threshold for species to be considered as benchmarks when computing time factors.
#' @param phi_prob Used to check that \code{phi_target} is not too low. A warning is generated if thee quantile of input
#'                 mean frequencies corresponding to \code{phi_prob} is larger than \code{phi_target}.
#'                 If \code{phi_target} is set to NA, the quantile corresponding to \code{phi_prob} is taken as the target.
#'                 Defaults to 0.985.
#' @param bench_exclude Vector of names of species not to be used as benchmarks when computing time factors.
#' @param cols  A named list or character vector with elements named site, spec, time,neigh, weight,
#'                   and values equal to the corresponding column names in data and weights.
#'                   Defaults to NULL in which case the order of the columns is used.
#'
#' @returns A frescalo object.
#'
#' @details
#' The implementation uses similar conventions to the original fortran program. E.g. small constants are added in strategic places
#' to avoid divisions by zero or other issues that can cause the algorithm to otherwise fail numerically.
#'
#'
#' @examples
#' weights = compute_weights(hectad_locations, vascular_plants)
#' fr = frescalo(bryophyte, weights)
#' summary(fr)
#' head(frequencies(fr))
#' head(timefactors(fr))
#' @export
frescalo = function(data, weights, phi_target = .74, Rstar = 0.27, phi_prob = .985, bench_exclude = NULL,
                    cols = NULL) {
  site_id = site = spec_id = time_id = neigh_id = phi0 = freq_obs = NULL # To avoid Notes in R CMD check

  if (is.null(cols)) {
    if (ncol(data) < 3 | ncol(weights) < 3) {
      stop("Arguments data and weights need to have three columns.")
    }
    data = data[,1:3]
    cols = c(colnames(data), colnames(weights)[2:3])
    names(cols) = c("site", "species", "time", "neigh", "weight")
    data = data.table::as.data.table(data) # Should copy, otherwise may be reordered on return!
    setnames(data, c("site", "species", "time"))
    weights = weights[,1:3]
    weights = data.table::as.data.table(weights)
    setnames(weights, c("site", "neigh", "weight"))
  } else {
    cind = match_cols2(data, cols, c("site", "species", "time"))
    data = data.table::as.data.table(data[cind])
    setnames(data, c("site", "species", "time"))
    cind = match_cols2(weights, cols, c("site", "neigh", "weight"))
    weights = data.table::as.data.table(weights[cind])
    setnames(weights, c("site", "neigh", "weight"))
  }
  data = data[!duplicated(data)]

  stopifnot(setequal(unique(weights$site), unique(weights$neigh))) # How handle this??
  sites = data.table(site = sort(unique(c(weights$site))))[, site_id := 1:.N]

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
  exclude_sites = setdiff(unique(data$site), sites$site)
  if (length(exclude_sites) > 0) {
    message(paste("Site(s)", paste(exclude_sites, collapse = ", "), "not present in weights, removed."))
  }
  data =  data[!(site %in% exclude_sites)]

  # Key tables for species and times
  species = data.table(species = sort(unique(data$species)))[, spec_id := 1:.N] # Note: species may have been removed, if only present in excluded sites.
  times = data.table(time = sort(unique(data$time)))[, time_id := 1:.N]

  if (length(bench_exclude) > 0) {
    bm = match(bench_exclude, species$species, nomatch = 0)
    if (sum(bm == 0) > 0) {
      writeLines(paste("Species", paste(bench_exclude[bm == 0], collapse = ", "), "to exclude from benchmarks not found, ignored."))
    }
    bench_exclude = bench_exclude[bm > 0]
  }

  # Add keys to weights
  weights[ , site_id := sites$site_id[pmatch(weights$site, sites$site ,duplicates.ok = TRUE)]]
  weights[ , neigh_id := sites$site_id[pmatch(weights$neigh, sites$site ,duplicates.ok = TRUE)]]

  # Add keys to data
  sp_id = species$spec_id[match(data$species, species$species)]
  data[ , spec_id := sp_id]
  data[ , site_id := sites$site_id[match(data$site, sites$site)]]
  data[ , time_id := times$time_id[match(data$time, times$time)]]

  # Compute frequency weighted mean frequencies
  freqs = nfcalc(data, weights, sites, species)
  if (is.na(phi_target)) {
    phi_target = freqs[ , list(phi0 = sum(freq_obs^2)/sum(freq_obs)), by = "site_id"][, quantile(phi0, prob = phi_prob, names = FALSE)]
  }
  nc = nrow(freqs)
  set(freqs, j = c("freq_est", "freq_est_se", "rank", "rank_scaled"),
      value = list(numeric(nc), numeric(nc), integer(nc), numeric(nc))) # rank_scaled = R´, Hill P200.
  setkey(freqs, site_id) # Not needed, minimal improvement if any?

  freqs[, c("freq_est", "freq_est_se", "rank", "rank_scaled") := frescaDT(.SD, sites, phi_target = phi_target, irepmax = 100), keyby = list(site_id), .SDcols = c("site_id", "freq_obs")]

  # Compute effort per site and time as proportion of benchmarks found.
  sampling_effort = benchmark_proportions0(data, freqs, species, Rstar = Rstar, bench_exclude = bench_exclude)

  # Compute time factors.
  tfs = tfcalc(data, freqs, species, sites, times, sampling_effort)

  freqs$species = species$species[match(freqs$spec_id, species$spec_id)]
  freqs$site = sites$site[match(freqs$site_id, sites$site_id)]
  setcolorder(freqs, c("site", "species",  "observed", "freq_obs", "freq_est", "freq_est_se", "rank", "rank_scaled"))
  setDF(freqs) # Copy on modify avoids risk of changing frescalo object via frequencies() etc

  tfs$species = species$species[match(tfs$spec_id, species$spec_id)]
  tfs$spec_id = NULL
  tfs$time = times$time[match(tfs$time_id, times$time_id)]
  tfs$time_id = NULL
  setcolorder(tfs, c("species", "time", "tf", "tf_se", "count", "occ_adj", "occ_est"))
  setDF(tfs)

  setcolorder(sites, c("site","species_count", "phi_obs", "alpha", "spnum_obs", "spnum_est", "wgt_n2", "phi_out", "iter_fresca", "conv_fresca"))
  setDF(sites)
  setDF(species)
  setDF(times)
  setDF(sampling_effort)
  out = list(call = match.call(),
             cols = cols,
             freqs = freqs, tfs = tfs, sites = sites, species = species, times = times,
             phi_target = phi_target, Rstar = Rstar, phi_prob = phi_prob, excluded_sites = exclude_sites,
             bench_exclude = bench_exclude, sampling_effort = sampling_effort,
             n_obs = nrow(data), n_weights = nrow(weights))
  class(out) = "frescalo"
  check_phi(out, prob = phi_prob, silent = TRUE)
  out
}

match_cols = function(data, col_names, expected) {
  parg = deparse1(substitute(col_names))
  if (length(col_names) != length(expected)) {
    stop(paste(parg, "should have length", length(expected)))
  }
  col_names = unlist(col_names)
  list_names = names(col_names)
  list_names = match.arg(list_names, expected, several.ok = TRUE)
  missing_cols = setdiff(expected, list_names)
  if (length(missing_cols)>0) {
    stop(paste0(paste0(missing_cols, collapse = ", "), " not found in ", parg))
  }
  data_names = match.arg(colnames(data), col_names, several.ok = TRUE)
  if (length(data_names) != length(expected)) {
    missing_cols = setdiff(expected, names(data_names))
    stop(paste0("Column(s) ", paste0(col_names[list_names %in% missing_cols], collapse = ", "), " not in found data."))
  }
  data_names
}

match_cols2 = function(data, col_names, required) {
  parg = deparse1(substitute(col_names))
  col_names = unlist(col_names)
  list_names = names(col_names)
  missing_names = setdiff(required, list_names)
  if (length(missing_names)>0) {
    stop(paste0(paste0(missing_names, collapse = ", "), " not found in ", parg))
  }
  data_names = colnames(data)
  missing_cols = setdiff(col_names[required], data_names)
  if (length(missing_cols) > 0) {
    stop(paste0("Column(s) ", paste0(missing_cols, collapse = ", ")," listed in ", parg," not found in ", deparse1(substitute(data))))
  }
  match(col_names[match(required, names(col_names))], data_names)
}

#' Summarize frescalo output
#'
#' @param object A frescalo object.
#'
#' @param ... Other arguments.
#'
#' @exportS3Method base::summary
summary.frescalo = function(object, ...) {
  out = list(call = object$call,
             nsp = nrow(object$species),
        nsite = nrow(object$sites),
        nt = nrow(object$times),
        n_obs = object$n_obs,
        n_weights = object$n_weights,
        phi_target = object$phi_target,
        phi_obs_quant = quantile(object$sites$phi_obs, probs = c(.25,.5,.75, object$phi_prob)),
        mean_nsp = c(mean(object$sites$species_count), mean(object$sites$spnum_obs), mean(object$sites$spnum_est)))
  class(out) = "summary.frescalo"
  out
}

#' @rdname summary.frescalo
#' @param x An object returned by \code{summary.frescalo}.
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
  cat("\n Quantiles of observed phi:\n")
  print(x$phi_obs_quant, digits = 2)
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
#' @returns A data frame with species frequencies across sites and the following columns and corresponding names in the original fortran output:
#'
#' \tabular{llll}{
#' \strong{column} \tab \strong{fortran} \tab \strong{meaning}\cr
#' \code{site} \tab \code{Location} \tab Name of site.\cr
#' \code{species} \tab \code{Species} \tab Species name.\cr
#' \code{observed} \tab \code{Pres} \tab Whether species was recorded at site.\cr
#' \code{freq_obs} \tab \code{Freq} \tab Observed neigbourhood frequency of species. \cr
#' \code{freq_est} \tab \code{Freq_1} \tab Neighbourhood frequency of species after rescaling. \cr
#' \code{freq_est_se} \tab \code{SD_Frq1} \tab Standard error of \code{freq_est}. \cr
#' \code{rank} \tab \code{Rank} \tab Rank of estimated species frequency in the neighbourhood. \cr
#' \code{rank_scaled} \tab \code{Rank_1} \tab Rescaled rank, rank divided by estimated number of species.\cr
#' }
#'
#' @examples
#' data(bryophyte)
#' weights = compute_weights(hectad_locations, vascular_plants)
#' fr = frescalo(bryophyte, weights)
#' head(frequencies(fr))
#' @export
frequencies = function(object) {
  # freqs = object$freqs # Results in additional copy of large table.
  #freqs$spec_id = NULL
  #freqs$samp_id = NULL
  #setDF(freqs)
  frq = object$freqs
  frq$spec_id = NULL
  frq$site_id = NULL
  frq
}


#' Extract time factors from a frescalo object.
#'
#' @param object An object as returned from the frescalo function.
#'
#' @returns A data frame with time factors across species and the following columns and corresponding names in the original fortran output:
#'
#' \tabular{llll}{
#' \strong{column} \tab \strong{fortran} \tab \strong{meaning}\cr
#' \code{species} \tab \code{Species} \tab Species name.\cr
#' \code{time} \tab \code{Time} \tab Time period.\cr
#' \code{tf} \tab \code{TFactor} \tab Estimated time factor.\cr
#' \code{tf_se} \tab \code{St_Dev} \tab Standard error of time factor. \cr
#' \code{count} \tab \code{Count} \tab Number of observed occurrences of the species at given time. \cr
#' \code{occ_adj} \tab \code{spt} \tab Number of occurrences after down-weighting site time combinations with low proportion of benchmark species found. \cr
#' \code{occ_est} \tab \code{est} \tab Estimated number of occurrences. Should be the same as \code{occ_adj} if algorithm converged. \cr
#' \code{occ_0} \tab \code{N>0.00} \tab Number of locations with non-zero occurrence probability.\cr
#' \code{occ_098} \tab \code{N>0.98} \tab Number of locations with occurrence probability greater than 0.98.\cr
#' \code{conv_tfcalc} \tab  \tab Binary variable indicating that the routine calculating time factors converged.\cr
#' }
#'
#' @examples
#' weights = compute_weights(hectad_locations, vascular_plants)
#' fr = frescalo(bryophyte, weights)
#' head(timefactors(fr))
#' @export
timefactors = function(object) {
  #if (is.null(object$tfs)) {
    #tfalc....
  #}
  tfs = object$tfs
  if (any(!tfs$conv_tfcalc)) {
    warning("Computation of time factors did no converge for all sites.")
  }
#  if (!full) {
#    tfs$occ_adj = tfs$occ_est = tfs$occ_0 = tfs$occ_098 = tfs$conv_tfcalc = NULL
#  }
  tfs$iter_tf = tfs$iter_tf_se = NULL
  tfs
}

#' Extract information about sites from a frescalo object.
#'
#' @param object An object as returned from the frescalo function.
#'
#' @returns A data frame with a row for each sites, and information about its neighbourhood.
#'
#' \tabular{llll}{
#' \strong{column} \tab \strong{fortran} \tab \strong{meaning}\cr
#' \code{site} \tab \code{Location} \tab Name of site.\cr
#' \code{species_count} \tab \code{No_spp} \tab Number of species observed at site.\cr
#' \code{phi_obs} \tab \code{Phi_in} \tab Values of phi (frequency-weighted mean frequency) computed from observed frequencies.\cr
#' \code{alpha} \tab \code{Alpha} \tab Sampling effort multiplier required to reach \code{phi_target}.\cr
#' \code{spnum_obs} \tab \code{Spnum_in} \tab Sum of observed neighbourhood frequencies, a raw estimate of number of species. \cr
#' \code{spnum_est} \tab \code{Spnum_out} \tab Sum of neighbourhood frequencies after rescaling, an estimate of number of species. \cr
#' \code{wgt_n2} \tab \code{Wgt_n2} \tab Ratio of square of sum of neighbourhood weights to the sum of squared neighbourhood weights.
#'                                       Higher values indicate higher similarity among sites in neighbourhood. \cr
#' \code{phi_out} \tab \code{Phi_out} \tab Values of phi after rescaling. Should be close to \code{phi_target} if algorithm converged. \cr
#' \code{iter_fresca} \tab \code{Iter} \tab Number of iterations required for the rescaling routine. \cr
#' \code{conv_fresca} \tab  \tab Indicates whether the rescaling routine converged.\cr
#' }
#'
#'
#' @examples
#' weights = compute_weights(hectad_locations, vascular_plants)
#' fr = frescalo(bryophyte, weights)
#' head(sites(fr))
#' @export
sites = function(object) {
  #if (is.null(object$tfs)) {
  #tfalc....
  #}
  sites = object$sites
  sites$site_id = NULL
  sites
}



#' Extract the species that are considered benchmarks for each site.
#'
#' @param object An object as returned from the frescalo function.
#'
#' @returns A dataframe with the species used as benchmarks for each site.
#' @export
#'
#' @examples
#' data(bryophyte)
#' weights = compute_weights(hectad_locations, vascular_plants)
#' fr = frescalo(bryophyte, weights)
#' head(benchmark_species(fr))
benchmark_species = function(object) {
  species = site = bwght = rank_scaled = bench = NULL
  sp = as.data.table(object$species)
  frq = as.data.table(object$freqs)
  # Below need to match computation in benchmark_proportions()
  bs = sp[frq, list(site, species, bench = bwght * (rank_scaled < object$Rstar | rank == 1)), on = "spec_id"][bench > .99]
  bs$bench = NULL
  setDF(bs)
  bs
}


#' Extract the proportion of benchmark species found for each site and time.
#' These are used for effort correction when computing time factors.
#'
#' @param object An object as returned from the frescalo function.
#'
#' @returns A dataframe with the proportion of benchmark species found.
#' @export
#'
#' @examples
#' data(bryophyte)
#' weights = compute_weights(hectad_locations, vascular_plants)
#' fr = frescalo(bryophyte, weights)
#' head(benchmark_proportions(fr))
benchmark_proportions = function(object) {
  s_eff = as.data.table(object$sampling_effort)
  times = as.data.table(object$times[, c("time", "time_id")])
  s_eff = times[s_eff, on = "time_id"]
  sites = as.data.table(object$site[, c("site", "site_id")])
  s_eff = sites[s_eff, on = "site_id"]
  setDF(s_eff)
  s_eff$time_id = s_eff$site_id = NULL
  s_eff
}




check_phi = function(object, prob = .985, plot = FALSE, silent = FALSE) {
  obs_p = stats::quantile(object$sites$phi_obs, probs = prob)
  phi_ok = obs_p < object$phi_target
  if (!silent) {
    message(paste0(round(100 * prob, 1), " percentile of input phi = ", round(obs_p,2), "\n",
                "Target phi = ", round(object$phi_target, 2)))
  }
  if (!phi_ok) {
    warning("Target value of phi may be too small.")
  }
  if (plot) {
    graphics::hist(object$sites$phi_obs, xlim = c(0, 1), xlab = "phi", main = "")
    graphics::abline(v = object$phi_target, col = "red")
  }
  if (!all(object$sites$conv_fresca)) {
    conv_fail = !object$sites$conv_fresca
    warning(paste("phi did not converge to phi_target for all sites, increasing max.iter might help. Convergence failed for site(s):",
                  paste(object$sites$site[conv_fail], collapse = ", ")))
  }
  phi_ok
}


# Is this meaningful?
check_r = function(object) {
  graphics::hist(object$freqs$rank_scaled, xlab = "R", main = "")
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
#' weights = compute_weights(hectad_locations, vascular_plants)
#' fr = frescalo(bryophyte, weights)
#' check_rescaling(fr, max_sites = 50)
#' @export
check_rescaling = function(object, max_sites = 500) {
  scaled = NULL
  if (nrow(object$sites) > max_sites) {
    keep = object$freqs[["site"]] %in% sample(object$sites[["site"]], max_sites)
  } else {
    keep = TRUE
  }
  pldat = object$freqs[keep, c("site", "rank", "rank_scaled", "freq_obs", "freq_est")]
  pldat = data.table::as.data.table(pldat)
  pldat[, rank := as.numeric(rank)]
  pldat = data.table::melt(pldat,
                           id.vars = "site", measure.vars = list(freq = c("freq_obs", "freq_est"), rank = c("rank", "rank_scaled")),
                           variable.name = "scaled", variable.factor = FALSE)
  pldat[, scaled := c("unscaled", "rescaled")[as.integer(scaled)]]
  ggplot(data = pldat, aes(x = .data[["rank"]], y = .data[["freq"]], group = .data[["site"]])) +
    geom_line(alpha = 100/min(nrow(object$sites), max_sites)) + facet_wrap("scaled", scales = "free_x") + theme_minimal()

}





#' Compute estimated recording probabilities for a given level of effort.
#'
#' @param object An object as returned from the frescalo function.
#' @param species A character vector of species names for which to compute probabilities.
#' @param s Level of effort.
#'
#' @returns A data.frame with estimated occupancy probabilities.
#'
#' @note
#' Estimated probability of recording species under standard effort, sit = 1 meaning all benchmarks found.
#' This depends on the proportion of benchmarks.
#' @export
recording_probs = function(object, species, s = 1) {
  freq_est = tf = NULL
  freqs = data.table::as.data.table(object$freqs)
  tfs = data.table::as.data.table(object$tfs)
  keep_f = which(freqs$species %in% species)
  keep_tfs = which(tfs$species %in% species)
  out = merge(freqs[keep_f, c("species", "site", "freq_est")],
                          tfs[keep_tfs, c("species", "time", "tf")], allow.cartesian = TRUE)
  out[, "p_rec" := 1 - (1 - s * freq_est)^tf]
  out$tf = NULL
  out$freq_est = NULL
  setDF(out)
  out
}



