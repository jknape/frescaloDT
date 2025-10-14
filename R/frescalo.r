#' Analyse species occurrence data with the frescalo algorithm of Hill 2012.
#'
#' @param data Data frame with the samples. First column should contain the name or id of the sampled site.
#'             Second column the observed species, and third column the time of observation.
#' @param weights Data frame with neighbourhood weights where the first column is the target site, second column is the
#'                neighbour, and third column is the weight in the neigbourhood of the target site.
#' @param phi_target Target value for adjusted frequency weighted mean frequencies.
#' @param Rstar Threshold for species to be considered as benchmarks when computing time factors.
#' @param bench_exclude Vector of names of species not to be used as benchmarks when computing time factors.
#'
#' @returns An object of class frescalo
#' @export
#'
#' @examples
frescalo = function(data, weights, phi_target = .74, Rstar = 0.27, bench_exclude = NULL) {
  samp_id = samp = spec_id = time_id = samp1_id = NULL # To avoid Notes in R CMD check

  setDT(data)
  setnames(data, c("samp", "spec", "time"))
  weights = weights[,1:3]
  setnames(weights, c("samp", "samp1", "wgt"))

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

  freqs[, c("Freq_1", "SD_Frq1", "rank", "rank1") := frescaDT2(.SD, sites, phi_target = phi_target), keyby = list(samp_id), .SDcols = c("samp_id", "freq")]

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
  out = list(freqs = freqs, tfs = tfs, sites = sites, species = species, times = times, excluded_sites = exclude_sites, phi_target = phi_target, Rstar = Rstar)
  class(out) = "frescalo"
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

