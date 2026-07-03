frescaloL = function(data, weights, phi_target = .74, Rstar = 0.27, phi_prob = .985, bench_exclude = NULL,
                    cols = NULL) {
  site_id = site = spec_id = time_id = neigh_id = phi0 = freq_obs = NULL # To avoid Notes in R CMD check

  if (is.null(cols)) {
    if (ncol(data) < 4 | ncol(weights) < 3) {
      stop("Arguments data and weights need to have three columns.")
    }
    data = data[,1:4]
    cols = c(colnames(data), colnames(weights)[2:3])
    names(cols) = c("site", "list", "species", "time", "neigh", "weight")
    data = data.table::as.data.table(data) # Should copy, otherwise may be reordered on return!
    setnames(data, c("site", "list", "species", "time"))
    weights = weights[,1:3]
    weights = data.table::as.data.table(weights)
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

  # Collapse lists for frequency computation
  fdata = data[, list(spec_id = unique(spec_id)), by = c("site_id", "time_id")]
  # Compute frequency weighted mean frequencies
  freqs = nfcalc(fdata, weights, sites, species)
  if (is.na(phi_target)) {
    phi_target = freqs[ , list(phi0 = sum(freq_obs^2)/sum(freq_obs)), by = "site_id"][, quantile(phi0, prob = phi_prob, names = FALSE)]
  }
  nc = nrow(freqs)
  set(freqs, j = c("freq_est", "freq_est_se", "rank", "rank_scaled"),
      value = list(numeric(nc), numeric(nc), integer(nc), numeric(nc))) # rank_scaled = R´, Hill P200.
  setkey(freqs, site_id) # Not needed, minimal improvement if any?
  freqs[, c("freq_est", "freq_est_se", "rank", "rank_scaled") := frescaDT(.SD, sites, phi_target = phi_target, irepmax = 100), keyby = list(site_id), .SDcols = c("site_id", "freq_obs")]


  # Compute effort per site and time as proportion of benchmarks found.
  sampling_effort = samp_effortL(data, freqs, species, Rstar = Rstar, bench_exclude = bench_exclude)

  # Compute time factors.
  tfs = tfcalcL(data, freqs, species, sites, times, sampling_effort)

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


samp_effortL = function(data, freqs, species, Rstar = .27, bench_exclude = NULL, method = "freq") {
  site_id = spec_id = time_id = rank_scaled = bench = nbench = freq_est = bwght = spec = NULL # To avoid Notes in R CMD check
  species[, bwght := 1]
  species[species %in% bench_exclude, bwght := .001]
  if (method == "bench") {
    bench_prop = species[freqs, list(site_id, spec_id, bench = bwght * (rank_scaled < Rstar | rank == 1)), on = "spec_id"][, nbench := sum(bench), by = "site_id"]
    bench_prop = bench_prop[data, list(list, site_id, spec_id, time_id, bench, nbench), on = c("site_id", "spec_id")][, list(samp_eff = sum(bench)/nbench[1]), by = list(site_id, list, time_id)]
  }
  if (method == "freq") {
    bench_prop = species[freqs, list(site_id, spec_id, freq_est = freq_est*bwght), on = "spec_id"][, fsum := sum(freq_est), by = "site_id"]
    bench_prop = bench_prop[data, list(list, site_id, spec_id, time_id, freq_est, fsum), on = c("site_id", "spec_id")][, list(samp_eff = sum(freq_est)/fsum[1]), by = list(site_id, list, time_id)]
  }
  browser()
  setorderv(bench_prop, cols = c("site_id", "time_id"))
  bench_prop
}




tfcalcL = function(data, freqs, species, sites, times, sampeff) {
  site_id = spec_id = time_id = rank_scaled = freq_est = spec = NULL # To avoid Notes in R CMD check

  iocc = data.table(spec_id = rep(species$spec_id, each = nrow(times)), time_id = rep(times$time_id, nrow(species)))
  iocc0 = data[, list(occ = list(list)), by = list(spec_id, time_id)]
  iocc = iocc0[iocc, on = list(spec_id, time_id)]
  setorder(iocc, spec_id, time_id)

  jind = iocc$spec_id
  tind = iocc$time_id
  iocc = iocc$occ
  f_l = freqs[order(spec_id), list(list(freq_est)), by = spec_id]$V1

  ntf = length(jind)
  tfs = data.table(spec_id = jind, time_id = tind, tf = numeric(ntf), tf_se = numeric(ntf),
                   count = integer(ntf), occ_adj = numeric(ntf), occ_est = numeric(ntf), occ_0 = integer(ntf), occ_098 = integer(ntf), iter_tf = integer(ntf), iter_tf_se = integer(ntf),  conv_tfcalc = logical(ntf), deviance = numeric(ntf))

  #sampeffOld = dcast(sampeff, site_id ~ time_id, value.var = "samp_eff", fill = 0)[sites[,list(site_id)], on = "site_id"]

  sampeff = dcast(sampeff, list ~ time_id, value.var = "samp_eff", fill = 0) #[sites[,list(site_id)], on = "site_id"]
  lists = unique(data[, c("site_id", "time_id","list")])[order(list)]
  lists[, weight := 1/.N, by = c("site_id", "time_id")]
  for (l in 1:length(jind)) {
    f_j = f_l[[jind[l]]][lists[time_id == tind[l]]$site_id]
    s_t = sampeff[[tind[l]+1]][lists[time_id == tind[l]]$list]
    # Temporary reindexing of iocc, should be implemented better
    iocc_l = rep(0, nrow(lists))
    iocc_l[iocc[[l]]] = 1
    iocc_l = iocc_l[lists[time_id == tind[l]]$list]
    stopifnot(sum(iocc_l) == length(iocc[[l]]))
    s_t[s_t == 0] = 1e-7 # ~ L37
    s_t[is.na(s_t)] = 1e-7 # ~ L37
    w_l = lists[time_id == tind[l]]$weight
    if (length(iocc[[l]]) < 1) {
      set(tfs, l, c("tf", "tf_se", "count", "occ_adj", "occ_est", "occ_0", "occ_098", "iter_tf", "iter_tf_se", "conv_tfcalc", "deviance"),
          list(0,0, length(iocc[[l]]), 0, 0,0,0, 0, 0, TRUE, 0))
    } else {
      set(tfs, l, c("tf", "tf_se", "count", "occ_adj", "occ_est", "occ_0", "occ_098",  "iter_tf", "iter_tf_se",  "conv_tfcalc", "deviance"),
          tfcalc0L(which(iocc_l == 1), s_t, f_j, w_l))
    }
    print(l)
  }
  tfs
}

tfcalc0L = function(iocc, s_t, f_j, w_l, kmax = 500) {
  wgt = rep(1, length(s_t))
  iw = which(s_t < .0995)
  wgt[iw] = 0.005  + 10 * s_t[iw]
  wgt = wgt * w_l
  sf = s_t * f_j
  ic1 = sum(sf > 0)
  ipf = which(sf > .98)
  sf[ipf] = 0.98
  ic2 = length(ipf)
  Q = -log(1 - sf)
  sptot = sum(wgt[iocc])
  n_obs = length(iocc)
  tf = 1
  conv = FALSE
  k = 1
  while(k <= kmax & !conv) {
    P = 1 - exp(-Q * tf)
    esttot = sum(wgt * P)
    estvar = sum(wgt^2 * P * (1-P))
    if (abs(sptot-esttot) < 0.0005) {
      conv = TRUE
    } else {
      tf=tf*sptot/(esttot+0.0000001)
    }
    k = k +1
  }
  dev = ifelse(length(iocc)>0, -2 * (sum(log(P[iocc])) + sum(log(1-P[-iocc]))), -2 * sum(log(1-P)))
  tf1=tf
  sptot1 = sptot + sqrt(estvar)
  k2 = 1
  conv2 = FALSE
  while(k2 <= kmax & !conv2) {
    P = 1 - exp(-Q * tf1)
    esttt1 = sum(wgt*P)
    if(abs(sptot1 - esttt1) < 0.0005) {
      conv2 = TRUE
    } else {
      tf1 = tf1 * sptot1 / (esttt1 + 0.0000001)
    }
    k2 = k2 + 1
  }
  list(tf = tf, tf_se = tf1 - tf, count = n_obs, occ_adj = sptot, occ_est = esttot, occ_0 = ic1, occ_098 = ic2, iter_tf1 = k, iter_tf_se = k2, conv_tfcalc = conv, deviance = dev)
}


