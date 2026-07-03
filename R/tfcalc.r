benchmark_proportions0 = function(data, freqs, species, Rstar = .27, bench_exclude = NULL) {
  site_id = spec_id = time_id = rank_scaled = bench = nbench = freq_est = bwght = spec = NULL # To avoid Notes in R CMD check
  species[, bwght := 1]
  species[species %in% bench_exclude, bwght := .001]
  bench_prop = species[freqs, list(site_id, spec_id, bench = bwght * (rank_scaled < Rstar | rank == 1)), on = "spec_id"][, nbench := sum(bench), by = "site_id"]
  bench_prop = bench_prop[data, list(site_id, spec_id, time_id, bench, nbench), on = c("site_id", "spec_id")][, list(samp_eff = sum(bench)/nbench[1]), by = list(site_id, time_id)]
  setorderv(bench_prop, cols = c("site_id", "time_id"))
  bench_prop
}

# Use sum of observed frequencies instead of proportion of benchmark species found. Avoids Rstar. Might fail if whole community is declining.
benchmark_proportions0f = function(data, freqs, species, Rstar = .27, bench_exclude = NULL) {
  site_id = spec_id = time_id = rank_scaled = bench = nbench = freq_est = bwght = spec = NULL # To avoid Notes in R CMD check
  species[, bwght := 1]
  species[species %in% bench_exclude, bwght := .001]
  bench_prop = species[freqs, list(site_id, spec_id, freq_est = freq_est*bwght), on = "spec_id"][, fsum := sum(freq_est), by = "site_id"]
  bench_prop = bench_prop[data, list(site_id, spec_id, time_id, freq_est, fsum), on = c("site_id", "spec_id")][, list(samp_eff = sum(freq_est)/fsum[1]), by = list(site_id, time_id)]
  setorderv(bench_prop, cols = c("site_id", "time_id"))
  bench_prop
}


tfcalc = function(data, freqs, species, sites, times, sampeff) {
  site_id = spec_id = time_id = rank_scaled = freq_est = spec = NULL # To avoid Notes in R CMD check

  iocc = data.table(spec_id = rep(species$spec_id, each = nrow(times)), time_id = rep(times$time_id, nrow(species)))
  iocc0 = data[, list(occ = list(site_id)), by = list(spec_id, time_id)]
  iocc = iocc0[iocc, on = list(spec_id, time_id)]
  setorder(iocc, spec_id, time_id)

  jind = iocc$spec_id
  tind = iocc$time_id
  iocc = iocc$occ
  f_l = freqs[order(spec_id), list(list(freq_est)), by = spec_id]$V1

  ntf = length(jind)
  tfs = data.table(spec_id = jind, time_id = tind, tf = numeric(ntf), tf_se = numeric(ntf),
                   count = integer(ntf), occ_adj = numeric(ntf), occ_est = numeric(ntf), occ_0 = integer(ntf), occ_098 = integer(ntf), iter_tf = integer(ntf), iter_tf_se = integer(ntf),  conv_tfcalc = logical(ntf), dev = numeric(ntf))

  sampeff = dcast(sampeff, site_id ~ time_id, value.var = "samp_eff", fill = 0)[sites[,list(site_id)], on = "site_id"] # joining by sites fills in missing values for sites missing in sampeff

  for (l in 1:length(jind)) {
    f_j = f_l[[jind[l]]]
    s_t = sampeff[[tind[l]+1]]
    s_t[s_t == 0] = 1e-7 # ~ L37
    s_t[is.na(s_t)] = 1e-7 # ~ L37
    if (length(iocc[[l]]) < 1) {
      set(tfs, l, c("tf", "tf_se", "count", "occ_adj", "occ_est", "occ_0", "occ_098", "iter_tf", "iter_tf_se", "conv_tfcalc", "deviance"),
          list(0,0, length(iocc[[l]]), 0, 0,0,0, 0, 0, TRUE, 0))
    } else {
      set(tfs, l, c("tf", "tf_se", "count", "occ_adj", "occ_est", "occ_0", "occ_098",  "iter_tf", "iter_tf_se",  "conv_tfcalc", "deviance"),
          tfcalc0(iocc[[l]], s_t, f_j))
    }
  }
  tfs
}

tfcalc0 = function(iocc, s_t, f_j, kmax = 500) {
  wgt = rep(1, length(s_t))
  iw = which(s_t < .0995)
  wgt[iw] = 0.005  + 10 * s_t[iw]
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

