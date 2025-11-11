benchmark_proportions = function(data, freqs, species, Rstar = .27, bench_exclude = NULL) {
  location_id = spec_id = time_id = rank_scaled = bench = nbench = Freq_1 = bwght = spec = NULL # To avoid Notes in R CMD check

  species[, bwght := 1]
  species[spec %in% bench_exclude, bwght := .001]

  bench_prop = species[freqs, list(location_id, spec_id, bench = bwght * (rank_scaled < Rstar | rank == 1)), on = "spec_id"][, nbench := sum(bench), by = "location_id"]

  bench_prop = bench_prop[data, list(location_id, spec_id, time_id, bench, nbench), on = c("location_id", "spec_id")][, list(samp_eff = sum(bench)/nbench[1]), by = list(location_id, time_id)]
  setorderv(bench_prop, cols = c("location_id", "time_id"))
  bench_prop
}

tfcalc = function(data, freqs, species, sites, times, sampeff) {
  location_id = spec_id = time_id = rank_scaled = Freq_1 = spec = NULL # To avoid Notes in R CMD check

  iocc = data.table(spec_id = rep(species$spec_id, each = nrow(times)), time_id = rep(times$time_id, nrow(species)))
  iocc0 = data[, list(occ = list(location_id)), by = list(spec_id, time_id)]
  iocc = iocc0[iocc, on = list(spec_id, time_id)]
  setorder(iocc, spec_id, time_id)

  jind = iocc$spec_id
  tind = iocc$time_id
  iocc = iocc$occ
  f_l = freqs[order(spec_id), list(list(Freq_1)), by = spec_id]$V1

  ntf = length(jind)
  tfs = data.table(spec_id = jind, time_id = tind, tf = numeric(ntf), tf_se = numeric(ntf),
                   n_obs = integer(ntf), sptot = numeric(ntf), esttot = numeric(ntf), ic1 = integer(ntf), ic2 = integer(ntf), iter_tf = integer(ntf), iter_tf_se = integer(ntf),  conv = logical(ntf))

  sampeff = dcast(sampeff, location_id ~ time_id, value.var = "samp_eff", fill = 0)[sites[,list(location_id)], on = "location_id"]

  for (l in 1:length(jind)) {
    f_j = f_l[[jind[l]]]
    s_t = sampeff[[tind[l]+1]]
    s_t[s_t == 0] = 1e-7 # ~ L37
    s_t[is.na(s_t)] = 1e-7 # ~ L37
    if (length(iocc[[l]]) < 1) {
      set(tfs, l, c("tf", "tf_se", "n_obs", "sptot", "esttot", "ic1", "ic2", "iter_tf", "iter_tf_se", "conv"),
          list(0,0, length(iocc[[l]]), 0, 0,0,0, 0, 0, TRUE))
    } else {
      #browser()
      set(tfs, l, c("tf", "tf_se", "n_obs", "sptot", "esttot", "ic1", "ic2",  "iter_tf", "iter_tf_se",  "conv"),
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
  list(tf = tf, tf_se = tf1 - tf, n_obs = n_obs, sptot = sptot, esttot = esttot, ic1 = ic1, ic2 = ic2, iter_tf1 = k, iter_tf_se = k2, conv = conv)
}

# # Seems not much faster, maybe slight gain in memory allocation.
# tfcalc0_Rcpp = function(iocc, s_t, f_j, kmax = 100) {
#   wgt = rep(1, length(s_t))
#   iw = which(s_t < .0995)
#   wgt[iw] = 0.005  + 10 * s_t[iw]
#   sf = s_t * f_j
#   ic1 = sum(sf > 0)
#   ipf = which(sf > .98)
#   sf[ipf] = 0.98
#   ic2 = length(ipf)
#   Q = -log(1 - sf)
#   sptot = sum(wgt[iocc])
#   n_obs = length(iocc)
#   tf = 1
#
#   res1 = tf_iter(tf, wgt,Q,sptot)
#
#   sptot1 = sptot + sqrt(res1[3])
#
#   res2 = tf_iter(res1[1], wgt, Q, sptot1)
#
#   list(tf = res1[1], se = res2[1] - res1[1], n_obs = n_obs, sptot = sptot, esttot = res1[2], ic1 = ic1, ic2 = ic2, iter_tf1 = res1[4], iter_tf_se = res2[4], conv = res1[4] < kmax)
# }
