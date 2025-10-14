tfcalc = function(data, freqs, species, sites, times, Rstar = .27, no_bench = NULL) {
  samp_id = spec_id = time_id = rank1 = bench = nbench = Freq_1 = bwght = spec = NULL # To avoid Notes in R CMD check
  no_bench = ""
  species[, bwght := 1]
  species[spec %in% no_bench, bwght := .001]

  sampeff = species[freqs, list(samp_id, spec_id, bench = bwght * (rank1 < Rstar | rank == 1)), on = "spec_id"][, nbench := sum(bench), by = "samp_id"]

  sampeff = sampeff[data, list(samp_id, spec_id, time_id, bench, nbench), on = c("samp_id", "spec_id")][, list(samp_eff = sum(bench)/nbench[1]), by = list(samp_id, time_id)]


  sampeff = dcast(sampeff, samp_id ~ time_id, value.var = "samp_eff", fill = 0)[sites[,list(samp_id)], on = "samp_id"]

  iocc = data.table(spec_id = rep(species$spec_id, each = nrow(times)), time_id = rep(times$time_id, nrow(species)))
  iocc0 = data[, list(occ = list(samp_id)), by = list(spec_id, time_id)]
  iocc = iocc0[iocc, on = list(spec_id, time_id)]
  setorder(iocc, spec_id, time_id)

  jind = iocc$spec_id
  tind = iocc$time_id
  iocc = iocc$occ
  fffl = freqs[order(spec_id), list(list(Freq_1)), by = spec_id]$V1

  ntf = length(jind)
  tfs = data.table(spec_id = jind, time_id = tind, tf = numeric(ntf), se = numeric(ntf),
                   jtot = integer(ntf), sptot = numeric(ntf), esttot = numeric(ntf), ic1 = integer(ntf), ic2 = integer(ntf), conv = logical(ntf))

  for (l in 1:length(jind)) {
    fff = fffl[[jind[l]]]
    smpint = sampeff[[tind[l]+1]]
    smpint[smpint == 0] = 1e-7 # ~ L37
    smpint[is.na(smpint)] = 1e-7 # ~ L37
    if (length(iocc[[l]]) < 1) {
      set(tfs, l, c("tf", "se", "jtot", "sptot", "esttot", "ic1", "ic2", "conv"),
          list(0,0, length(iocc[[l]]), 0, 0,0,0, TRUE))
    } else {
      set(tfs, l, c("tf", "se", "jtot", "sptot", "esttot", "ic1", "ic2", "conv"),
          tfcalc0(iocc[[l]], smpint, fff))
    }
  }
  tfs
}

tfcalc0 = function(iocc, smpint, fff, kmax = 100) {
  wgt = rep(1, length(smpint))
  iw = which(smpint < .0995)
  wgt[iw] = 0.005  + 10 * smpint[iw]
  pfac = smpint * fff
  ic1 = sum(pfac > 0)
  ipf = which(pfac > .98)
  pfac[ipf] = 0.98
  ic2 = length(ipf)
  plog = -log(1 - pfac)
  sptot = sum(wgt[iocc])
  jtot = length(iocc)
  tf = 1
  conv = FALSE
  k = 1
  mrel = 0
  iflag = 0
  while(k <= kmax & !conv) {
    estval = 1 - exp(-plog * tf)
    esttot = sum(wgt * estval)
    estvar = sum(wgt^2 * estval * (1-estval))
    if (abs(sptot-esttot) < 0.0005) {
      conv = TRUE
    } else {
      tf=tf*sptot/(esttot+0.0000001)
    }
    k = k +1
  }
  tf1=tf
  sptot1 = sptot + sqrt(estvar)
  k = 1
  conv2 = FALSE
  while(k <= kmax & !conv2) {
    estval = 1 - exp(-plog * tf1)
    esttt1 = sum(wgt*estval)
    if(abs(sptot1 - esttt1) < 0.0005) {
      conv2 = TRUE
    } else {
      tf1 = tf1 * sptot1 / (esttt1 + 0.0000001)
    }
    k = k + 1
  }
  list(tf = tf, se = tf1 - tf, jtot = jtot, sptot = sptot, esttot = esttot, ic1 = ic1, ic2 = ic2, conv = conv)
}




