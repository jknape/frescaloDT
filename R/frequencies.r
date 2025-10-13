# Temporary function, should not be needed later.
fwf = function(data, weights, return_type = "dt") {
  setDT(data)
  setnames(data, c("samp", "spec", "time"))

  weights = weights[,1:3]
  setnames(weights, c("samp", "samp1", "wgt"))


  stopifnot(setequal(unique(weights$samp), unique(weights$samp1))) # How handle this??

  sites = data.table(samp = sort(unique(c(weights$samp))))[, samp_id := 1:.N]


  # Filter sites in data not present in weights
  exclude_sites = setdiff(unique(data$samp), sites$samp)
  if (length(exclude_sites) > 0) {
    message(paste("Site(s)", paste(exclude_sites, collapse = ", "), "not present in weights, removed."))
  }
  data =  data[!(samp %in% exclude_sites)]

  times = data.table(time = sort(unique(data$time)))[, time_id := 1:.N]

  species = data.table(spec = sort(unique(data$spec)))[, spec_id := 1:.N] # Note: species may have been removed, if only present in excluded sites.

  weights[ , samp_id := sites$samp_id[pmatch(weights$samp, sites$samp ,duplicates.ok = TRUE)]]
  weights[ , samp1_id := sites$samp_id[pmatch(weights$samp1, sites$samp ,duplicates.ok = TRUE)]]

  data[ , spec_id := species$spec_id[match(data$spec, species$spec)]]
  data[ , samp_id := sites$samp_id[match(data$samp, sites$samp)]]
  data[ , time_id := times$time_id[match(data$time, times$time)]]

  fwf0(data,weights, sites, species, return_type)
}

# Two functions only necessary
fwf0 = function(data, weights, sites, species, return_type = "dt") {

  # L254
  #swgt = weights[, .(wgttot = sum(wgt), wgt2 = sum(wgt^2)),by = samp]
  swgt = weights[, .(wgttot = sum(wgt), wgt2 = sum(wgt^2)), by = samp_id]
  setorder(swgt, samp_id)

  set(sites, j = "wn2", value = swgt$wgttot^2 / (swgt$wgt2 + 1e-12)) # Save for computation in fresca

  swgt_inv = weights[, .(neighs = list(samp_id), neighw = list(wgt)), by = samp1_id]
  setorder(swgt_inv, samp1_id)

  stopifnot(identical(sites$samp_id, swgt$samp_id))
  stopifnot(identical(sites$samp_id, swgt_inv$samp1_id))

  occ = data[, .(occ_ind = list(unique(spec_id))), by = samp_id]
  occ = occ[sites, on = c("samp_id")]


  stopifnot(identical(sites$samp_id, occ$samp_id))

  # L267
  #idat = data[ , .(idat = .N, iitot = uniqueN(spec)), by = samp]

  #idata = matrix(0L, nrow = uniqueN(data$samp), ncol = uniqueN(data$spec))

  # Number of detections
  #idata = dcast(data,  samp ~ spec, fun.aggregate = length)[, -1]


  ffij = matrix(0, nrow = nrow(sites), ncol = nrow(species))


  ineighs = swgt_inv$neighs
  iweigh = swgt_inv$neighw
  wgttot = swgt$wgttot
  occ_ind = occ$occ_ind
  for (ii in 1:nrow(ffij)) {
    occs = occ_ind[[ii]]
    nd = length(occs)
    if (nd > 0) {
      ineis <- ineighs[[ii]]
      ffij[ineis, occs] <- ffij[ineis, occs, drop = FALSE] + (iweigh[[ii]] / (wgttot[ineis] + 1e-10)) %o% (rep(1, nd))
    }
  }


  if (identical(return_type, "matrix")) {
      set(sites, j = "occ_ind", value = occ_ind)
      return(ffij)
  }
  freqs = data.table(samp_id = rep(1:nrow(ffij), ncol(ffij)), spec_id = rep(1:ncol(ffij), each = nrow(ffij)), freq = as.vector(ffij), pres = 0L)
  pind = do.call(c, mapply(function(j,i) {i + (j-1L)*nrow(ffij)} , occ_ind, 1:nrow(ffij)))
  set(freqs, i = pind, j = "pres", 1L)
  stopifnot(identical(freqs[samp_id == 114 & pres]$spec_id, occ_ind[[114]]))
  freqs
}





