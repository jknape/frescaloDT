# Compute neighbourhoud frequencies
nfcalc = function(data, weights, sites, species, return_type = "dt") {
  weight = site_id = spec_id = neigh_id = pres = NULL # To avoid Notes in R CMD check
  # L254
  #swgt = weights[, .(wgttot = sum(wgt), wgt2 = sum(wgt^2)),by = samp]
  swgt = weights[, list(wgttot = sum(weight), wgt2 = sum(weight^2)), by = site_id]
  setorder(swgt, site_id)

  set(sites, j = "wgt_n2", value = swgt$wgttot^2 / (swgt$wgt2 + 1e-12)) # Save for computation in fresca

  swgt_inv = weights[, list(neighs = list(site_id), neighw = list(weight)), by = neigh_id]
  setorder(swgt_inv, neigh_id)

  stopifnot(identical(sites$site_id, swgt$site_id))
  stopifnot(identical(sites$site_id, swgt_inv$neigh_id))

  occ = data[, list(occ_ind = list(unique(spec_id))), by = site_id]
  occ = occ[sites, on = c("site_id")]

  set(sites, j = "species_count", value = sapply(occ$occ_ind, length))

  stopifnot(identical(sites$site_id, occ$site_id))

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
  freqs = data.table(site_id = rep(1:nrow(ffij), ncol(ffij)), spec_id = rep(1:ncol(ffij), each = nrow(ffij)), freq_obs = as.vector(ffij), observed = FALSE)
  pind = do.call(c, mapply(function(j,i) {i + (j-1L)*nrow(ffij)} , occ_ind, 1:nrow(ffij), SIMPLIFY = FALSE))
  set(freqs, i = pind, j = "observed", TRUE)

  freqs
}

