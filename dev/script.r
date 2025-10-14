
library(data.table)
setDTthreads(threads = 1)
#profvis::profvis({for (kk in 1:100){

#weights = read.table("../input/weights.txt")
weights = fread("tests/testthat/testdata/weights.txt")
weights = weights[,1:3]
setnames(weights, c("samp", "samp1", "wgt"))

#data = read.table("../input/occ.txt")
data = fread("tests/testthat/testdata/Test.txt")

fr = frescalo(data, weights)

check_phi(fr)

frescalo(data, weights, phi_target = 1.01)

#######################################
# Format data
#######################################

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


# Variant 1. ffij as matrix
system.time({

  ffij = fwf0(data, weights, sites, species, "matrix")

  nc = prod(dim(ffij))
  fre = data.frame(samp = vector(mode = mode(sites$samp), length = nc), spec = vector(mode = mode(species$spec), length = nc), Pres = integer(nc), Freq = numeric(nc), Freq_1 = numeric(nc),
                   SD_Frq1 = numeric(nc), rank = integer(nc), rank1 = numeric(nc)) # rank1 = R´, Hill P200.
  setDT(fre)

  #ns = nrow(sites)
  #fre.0 = data.frame(site = sites$samp, alpha = numeric(ns), spnum1 = numeric(ns), phi1 = numeric(ns), iter = integer(ns), conv = logical(ns))

  wn2 = sites$wn2
  occ_ind = sites$occ_ind
  for (i in 1:nrow(ffij)) {
    f = ffij[i,]
    jocc = integer(ncol(ffij))
    jocc[occ_ind[[i]]] <- 1L
    samp = sites$samp_id[i]
    itot = length(occ_ind[[i]])
    #fre[1:ncol(ffij) + (i-1)*ncol(ffij),] = fresca(f, wn2[i])
    set(fre, 1:ncol(ffij) + (i-1L) * ncol(ffij), "samp", sites$samp[i])
    fresca.DT(fre, i,f , wn2[i])
    #fre.0[i, 2:ncol(fre.0)] = fresca.0(f, wn2[i])
  }
})


# attempt data.table variant. Similar runtimes...
system.time({
  fre1 = fwf0(data, weights, sites, species, "dt")
  nc = nrow(fre1)
  set(fre1, j = c("Freq_1", "SD_Frq1", "rank", "rank1") ,
      value = list(numeric(nc), numeric(nc), integer(nc), numeric(nc))) # rank1 = R´, Hill P200.
  setkey(fre1, samp_id) # Not needed, minimal improvement if any?

  fre1[, c("Freq_1", "SD_Frq1", "rank", "rank1") := frescaDT2(.SD, sites), keyby = .(samp_id), .SDcols = c("samp_id", "freq")] # Note: info is added to 'sites' as side effect.
})


setorder(fre, samp_id, -freq, spec_id)


# Time factors
system.time({

  Rstar = .27

  no_bench = ""
  species[, bwght := 1]
  species[spec %in% no_bench, bwght := .001]

  #fre[ , bench := rank1 < Rstar][, nbench := sum(bench), by = "samp"]
  # sampeff = fre[data, ...
  sampeff = species[fre, .(samp, spec, bench = bwght * (rank1 < Rstar | rank == 1)), on = "spec"][, nbench := sum(bench), by = "samp"] # Slower and more memory (?), but necessary with excluded species
  sampeff = sampeff[data, .(samp, spec, time, bench, nbench), on = c("samp", "spec")][, .(samp_eff = sum(bench)/nbench[1]), by = .(samp, time)]

  setorder(sampeff, samp, time)

  setorder(data, spec, time, samp)


  sampeff = dcast(sampeff, samp ~ time, value.var = "samp_eff", fill = 0)[sites[,.(samp)], on = "samp"]

  #  lendat = dcast(data[, .N, by = .(spec, time)], spec ~ time, fill = 0) # Number of sites where species found at given time

  iocc = data.table(spec_id = rep(species$spec_id, each = nrow(times)), time_id = rep(times$time_id, nrow(species)))
  iocc0 = data[, .(occ = list(samp_id)), by = .(spec_id, time_id)]
  iocc = iocc0[iocc, on = .(spec_id, time_id)]
  setorder(iocc, spec_id, time_id)

  jind = iocc$spec_id
  tind = iocc$time_id
  iocc = iocc$occ
  fffl = fre[order(spec), .(list(Freq_1)), by = spec]$V1

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
          list(0,0, length(iocc[[l]]), 0, 0,0,0, TRUE)) # ic1 is not correct
    } else {
      set(tfs, l, c("tf", "se", "jtot", "sptot", "esttot", "ic1", "ic2", "conv"),
          tfcalc0(iocc[[l]], smpint, fff)) # ic1 is not correct
    }
  }
})


# Time factors
system.time({
  tfs = tfcalc(data, fre1, species, sites, Rstar = .27, no_bench = NULL)
})

