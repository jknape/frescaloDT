# Note: info is added to 'sites' as side effect.
frescaDT2 = function(dat, sites, phi_target = .74, irepmax = 100, fmax=0.99999, fmin=1.0E-10,tol=0.0003) {
  f = dat$freq
  sa_id = dat$samp_id[1]
  wn2 = sites$wn2[sa_id]
  alpha = 1
  conv = FALSE
  i = 1L
  ifz = which(f == 0)
  f[f> fmax] = fmax
  f[f < fmin] = fmin
  while (i <= irepmax & ! conv) {
    ff = -log(1-f)
    ffn = 1 - exp(-ff  * alpha)
    tot = sum(ffn)
    tot2 = sum(ffn ^ 2)
    phi = tot2 / tot
    an2 = tot^2 / tot2
    spnum = tot
    if (i < 20) {
      alpha = alpha * exp(1.86 * (log(1-phi) - log(1-phi_target)))
    } else {
      alpha = alpha * phi_target / phi
    }
    if(i == 1) {
      phi0=phi
      an21=an2
      spnum0=tot
    }
    if(abs(phi-phi_target) < tol) {
      conv = TRUE
    }
    i = i + 1
  }
  if (!conv) {
    warning("fresca failed to converge")
  }
  ffrs = 1-exp(alpha*log(1-f))
  set(sites, i = sa_id, j = c("alpha", "spnum", "phi0", "spnum0", "iter.fresca", "conv.fresca"), value = list(alpha, "spnum", phi0, spnum0, i, conv)) # Could preallocate, but not necessary?
  # All of below could be done later, depends only on f, spnum and alpha.
  jrank = order(-f + 1:length(f) * 1e-12)
  ffrs = 1-exp(alpha*log(1-f))
  sdffrs = sqrt(f * (1-f) / wn2)
  ffff = f + sdffrs
  ffff[1-ffff < 1e-12] = 1 - 1e-12
  fffff = f - sdffrs
  ffsd = 1 - exp(alpha * log(1 - ffff))
  fffsd = 1 - exp(alpha * log(1 - fffff))
  sdij = 0.5 * (ffsd - fffsd)
  ffrs[ifz] = 0 # Line 333
  list(ffrs, sdij, order(jrank), order(jrank)/spnum)
}


fresca.DT = function(dat, row_ind, f, wn2, phi_target = .74, irepmax = 100, fmax=0.99999, fmin=1.0E-10,tol=0.0003) {
  ns = length(f)
  alpha = 1
  conv = FALSE
  i = 1L
  ifz = which(f == 0)
  f[f> fmax] = fmax
  f[f < fmin] = fmin
  while (i <= irepmax & ! conv) {
    ff = -log(1-f)
    ffn = 1 - exp(-ff  * alpha)
    tot = sum(ffn)
    tot2 = sum(ffn ^ 2)
    phi = tot2 / tot
    an2 = tot^2 / tot2
    spnum = tot
    if (i < 20) {
      alpha = alpha * exp(1.86 * (log(1-phi) - log(1-phi_target)))
    } else {
      alpha = alpha * phi_target / phi
    }
    if(i == 1) {
      phi0=phi
      an21=an2
      spnum0=tot
    }
    if(abs(phi-phi_target) < tol) {
      conv = TRUE
    }
    i = i + 1
  }
  if (!conv) {
    warning("fresca failed to converge")
  }
  # All of below could be done later, depends only on f, spnum and alpha.
  jrank = order(-f + 1:ns * 1e-12)
  ffrs = 1-exp(alpha*log(1-f))
  sdffrs = sqrt(f * (1-f) / wn2)
  ffff = f + sdffrs
  ffff[1-ffff < 1e-12] = 1 - 1e-12
  fffff = f - sdffrs
  ffsd = 1 - exp(alpha * log(1 - ffff))
  fffsd = 1 - exp(alpha * log(1 - fffff))
  sdij = 0.5 * (ffsd - fffsd)
  ffrs[ifz] = 0 # Line 333
  set(dat, 1:ns + (row_ind-1L)*ns, c("spec",     "Pres", "Freq", "Freq_1", "SD_Frq1",       "rank",            "rank1"),
                                list(species$spec, jocc,      f,     ffrs,      sdij, order(jrank), order(jrank)/spnum))
} # Note species$spec does not seem to work if species not defined when frescalo.DT defined.

