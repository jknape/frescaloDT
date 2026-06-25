# Compare against fortran output on Test.txt
library(data.table)

freqs_fortran = read.table("testdata/freqs.txt", header = TRUE)
setDT(freqs_fortran)


data("bryophyte")


bryophyte_weights = read.table("testdata/weights.txt", header = FALSE)

bryophyte_weights2 = compute_weights(hectad_locations, vascular_plants)

fr = frescalo(bryophyte, bryophyte_weights)

fr2 = frescalo(bryophyte, bryophyte_weights2)

# Fortran program is truncating frequencies at 5e-5
freqs = frequencies(fr)
freqs2 = frequencies(fr2)
setDT(freqs)
setDT(freqs2)




comp_frq = freqs[freqs_fortran, on = c("location" = "Location", "species" = "Species")]
comp_frq2 = freqs2[freqs_fortran, on = c("location" = "Location", "species" = "Species")]

expect_lt(max(abs(comp_frq$freq - comp_frq$Freq__)), 1e-4)
expect_lt(max(abs(comp_frq$Freq_1 - comp_frq$i.Freq_1)), 1e-4)
expect_lt(max(abs(comp_frq$SD_Frq1 - comp_frq$i.SD_Frq1)), 1e-4)
expect_lt(max(abs(comp_frq$rank - comp_frq$Rank)), 3)
expect_lt(max(abs(comp_frq$rank_scaled - comp_frq$Rank_1)), .02)


expect_lt(max(abs(comp_frq2$freq - comp_frq2$Freq__)), 1e-2)
expect_lt(max(abs(comp_frq2$Freq_1 - comp_frq2$i.Freq_1)), 5e-2)
expect_lt(max(abs(comp_frq2$SD_Frq1 - comp_frq2$i.SD_Frq1)), 5e-2)
expect_lt(max(abs(comp_frq2$rank_scaled - comp_frq2$Rank_1)), 5e-1)

# Time factors

tfs_fortran = read.table("testdata/trends.txt", header = TRUE)
setDT(tfs_fortran)
tfs = timefactors(fr)
tfs2 = timefactors(fr2)
setDT(tfs)
setDT(tfs2)

comp_tfs = tfs[tfs_fortran, on = c("time" = "Time______", "species" = "Species__")]
comp_tfs2 = tfs2[tfs_fortran, on = c("time" = "Time______", "species" = "Species__")]

# The fortran printing statement f7.3 cannot print tf > 10000, these are therefore ****** -> NA in the fortran version.
na_ind = which(is.na(comp_tfs$TFactor))
na_ind2 = which(is.na(comp_tfs2$TFactor))

expect_gt(min(comp_tfs[na_ind]$tf) , 10000)
expect_gt(min(comp_tfs[na_ind]$tf_se) , 10000)


comp_tfs = comp_tfs[-na_ind]
comp_tfs2 = comp_tfs2[-na_ind2]


na_ind = which(is.na(comp_tfs$St_Dev))
expect_lt(comp_tfs$tf[na_ind] / comp_tfs$tf_se[na_ind], .02)

na_ind2 = which(is.na(comp_tfs2$St_Dev))

comp_tfs = comp_tfs[-na_ind]
comp_tfs2 = comp_tfs2[-na_ind2]

expect_lt(max(abs(comp_tfs$tf - comp_tfs$TFactor)), .01)
expect_lt(max(abs(comp_tfs$tf_se - comp_tfs$St_Dev)), .02)
expect_equal(max(abs(comp_tfs$n_obs - comp_tfs$X_Count)), 0)
expect_lt(max(abs(comp_tfs$sptot - comp_tfs$X___spt)), .1)
expect_equal(max(abs(comp_tfs$ic1 - comp_tfs$N.0.00)), 0)
expect_equal(max(abs(comp_tfs$ic2 - comp_tfs$N.0.98)), 0)
expect_lt(max(abs(comp_tfs$esttot - comp_tfs$X___est)), .1)

expect_lt(max(abs(comp_tfs2$tf - comp_tfs2$TFactor)), .5)
expect_lt(max(abs(comp_tfs2$tf_se - comp_tfs2$St_Dev)), .8)
expect_equal(max(abs(comp_tfs2$n_obs - comp_tfs2$X_Count)), 0)
expect_lt(max(abs(comp_tfs2$sptot - comp_tfs2$X___spt)), .3)
expect_lt(max(abs(comp_tfs2$esttot - comp_tfs2$X___est)), .4)


# TODO: check against output in stats file


# TODO: check what happens when there are duplicate records in input

#fr3 = frescalo(bryophyte[c(1, 1:nrow(bryophyte)), ], bryophyte_weights)

# TODO: check that order of input doesn't matter

#fr4 = frescalo(bryophyte[sample(1:nrow(bryophyte)), ], bryophyte_weights[sample(1:nrow(bryophyte_weights)),])

# TODO: check varying names of input columns


# Check that helper functions work

expect_no_error(summary(fr))
expect_no_error(recording_probs(fr, species = c("Bry_1079")))
expect_no_error(recording_probs(fr, species = c("Bry_11", "Bry_126")))
expect_no_error(benchmark_species(fr))


