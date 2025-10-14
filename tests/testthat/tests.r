# Compare against fortran output on Test.txt
library(data.table)

freqs_fortran = read.table("testdata/freqs.txt", header = TRUE)
setDT(freqs_fortran)

weights = fread("testdata/weights.txt")
weights = weights[,1:3]
setnames(weights, c("samp", "samp1", "wgt"))

data = fread("testdata/Test.txt")

fr = frescalo(data, weights)

# Fortran program is truncating frequencies at 5e-5
freqs = frequencies(fr)
setDT(freqs)


comp_frq = freqs[freqs_fortran, on = c("samp" = "Location", "species" = "Species")]


expect_lt(max(abs(comp_frq$freq - comp_frq$Freq__)), 1e-4)
expect_lt(max(abs(comp_frq$Freq_1 - comp_frq$i.Freq_1)), 1e-4)
expect_lt(max(abs(comp_frq$SD_Frq1 - comp_frq$i.SD_Frq1)), 1e-4)
expect_lt(max(abs(comp_frq$rank - comp_frq$Rank)), 3)
expect_lt(max(abs(comp_frq$rank1 - comp_frq$Rank_1)), .02)


# Time factors

tfs_fortran = read.table("testdata/trends.txt", header = TRUE)
setDT(tfs_fortran)
tfs = timefactors(fr)
setDT(tfs)

comp_tfs = tfs[tfs_fortran, on = c("time" = "Time______", "species" = "Species__")]

# The fortran printing statement f7.3 cannot print tf > 10000, these are therefore ****** -> NA in the fortran version.
na_ind = which(is.na(comp_tfs$TFactor))
expect_gt(min(comp_tfs[na_ind]$tf) , 10000)
expect_gt(min(comp_tfs[na_ind]$se) , 10000)



comp_tfs = comp_tfs[-na_ind]


na_ind = which(is.na(comp_tfs$St_Dev))
expect_lt(comp_tfs$tf[na_ind] / comp_tfs$se[na_ind], .02)

comp_tfs = comp_tfs[-na_ind]

expect_lt(max(abs(comp_tfs$tf - comp_tfs$TFactor)), .01)
expect_lt(max(abs(comp_tfs$se - comp_tfs$St_Dev)), .02)
expect_equal(max(abs(comp_tfs$jtot - comp_tfs$X_Count)), 0)
expect_lt(max(abs(comp_tfs$sptot - comp_tfs$X___spt)), .1)
expect_equal(max(abs(comp_tfs$ic1 - comp_tfs$N.0.00)), 0)
expect_equal(max(abs(comp_tfs$ic2 - comp_tfs$N.0.98)), 0)
expect_lt(max(abs(comp_tfs$esttot - comp_tfs$X___est)), .1)

# TODO: check against output in stats file


# TODO: check what happens when there are duplicate records in input

fr = frescalo(data[c(1, 1:nrow(data))], weights)

# TODO: check varying names of input columns
