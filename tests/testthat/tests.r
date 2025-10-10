freqs_fortran = read.table("freqs.txt", header = TRUE)

# Fortran program is truncating frequencies at 5e-5
nrow(freqs)
nrow(freqs_fortran)
max(abs(freqs[freq > 5e-5][,freq] - freqs_fortran[,"Freq__"]))

identical(freqs[freq > 5e-5][,freq],freqs_fortran[,"Species"]) # Some rankings not identical due to numerical difference


plot(sort(ffij[1,], decreasing = TRUE), type = "l", col = adjustcolor("black", alpha = .1))
for (i in 2:nrow(ffij))
  points(sort(ffij[i,], decreasing = TRUE), type = "l", col = adjustcolor("black", alpha = .1))




# TFcalc
#setorder(fre, Samp, -Freq)
# Verify against fortran
# fre = fre[fre$Freq > 5e-5, ]
# max(abs(fre$Freq - freqs_fortran$Freq__))
# max(abs(fre$Freq_1 - freqs_fortran$Freq_1))
# max(abs(fre$SD_Frq1 - freqs_fortran$SD_Frq1))
# max(abs(fre$rank - freqs_fortran$Rank))
# max(abs(fre$rank1 - freqs_fortran$Rank_1))




tfs_fortran = read.table("trends.txt", header = TRUE)

plot(tfs_fortran$TFactor[is.finite(tfs_fortran$TFactor)], tfs$tf[is.finite(tfs_fortran$TFactor)])
tfdiffs = which(abs(tfs_fortran$TFactor - tfs$tf) >.005)
cbind(tfs[tfdiffs], tfs_fortran[tfdiffs, ])

max(abs(tfs_fortran$TFactor - tfs$tf), na.rm = TRUE)

# Note: The fortran printing statement f7.3 cannot print tf > 10000, these are therefore ****** in the trends file.
tfdiffs = which(tfs$tf > max(tfs_fortran$TFactor, na.rm = TRUE))
cbind(tfs[tfdiffs], tfs_fortran[tfdiffs, ])
