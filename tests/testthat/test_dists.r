# Compare against fortran output on Test.txt
#library(data.table)

dists_fortran = read.table("testdata/dist.txt", header = FALSE)

sample_locs = read.table("testdata/Samp_locations.txt", header = FALSE) # Locations of samples

dists = eco_dist(sample_locs[,c("V2", "V3")], labels = sample_locs$V1, method = "euclidean")

setDT(dists)
setorder(dists, "site", "dist_rank", "neigh")
setDT(dists_fortran)
setorder(dists_fortran, "V1", "V3","V2")

expect_equal(dists$site, dists_fortran$V1)
#expect_equal(dists$neigh, dists_fortran$V2) # Does not hold, ties are broken differently
expect_lt(max(abs(dists$dist - dists_fortran$V4)), 1)




train_vasc = read.table("testdata/Training_vasc.txt", header = FALSE) # Locations of samples
setDT(train_vasc)

tv = dcast(train_vasc, V1~V2)
