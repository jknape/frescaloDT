# Compare against fortran output on Test.txt
#library(data.table)


# Verify euclidean distance calculation
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



# Verify habitat/ecological distance calculation
train_vasc = read.table("testdata/Training_vasc.txt", header = FALSE) # Locations of samples
setDT(train_vasc)

tv = dcast(train_vasc, V1~V2, fun.aggregate = length)

edist = eco_dist(tv[, -1], labels = tv$V1, subset = dists_fortran[,c("V1", "V2")], max_neigh = 100, method = "bray")

sim_fortran = read.table("testdata/sim.txt", header = FALSE)

expect_lt(max(abs(1 - sim_fortran$V4 - edist$dist)), 1e-3)

# Verify weights calculation

# This is how the fortran routine calculates weights. Weights less than 5e-5 are discarded.
#weights = dists_fortran[sim_fortran,  on = c("V1", "V2")][, weight := (1 - (i.V3 -1)^2 / 100^2)^4 * (1 - (V3 -1)^2 / 200^2)^4][weight > 5e-5]
#plot(weights$weight, bryophyte_weights$weight)


weights = compute_weights(sample_locs, tv[sample(1:404),], max_sp = 200, max_neigh = 100, cols = c(y = "V3", x = "V2", site = "V1"))
setDT(weights)

weights_fortran = read.table("testdata/weights.txt", header = FALSE)
setDT(weights_fortran)

joint_weights = weights[weights_fortran, on = c(site = "V1", neigh = "V2")]

# Match is not exact, because ties can lead to same site getting different ranks.
expect_lt(max(abs(joint_weights$V3 - joint_weights$weight)), .1)
expect_lt(abs(mean(joint_weights$V3 - joint_weights$weight)), .001)

