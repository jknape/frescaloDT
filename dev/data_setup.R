#dists = read.table("tests/testthat/testdata/dist.txt", header = FALSE)
#colnames(dists) = c("hectad", "neighbour", "rank", "distance")


hectad_locations = read.table("tests/testthat/testdata/Samp_locations.txt", header = FALSE) # Locations of samples
colnames(hectad_locations) = c("hectad", "x", "y")


# Verify habitat/ecological distance calculation
train_vasc = read.table("tests/testthat/testdata/Training_vasc.txt", header = FALSE) # Locations of samples
colnames(train_vasc) = c("hectad", "species")
setDT(train_vasc)

vascular_plants = dcast(train_vasc, hectad~species, fun.aggregate = length)
setDF(vascular_plants)


bryophyte = read.table("tests/testthat/testdata/Test.txt", header = FALSE) # Locations of samples
colnames(bryophyte) = c("hectad", "species", "year")


usethis::use_data(bryophyte, hectad_locations, vascular_plants)
