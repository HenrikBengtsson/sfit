filename <- system.file("data-ex/dye2.dat", package="sfit")
Y <- as.matrix(read.table(filename, sep="\t"))

Ms <- sfit::cfit(Y, dump=2L)
print(Ms)

library("sfit")
Ms <- cfit(Y, dump=2L)
print(Ms)
