library("sfit")

filename <- system.file("data-ex/dye2.dat", package="sfit")
Y <- as.matrix(read.table(filename, sep="\t"))
Ms <- cfit(Y, dump=2L)
print(Ms)
