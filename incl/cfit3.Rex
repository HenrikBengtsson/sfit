library(R.basic);  # plot3d(), lines3d()

filename <- system.file("data-ex/dye2.dat", package="sfit")
Y <- as.matrix(read.table(filename, sep="\t"))
#names(Y) <- c("x", "y", "z", "w")

Ms <- cfit(Y, dump=2)

layout(matrix(1:4, nrow=2, byrow=TRUE))
opar <- par(mar=c(1,2,2,1)+0.1, xpd=TRUE)

phi <- c(35, 35, 70, 70)
theta <- c(35, -35, 35, -35)

cols <- terrain.colors(length(Ms))
cols <- rainbow(length(Ms))
cols[length(Ms)] <- "blue";

for (view in 1:4) {
  plot3d(Y[,1:3], pch=".", xlab="dye 1", ylab="dye 2", zlab="dye 3", phi=phi[view], theta=theta[view], box=FALSE)
    
  for (iter in seq(Ms)) {
    M <- Ms[[iter]]
#    colnames(M) <- c("x", "y", "z", "w")
    
    points3d(M, pch="o", col=cols[iter])
    lines3d(M);
  }
}
  
par(opar)
