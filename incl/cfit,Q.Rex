# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
set.seed(0xbeef)
N <- 1000

# Simulate genotypes
g <- sample(c("AA", "AB", "AB", "BB"), size=N, replace=TRUE)

# Simulate concentrations of allele A and allele B
X <- matrix(rexp(N), nrow=N, ncol=2)
colnames(X) <- c("A", "B")
X[g == "AA", "B"] <- 0
X[g == "BB", "A"] <- 0
X[g == "AB",] <- X[g == "AB",] / 2

# Transform noisy X
xi <- matrix(rnorm(2*N, mean=0, sd=0.05), ncol=2)
a0 <- c(0,0)+0.3
A <- matrix(c(0.9, 0.1, 0.1, 0.8), nrow=2, byrow=TRUE)
A <- apply(A, MARGIN=2, FUN=function(u) u / sqrt(sum(u^2)))
Z <- t(a0 + A %*% t(X + xi))

# Add noise to Y
eps <- matrix(rnorm(2*N, mean=0, sd=0.05), ncol=2)
Y <- Z + eps

layout(matrix(1:4, ncol=2, byrow=TRUE))
par(mar=c(5,4,3,2)+0.1)
xlab <- "Allele A"
ylab <- "Allele B"
lim <- c(-0.5,8)
plot(X, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim)
points(Z, col="blue")
points(Y, col="red")

legend("topright", pch=19, pt.cex=2, legend=c("X", "Z", "Y"), 
       col=c("black", "blue", "red"), title="Variables:", bg="#eeeeee")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
plot(Y, cex=0.8, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim, main="Y")

Qs <- 100-c(0,1,2,5,10,25)
qs <- rep(2, length(Qs))
col <- terrain.colors(length(qs))
fits <- list()
for (kk in seq(along=qs)) {
  q <- qs[kk]
  Q <- Qs[kk]
  fit <- cfit(Y, q=q, Q=Q)
  fits[[kk]] <- fit

  M <- fit$M
  points(M, pch=19, cex=2.5, col=col[kk])
  lines(M, col=col[kk], lwd=2)
  text(M, cex=0.8, labels=kk)
}

legend("topright", pch=19, pt.cex=2, legend=sprintf("(%d,%d)", qs, Qs),
       col=col, title="(q,Q):", bg="#eeeeee")



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Backtransform
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
M <- fits[[1]]$M
apex <- which.min(apply(M, MARGIN=1, FUN=function(u) sum(u^2)))
a0hat <- M[apex,]
Ahat <- M[-apex,]
Ahat <- apply(Ahat, MARGIN=2, FUN=function(u) u / sqrt(sum(u^2)))
if (sum(Ahat[c(1,4)]^2) < sum(Ahat[c(2,3)]^2)) {
  Ahat <- matrix(Ahat[c(2,1,4,3)], nrow=2)
}
Ainv <- solve(Ahat)
Xhat <- t(Ainv %*% (t(Y) - a0hat))

cat("Estimated apex:\n")
print(a0hat);

cat("True A:\n")
print(A)

cat("Estimated A:\n")
print(Ahat)

plot(Xhat, cex=0.8, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim, main=expression(hat(X)))
x1 <- par("usr")[2]
y1 <- par("usr")[4]
lines(x=c(0,x1), y=c(0,0), col="red", lwd=2)
lines(x=c(0,0), y=c(0,y1), col="red", lwd=2)
lines(x=c(0,x1), y=c(0,y1), col="blue", lwd=2)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Backtransform
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
M <- fits[[length(fits)]]$M
apex <- which.min(apply(M, MARGIN=1, FUN=function(u) sum(u^2)))
a0hat <- M[apex,]
Ahat <- M[-apex,]
Ahat <- apply(Ahat, MARGIN=2, FUN=function(u) u / sqrt(sum(u^2)))
if (sum(Ahat[c(1,4)]^2) < sum(Ahat[c(2,3)]^2)) {
  Ahat <- matrix(Ahat[c(2,1,4,3)], nrow=2)
}
Ainv <- solve(Ahat)
Xhat <- t(Ainv %*% (t(Y) - a0hat))

cat("Estimated apex:\n")
print(a0hat);

cat("True A:\n")
print(A)

cat("Estimated A:\n")
print(Ahat)

plot(Xhat, cex=0.8, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim, main=expression(hat(X)))
x1 <- par("usr")[2]
y1 <- par("usr")[4]
lines(x=c(0,x1), y=c(0,0), col="red", lwd=2)
lines(x=c(0,0), y=c(0,y1), col="red", lwd=2)
lines(x=c(0,x1), y=c(0,y1), col="blue", lwd=2)
