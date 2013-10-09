library("sfit")
set.seed(0xBEEF)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N <- 1000L

# Simulate genotypes
g <- sample(c("AA", "AB", "AB", "BB"), size=N, replace=TRUE)

# Simulate concentrations of allele A and allele B
X <- matrix(rexp(N), nrow=N, ncol=2L)
colnames(X) <- c("A", "B")
X[g == "AA", "B"] <- 0
X[g == "BB", "A"] <- 0
X[g == "AB",] <- X[g == "AB",] / 2

# Transform noisy X
xi <- matrix(rnorm(2*N, mean=0, sd=0.05), ncol=2L)
a0 <- c(0,0) + 0.3
A <- matrix(c(0.9, 0.1, 0.1, 0.8), nrow=2L, byrow=TRUE)
A <- apply(A, MARGIN=2L, FUN=function(u) u / sqrt(sum(u^2)))
Z <- t(a0 + A %*% t(X + xi))

# Add noise to Y
eps <- matrix(rnorm(2*N, mean=0, sd=0.05), ncol=2L)
Y <- Z + eps

layout(matrix(1:4, ncol=2L, byrow=TRUE))
par(mar=c(5,4,3,2) + 0.1)
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

qs <- c(0, 1, 2, 5, 10, 25)
Qs <- 100-qs
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
# Fit model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(Y, cex=0.8, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim, main="Y")

qs <- c(0, 1, 2, 5, 10, 25)
Qs <- rep(98, times=length(qs))
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
# Fit model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(Y, cex=0.8, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim, main="Y")

Qs <- 100-c(0, 1, 2, 5, 10, 25)
qs <- rep(2, times=length(Qs))
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

