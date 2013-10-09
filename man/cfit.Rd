%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  cfit.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{cfit}
\alias{cfit}

\alias{cfit.matrix}

\title{Fits a K-dimensional simplex in M dimensions}

\description{
  Fits a K-dimensional simplex in M dimensions.
  A K-dimensional simplex is the K-dimensional generalization of a
  triangle.
}

\usage{
\method{cfit}{matrix}(y, k=ncol(y) + 1, dump=1, chopless=NULL, chopmore=NULL, maxiter=NULL, ...,
  retX=FALSE, cfit=getOption("cfit"), verbose=FALSE)
}

\arguments{
  \item{y}{Matrix or data frame of size IxN containing I rows of
    vectors in \eqn{R^N}.}
  \item{k}{The number of vertices of the fitted simplex. By default, the
    number of vertices is equal to the number of dimension (N) + 1.}
  \item{dump}{The output format.}
  \item{chopless, chopmore}{Lower and upper percentile thresholds at
    which extreme data points are assigned zero weights.}
  \item{maxiter}{"maximum number of REX steps". Default value is 60.}
  \item{...}{Named argument passed to the external 'cfit' program.}
  \item{retX}{If \code{\link[base:logical]{TRUE}}, an estimate of \code{X} is returned,
    otherwise not.}
  \item{cfit}{Shell command to call the 'cfit' executable.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, verbose output is displayed, otherwise not.}
}

\value{
  Returns a named \code{\link[base]{list}} structure elements:
   \item{\code{M}}{IxN \code{\link[base]{matrix}} where each rows is the coordinate for
     one of the vertices.}
   \item{\code{X}}{(optional) the IxN \code{\link[base]{matrix}} \eqn{X}.}
}

\details{
  Let \eqn{Y=(y_1, \ldots, y_I)} where \eqn{y_i=(y_{i1},\ldots,y_{iN})}
  is an observation in \eqn{N} dimensions.
  Let \eqn{M=(\mu_1,\ldots,\mu_K)} be the \eqn{K}-dimensional simplex
  where \eqn{mu_k} is a vertex in \eqn{N} dimensions.
  Let \eqn{X=(x_1,\ldots,X_I)} where \eqn{x_i=(x_{i1},\ldots,x_{iN})}.
  The simplex fitting algorithm decompose \eqn{Y} into:
  \deqn{
    Y \approx MX
  }
  such that \eqn{\sum_i x_{ik} = 1}.
}

\examples{# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
Z <- t(a0 + A \%*\% t(X + xi))

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
alpha <- c(0.10, 0.075, 0.05, 0.03, 0.01, 0.001)
fit <- cfit(Y, dump=2, alpha=alpha, q=2, Q=98)
Ms <- fit$M
col <- terrain.colors(length(Ms))
col[length(Ms)] <- "red"

plot(Y, cex=0.8, xlab=xlab, ylab=ylab, xlim=lim, ylim=lim, main="Y")

for (kk in seq(along=Ms)) {
  M <- Ms[[kk]]
  points(M, pch=19, cex=2.5, col=col[kk])
  lines(M, col=col[kk], lwd=2)
  text(M, cex=0.8, labels=kk)
}

legend("topright", pch=19, pt.cex=2, legend=c(alpha, "final"),
       col=col, title=expression(alpha), bg="#eeeeee")

apex <- which.min(apply(M, MARGIN=1, FUN=function(u) sum(u^2)))
a0hat <- M[apex,]
Ahat <- M[-apex,]
Ahat <- apply(Ahat, MARGIN=2, FUN=function(u) u / sqrt(sum(u^2)))
if (sum(Ahat[c(1,4)]^2) < sum(Ahat[c(2,3)]^2)) {
  Ahat <- matrix(Ahat[c(2,1,4,3)], nrow=2)
}
Ainv <- solve(Ahat)
Xhat <- t(Ainv \%*\% (t(Y) - a0hat))

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

plot(X[,1], Xhat[,1], cex=0.8, xlab=expression(X), ylab=expression(hat(X)), xlim=lim, ylim=lim)
points(X[,2], Xhat[,2], cex=0.8, col="red")
abline(a=0, b=1, lwd=2)
}

\author{
  Algorithm and C code/binary by Pratyaksha J. Wirapati,
  \email{wirapati@wehi.edu.au}.
  R wrapper by Henrik Bengtsson, \email{hb@maths.lth.se}.
}

\references{
 [1] P. Wirapati, & T. Speed, \emph{Fitting polyhedrial cones and
    simplices to multivariate data points}, Walter and Eliza Hall Institute
    of Medical Research, December 30, 2001.\cr
 [2] P. Wirapati and T. Speed, \emph{An algorithm to fit a simplex
    to a set of multidimensional points}, Walter and Eliza Hall Institute
    of Medical Research, January 15, 2002.\cr
}
\keyword{methods}
