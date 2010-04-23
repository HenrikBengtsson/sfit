#########################################################################/**
# @set "class=matrix"
# @RdocMethod cfit
#
# @title "Fits a K-dimensional simplex in M dimensions"
#
# \description{
#   @get "title".  
#   A K-dimensional simplex is the K-dimensional generalization of a
#   triangle.
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{Matrix or data frame of size IxN containing I rows of 
#     vectors in \eqn{R^N}.}
#   \item{k}{The number of vertices of the fitted simplex. By default, the
#     number of vertices is equal to the number of dimension (N) + 1.}
#   \item{dump}{The output format.}
#   \item{chopless, chopmore}{Lower and upper percentile thresholds at
#     which extreme data points are assigned zero weights.}
#   \item{maxiter}{"maximum number of REX steps". Default value is 60.}
#   \item{...}{Named argument passed to the external 'cfit' program.}
#   \item{retX}{If @TRUE, an estimate of \code{X} is returned, 
#     otherwise not.}
#   \item{cfit}{Shell command to call the 'cfit' executable.}
#   \item{verbose}{If @TRUE, verbose output is displayed, otherwise not.}
# }
#
# \value{
#   Returns a named @list structure elements:
#    \item{\code{M}}{IxN @matrix where each rows is the coordinate for
#      one of the vertices.}
#    \item{\code{X}}{(optional) the IxN @matrix \eqn{X}.}
# }
#
# \details{
#   Let \eqn{Y=(y_1, \ldots, y_I)} where \eqn{y_i=(y_{i1},\ldots,y_{iN})}
#   is an observation in \eqn{N} dimensions.  
#   Let \eqn{M=(\mu_1,\ldots,\mu_K)} be the \eqn{K}-dimensional simplex
#   where \eqn{mu_k} is a vertex in \eqn{N} dimensions.
#   Let \eqn{X=(x_1,\ldots,X_I)} where \eqn{x_i=(x_{i1},\ldots,x_{iN})}.
#   The simplex fitting algorithm decompose \eqn{Y} into:
#   \deqn{
#     Y \approx MX
#   }
#   such that \eqn{\sum_i x_{ik} = 1}.
# }
#
# \examples{@include "..\incl\cfit.Rex"}
#
# \author{
#   Algorithm and C code/binary by Pratyaksha J. Wirapati, 
#   \email{wirapati@wehi.edu.au}.
#   R wrapper by Henrik Bengtsson, \email{hb@maths.lth.se}.
# }
#
# \references{
#  [1] P. Wirapati, & T. Speed, \emph{Fitting polyhedrial cones and
#     simplices to multivariate data points}, Walter and Eliza Hall Institute 
#     of Medical Research, December 30, 2001.\cr
#  [2] P. Wirapati and T. Speed, \emph{An algorithm to fit a simplex 
#     to a set of multidimensional points}, Walter and Eliza Hall Institute 
#     of Medical Research, January 15, 2002.\cr
# }
#
#*/#########################################################################
setMethodS3("cfit", "matrix", function(y, k=ncol(y)+1, dump=1, chopless=NULL, chopmore=NULL, maxiter=NULL, ..., retX=FALSE, cfit=getOption("cfit"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments:
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'k'
  if (!is.numeric(k)) {
    throw("Argument 'k' (number of simplex vertices) is not numeric: ", k);
  }

  if (k <= 0 || k %% 1 != 0) {
    throw("Argument 'k' (number of simplex vertices) must be a positive integer: ", k);
  }

  # Argument 'cfit'
  if (is.null(cfit)) {
    cfit <- "cfit";
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the arguments passed to the external 'cfit' software
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the system call
  args <- "";
  args <- paste(args, k);

  # -M, --outM <filename>:
  #  where M is dumped [ stdout ]
  args <- paste(args, "--outM", "stdout");

  # -X, --outX <filename>:
  #  where X is dumped [ stdout ]
  if (retX) {
    fileX <- tempfile();
    on.exit(file.remove(fileX), add=TRUE);
    args <- paste(args, "--outX", fileX);
  }

  # -d, --dump <int>:
  #  dump mode (0 none, 1 final, 2 each alpha, 3 each REX)
  if (!is.numeric(dump))
    stop("Argument 'dump' is not numeric.");
  if (dump < 0 || dump > 3)
    stop(paste("Argument 'dump' must be between 0 and 3: ", dump, sep=""));
  args <- paste(args, "--dump", dump);

  # -q, --chopless <num>:
  #  assign zero weight to points if u^T y is less
  if (!is.null(chopless)) {
    if (!is.numeric(chopless))
      stop("Argument 'chopless' is not numeric.");
    args <- paste(args, "--chopless", chopless);
  }

  # -Q, --chopmore <num>:
  #  or more than the specified percentiles.
  if (!is.null(chopmore)) {
    if (!is.numeric(chopmore))
      stop("Argument 'chopmore' is not numeric.");
    args <- paste(args, "--chopmore", chopmore);
  }

  # -T, --maxiter <int>:
  #  maximum number of REX steps
  if (!is.null(maxiter)) {
    if (!is.numeric(maxiter))
      stop("Argument 'maxiter' is not numeric.");
    maxiter <- as.integer(maxiter);
    args <- paste(args, "--maxiter", maxiter);
  }

  args0 <- list(...);
  if (length(args0) > 0) {
    ndashes <- sapply(nchar(names(args0)), FUN=min, 2);
    names(args0) <- paste(c("-", "--")[ndashes], names(args0), sep="");
    args0 <- lapply(args0, FUN=function(values) {
      values <- paste(values, collapse=",");
    })
    args0 <- paste(names(args0), unlist(args0), sep=" ");
    args0 <- paste(args0, collapse=" ");
    args <- paste(args, args0);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Export data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write data to file
  infile <- tempfile();
  write.table(y, file=infile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE);
  on.exit(file.remove(infile), add=TRUE);

  fi <- file.info(infile);
  if (verbose) {
    cat("Temporary data file written:\n");
    print(fi);
    bfr <- readLines(infile);
    n <- length(bfr);
    cat("Number of lines: ", n, "\n");
    cat("First 5 lines of temporary data file:\n");
    rows <- intersect(seq(length=n), 1:5);
    print(bfr[rows]);
    cat("Last 5 lines of temporary data file:\n");
    rows <- intersect(seq(length=n), (n-4):n);
    print(bfr[rows]);
    rm(bfr, n, rows);
  }

  if (is.na(fi$size) || fi$size == 0) {
    t <- capture.output(print(fi));
    throw("Cannot fit simplex: Failed to write the temporary data file: ", t);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (FALSE) {
    path <- dirname(cfit);
    opwd <- getwd();
    on.exit(setwd(opwd), add=TRUE);
    setwd(path);
    cfit <- basename(cfit);
  }
  cmd <- paste(cfit, args, infile);
  if (verbose)
    cat(cmd, "\n");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Launch algorithm
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  con <- pipe(cmd, open="");
  if (verbose) {
    cat("Pipe connection:\n");
    str(summary(con));
  }

  # Parse results
  colClasses <- rep("double", k-1);
  M <- read.table(file=con, colClasses=colClasses, quote="", comment.char="");
  M <- as.matrix(M);

  if (retX) {
    colClasses <- rep("double", ncol(y));
    X <- read.table(fileX, colClasses=colClasses, quote="", comment.char="");
    X <- as.matrix(X);
  }

  if (dump == 2) {
    l <- list();
    for (i in 1:(nrow(M)/k)) {
      offset <- (i-1)*k;

      tmp <- M[offset+1:k,];
      class(tmp) <- "cfit";
      l <- c(l, list(tmp));
    }
    M <- l;

    if (retX) {
      l <- list();
      for (i in 1:length(M)) {
        offset <- (i-1)*nrow(y);
  
        tmp <- X[offset+1:nrow(y),];
        l <- c(l, list(tmp));
      }
      X <- l; 
    }
  } else {
    class(M) <- "cfit";
  }

  fit <- list(M=M);
  if (retX) 
    fit$X <- X;

  fit;
})


###########################################################################
# HISTORY:
# 2010-04-23
# o Added Rd help to argument 'verbose' of cfit() for the matrix class.
# 2008-02-14
# o Added more verbose output.
# o Added some validation that the temporary data was written.
# 2007-05-20
# o WORKAROUND: Now the 'cfit' executable is called by its absolute 
#   pathname. To avoid problems with spaces in the pathname, the command
#   should be put within qoutation marks (as in zzz.R). Some Unix setups
#   would not recognize the command 'cfit' even if it was in the current
#   directory, but only './cfit'.  However, Windows requires '.\\cfit'
#   in such cases.  The only cross-platform solution is to use the
#   absolute pathname.
# 2006-07-21
# o Forgot the usage in Rdoc comments.
# 2006-05-20
# o Added more Rdoc comments.
# 2006-05-07
# o Updated.
# 2003-03-03
# o Created.
###########################################################################
