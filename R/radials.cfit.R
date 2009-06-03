setMethodS3("radials", "cfit", function(fit, ...) {
  usr <- matrix(par("usr"), nrow=2);
  stretch <- max(apply(usr, MARGIN=2, FUN=diff));
  X <- t(fit);
  xy <- stretch*(X[,c(3,1,2)] - X[,1]) + X[,1];
  xy <- t(xy);
  lines(xy, ...);
})


###########################################################################
# HISTORY:
# 2008-08-31
# o Added radials() for class cfit to be compatible with the newer 
#   expectile package.
###########################################################################
