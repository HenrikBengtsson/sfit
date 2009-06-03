setMethodS3("drawApex", "cfit", function(fit, ...) {
  X <- t(fit);
  X <- X[,1];
  X <- t(X);
  points(X, ...);
})

###########################################################################
# HISTORY:
# 2008-08-31
# o Added drawApex() for class cfit to be compatible with the newer 
#   expectile package.
###########################################################################
