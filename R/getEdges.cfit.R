setMethodS3("getEdges", "cfit", function(object, ...) {
  u <- list();
  for (ii in 1:(nrow(object)-1)) {
    for (jj in (ii+1):nrow(object)) {
      idx <- c(ii,jj);;
      u <- c(u, list(object[idx,]));
    }
  }

  u;
}, private=TRUE) # getEdges.cfit()


###########################################################################
# HISTORY:
# 2007-06-10
# o BUG FIX: lines3d() for 'cfit' queried non-existing objects.
# 2006-05-07
# o Created.  For now, these functions are only for internal use and
#   for the examples.
###########################################################################
