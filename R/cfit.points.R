setMethodS3("points", "cfit", function(x, dim=c(1,2), ...) {
  # To please R CMD check
  object <- x;

  xyz <- object[,dim];
  points(xyz, ...);
}) # points.cfit()


###########################################################################
# HISTORY:
# 2007-06-10
# o BUG FIX: lines3d() for 'cfit' queried non-existing objects.
# 2006-05-07
# o Created.  For now, these functions are only for internal use and
#   for the examples.
###########################################################################
