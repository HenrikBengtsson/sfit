setMethodS3("lines3d", "cfit", function(object, dim=c(1,2,3), ...) {
  u <- getEdges(object);

  for (ii in seq(u)) {
    xyz <- u[[ii]][,dim];
    lines3d(xyz, ...);
  }
}) # lines3d.cfit()


###########################################################################
# HISTORY:
# 2007-06-10
# o BUG FIX: lines3d() for 'cfit' queried non-existing objects.
# 2006-05-07
# o Created.  For now, these functions are only for internal use and
#   for the examples.
###########################################################################
