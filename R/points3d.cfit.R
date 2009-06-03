setMethodS3("points3d", "cfit", function(object, dim=c(1,2,3), ...) {
  xyz <- object[,dim];
  points3d(xyz, ...);
}) # points3d.cfit()



###########################################################################
# HISTORY:
# 2007-06-10
# o BUG FIX: lines3d() for 'cfit' queried non-existing objects.
# 2006-05-07
# o Created.  For now, these functions are only for internal use and
#   for the examples.
###########################################################################
