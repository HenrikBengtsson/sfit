cfitTests <- function(clean=TRUE, ...) {
  # Assert that 'cfit' can be called
  cfitOptions();

  # Run example
  example("cfit", package="sfit", ask=FALSE);
  if (clean) dev.off();
} # cfitTests()

###########################################################################
# HISTORY:
# 2013-10-08
# o Added cfitTests().
###########################################################################

