cfitOptions <- function() {
  bin <- findCfitBinary(quote=FALSE);
  suppressWarnings({
    bfr <- system2(bin, stderr=TRUE)
  });
  bfr <- paste(bfr, collapse="\n");
  message(bfr);
} # cfitOptions()


###########################################################################
# HISTORY:
# 2013-10-08
# o Added cfitOptions().
###########################################################################
