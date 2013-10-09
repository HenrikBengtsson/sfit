# Allows conflicts. For more information, see library() and
# conflicts() in R base.
.conflicts.OK <- TRUE

.onLoad <- function(libname, pkgname) {
  findCfitBinary(pkgname=pkgname);
}

.onAttach <- function(libname, pkgname) {
  pkg <- utils::packageDescription(pkgname);
  pkgStartupMessage(pkgname, " v", pkg$Version, " (", pkg$Date, ")",
      " successfully loaded. See ?", pkgname, " for help.\n", sep="");
}


###########################################################################
# HISTORY:
# 2012-10-07
# o Now utilizing pkgStartupMessage() of R.methodsS3.
# 2012-03-23
# o Now .onAttach() uses packageStartupMessage() instead of cat().
# 2011-05-15
# o Now .onAttach() utilizes new findCfitBinary() to locate and set
#   the 'cfit' binary option.
# 2007-06-12
# o Replaced .First.lib() with .onAttach().
# 2007-05-20
# o WORKAROUND: Put quotation marks around default 'cfit' command.
###########################################################################
