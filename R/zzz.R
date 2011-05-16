# Allows conflicts. For more information, see library() and
# conflicts() in R base.
.conflicts.OK <- TRUE


.onAttach <- function(libname, pkgname) {
## .First.lib <- function(libname, pkgname) {
  pkg <- utils::packageDescription(pkgname);

  findCfitBinary(pkgname=pkgname);

  cat(pkgname, " v", pkg$Version, " (", pkg$Date, ")",
      " successfully loaded. See ?", pkgname, " for help.\n", sep="");
}


###########################################################################
# HISTORY:
# 2011-05-15
# o Now .onAttach() utilizes new findCfitBinary() to locate and set
#   the 'cfit' binary option.
# 2007-06-12
# o Replaced .First.lib() with .onAttach().
# 2007-05-20
# o WORKAROUND: Put quotation marks around default 'cfit' command.
###########################################################################
