# Allows conflicts. For more information, see library() and
# conflicts() in R base.
.conflicts.OK <- TRUE


.onAttach <- function(libname, pkgname) {
## .First.lib <- function(libname, pkgname) {
  pkg <- utils::packageDescription(pkgname);

  # Assert that the binaries install successfully
  if (.Platform$OS.type == "windows") {
    binfile <- "cfit.exe";
  } else {
    binfile <- "cfit";
  }
  pathname <- system.file("bin", binfile, package=pkgname);
  if (is.null(pathname)) {
    stop("Could not locate 'bin/", binfile, "' in package '", pkgname, "' (v", pkg$Version, "). It seems like the installation of the package failed.  Please report this to ", pkg$Maintainer, ".");
  }

  # cfit system command
  cfitPath <- system.file("bin", package=pkgname);
  pathname <- file.path(cfitPath, "cfit");

  # Set the default 'cfit' command (within quotation marks)
  cmd <- sprintf("\"%s\"", pathname);
  options("cfit"=cmd);

  cat(pkgname, " v", pkg$Version, " (", pkg$Date, ")",
      " successfully loaded. See ?", pkgname, " for help.\n", sep="");
}


###########################################################################
# HISTORY:
# 2007-06-12
# o Replaced .First.lib() with .onAttach().
# 2007-05-20
# o WORKAROUND: Put quotation marks around default 'cfit' command.
###########################################################################
