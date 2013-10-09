#########################################################################/**
# @RdocFunction findCfitBinary
#
# @title "Located the cfit binary"
#
# \description{
#   @get "title" and sets option 'cfit'.
# }
#
# @synopsis
#
# \arguments{
#   \item{pkgname}{A @character string specifying the name of the
#     package to be searched.}
#   \item{quote}{If @FALSE, the returned string is not quoted.}
#   \item{force}{If @FALSE and a previously located binary exists,
#     then it is used.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) a @character string.
# }
#
# @keyword internal
#*/#########################################################################
findCfitBinary <- function(pkgname="sfit", quote=TRUE, force=FALSE, ...) {
  # Already located?
  cmd <- getOption("cfit");
  if (!force && !is.null(cmd)) {
    if (!quote) {
      cmd <- gsub("^\"", "", cmd);
      cmd <- gsub("\"$", "", cmd);
    }
    return(invisible(cmd));
  }

  # (1) Assert that the binaries exists
  path <- system.file("bin", package=pkgname);
  if (is.null(path)) {
    pkg <- utils::packageDescription(pkgname);
    stop("The 'bin/' directory does not exist in package '", pkgname, "' (v", pkg$Version, "). It seems like the installation of the package failed.  Please report this to ", pkg$Maintainer, ".");
  }

  # (2) Assert that the binaries install successfully
  onWindows <- (.Platform$OS.type == "windows");
  if (onWindows) {
    filename <- "cfit.exe";
  } else {
    filename <- "cfit";
  }
  pathname <- file.path(path, filename);
  if (is.null(pathname)) {
    stop("Could not file '", filename, "' in directory '", path, "' of package '", pkgname, "' (v", pkg$Version, "). It seems like the installation of the package failed.  Please report this to ", pkg$Maintainer, ".");
  }

  # cfit system command (independent of platform/OS)
  pathname <- file.path(path, "cfit");

  # Special cases?
  if (onWindows) {
    # If on Windows, and there are commas in the pathname, those
    # needs to be escaped.  Details: cfit.matrix() uses pipe(cmd)
    # which in turn uses shell().  On Windows commas in shell
    # commands has a special meaning and needs to be escaped
    # using '^' [1].
    # REFERENCES:
    # [1] Command shell overview, Microsoft, 2011-05-15
    #     http://technet.microsoft.com/en-us/library/cc737438(WS.10).aspx
    pathname <- gsub(",", "^,", pathname, fixed=TRUE);
  }


  # Set the default 'cfit' command (within quotation marks)
  cmd <- sprintf("\"%s\"", pathname);
  options("cfit"=cmd);

  if (!quote) {
    cmd <- gsub("^\"", "", cmd);
    cmd <- gsub("\"$", "", cmd);
  }

  # Return the located pathname
  invisible(cmd);
} # findCfitBinary()


###########################################################################
# HISTORY:
# 2013-10-08
# o Turned into a plain function.
# o Added argument 'quote' to findCfitBinary().
# 2011-05-15
# o Now commas in the located pathname are escape if on Windows.
# o Created from zzz.R.  Previous history below.
# 2007-05-20
# o WORKAROUND: Put quotation marks around default 'cfit' command.
###########################################################################
