# Version (development version)

 * ...
 

# Version 0.3.2 [2022-11-07]

## Bug Fixes

 * Native code would try to read from file even where there was
   nothing to read.  This was caught by the `-Wunused-result` compiler
   option.
   
 * Native code had a non-guarded else clause, i.e. it was not
   surrounded by `{ ... }`.
 

# Version 0.3.1 [2017-03-22]

## Software Quality

 * Now explicitly importing **graphics**, **grDevices**, and **utils**
   functions.
 
 * Removed some clang compiler warnings and a harmless bug in a part
   of the C code that was never used by this package.
 
 
# Version 0.3.0 [2012-10-07]

## New Features

 * Now `sfit::cfit()` works (without having to attach the package).
 
 * Now `library("sfit", quietly = TRUE)` attaches the package
   completely silently without any output message.

## Refactoring

 * Now only importing **R.methodsS3** (no longer attaching it).

 * Updated details on authors.
 
## Software Quality

 * Added `cfitTests()` to test the installation.
 
 * Added several system tests.

## Bug Fixes

 * Now declaring all S3 methods in NAMESPACE.
 
 
# Version 0.2.2 [2012-07-16]

## Software Quality

 * Now declaring S3 methods in NAMESPACE.

 * Added an `Authors@R` field to DESCRIPTION.

 * CRAN POLICY: Now all Rd example lines are at most 100 characters
   long.
 
 
# Version 0.2.1 [2012-03-23]

## New Features

 * Now package uses `packageStartupMessage()` instead of `cat()` when
   loaded.
 
 
# Version 0.2.0 [2011-05-15]
 
## New Features

 * If argument `cfit` of `cfit()` is NULL, which happens if the `cfit`
   option is not set, then `cfit()` will generate a warning explaining
   that it will rely on the operating system to find the `cfit`
   binary, which may not work.

## Software Quality

 * Removed warning reporting on an uninitialized `f` when building
   executable `cli`.  This warning was unharmful and the "fix" makes
   no difference.

## Bug Fixes

 * BUILD BUG FIX: The makefiles for building package binaries assumed
   that a directory `../incl/bin/` existed.  Starting with R v2.13.0,
   empty ("skeleton") directories are automatically dropped (already
   when building the `*.tar.gz`), which caused error `../inst/bin/: No
   such file or directory` during the binary builds.
 
 * If the `cfit` binary had a comma in its path, `cfit()` would, if on
   Windows, throw an error reporting that the command "is not
   recognized as an internal or external command, operable program or
   batch file". Now such commas are escaped if on Windows.  See source
   code of `findCfitBinary()` for all details.
 
 
# Version 0.1.9 [2010-04-23]

## Documentation

 * Added Rd help to argument `verbose` of `cfit()` for the `matrix`
   class.
 
 
# Version 0.1.8 [2008-05-09]

## Deprecated and Defunct

 * Removed obsolete `SaveImage` from DESCRIPTION.
 
 
# Version 0.1.7 [2008-08-31]

## New Features

 * Added `drawApex()` and `radials()` for class `cfit`.

## Refactoring

 * Renames old HISTORY file to NEWS.
 
 
# Version 0.1.6 [2008-02-14]

## Software Quality

 * `cfit.matrix()` now tests if the temporary data file was
   written/exists.

## Refactoring

 * Replaced dependency on **R.oo** with dependency on **R.methodsS3**.
 
 
# Version 0.1.5 [2007-06-12]

## Software Quality

 * Added a NAMESPACE to the package.
 
 
# Version 0.1.4 [2007-06-10]
 
## Software Quality

 * Package pass `R CMD check` R v2.6.0.

## Bug Fixes

 * Internal `lines3d()` for `cfit` queried non-existing objects.
 
 
# Version 0.1.3 [2007-05-20]

## Software Quality

 * WORKAROUND: On some Unix systems, `cfit()` would give output `sh:
   cfit: command not found`.  This is because the sh shell has not be
   setup to identify executable in the current directory, i.e. they
   have to be called with `./cfit`, but that does not work on Windows.
   Instead, the `cfit()` is now calling `pipe()` with the absolute
   pathname to the `cfit` executable within quotation marks (to avoid
   problems with spaces) in the command string.
 
 
# Version 0.1.2 [2006-07-21]

## Documentation

 * Updated the help pages.
 
 
# Version 0.1.1 [2006-05-20]

## Software Quality

 * First version of the package where the `cfit` binary is
   automatically build upon `R CMD INSTALL`.  It is still no shared
   library that is build, i.e. the approach to dump data to file, call
   `cfit`, and then let R parse the cfit result is still used.  But at
   least it should not (in theory) install and run on all platforms.
 
 
# Version 0.1.0 [2006-05-07]
 
## New release

 * Created.  Had a similar version in 2003.

