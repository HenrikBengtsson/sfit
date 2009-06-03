The interface from R is very ad hoc; it dumps all data to file(s), calls 'cfit' with the correct parameters via system(), and parses the output files from 'cfit'.  Ideally, we would link the cfit code to R, but that is for the future.

The 'cfit' binary is automatically built by R CMD INSTALL (when a package is installed) and then moved to the ../inst/bin/ directory which will be copied to the R package directory.

