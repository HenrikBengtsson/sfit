## Installing

The **sfit** package is _neither_ on CRAN nor on Bioconductor.  To install it, it can be installed from [HenrikBengtsson's personal DRAT](https://github.com/HenrikBengtsson/drat) by calling:

```r
install.packages("sfit", repos = "https://henrikbengtsson.github.io/drat")
```

This should work on Linux and MS Windows.  On macOS, you need to have compiler tools installed.


### Installing from source

Alternatively, you can install it directly from source.  Because this package also compiles native code, Windows users need to have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed and macOS users need to have [Xcode](https://developer.apple.com/xcode/) installed.  To install the latest stable version from source, call:
```r
remotes::install_github("HenrikBengtsson/sfit", ref="master")
```
To install the pre-release version that is available in Git branch `develop` on GitHub, use:
```r
remotes::install_github("HenrikBengtsson/sfit", ref="develop")
```
