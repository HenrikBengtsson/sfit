


# sfit: Multidimensional Simplex Fitting


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


## Contributions

This Git repository uses the [Git Flow](https://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/sfit/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/sfit) branch contains the code of the latest release.

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [sfit repository](https://github.com/HenrikBengtsson/sfit).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/sfit">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/sfit">AppVeyor CI</a> when the PR is submitted.

We abide to the [Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/) of Contributor Covenant.


## Software status

| Resource      | GitHub        | GitHub Actions      | Travis CI       | AppVeyor CI      |
| ------------- | ------------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   |  |        | <a href="https://travis-ci.org/HenrikBengtsson/sfit"><img src="https://travis-ci.org/HenrikBengtsson/sfit.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/sfit"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/sfit?svg=true" alt="Build status"></a> |
| Test coverage |                     |                     | <a href="https://codecov.io/gh/HenrikBengtsson/sfit"><img src="https://codecov.io/gh/HenrikBengtsson/sfit/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |
