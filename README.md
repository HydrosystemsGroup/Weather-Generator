weathergen R Package
====================

[Jeffrey D. Walker, PhD](http://walkerjeff.com)

This package will contain functions for generating synthetic daily timeseries of precipitation and air temperature under future climate change scenarios. These timeseries can then be used for performing a climate stress test.

This package is under heavy development and should not be used until its official release.

## Installation

Before installing the `weathergen` package, you must install the `hydromad` package ([homepage](http://hydromad.catchment.org/)) using the following commands:

```r
install.packages(c("zoo", "latticeExtra", "polynom", "car", "Hmisc","reshape"))
install.packages("hydromad", repos="http://hydromad.catchment.org")
```

The `weathergen` package can then be installed from this github repo using the `install_github()` function from the `devtools` package.

```R
devtools::install_github('walkerjeffd/weathergen')
library(weathergen)
```

## Building From Source

To build the package from source:

```R
devtools::document()
devtools::build()
```

To build only vignettes, run:

```R
devtools::build_vignettes()
```