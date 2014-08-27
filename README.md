weatherGen Package
=========================

Jeffrey D. Walker, PhD

This package will contain functions for generating synthetic daily timeseries of precipitation and air temperature under future climate change scenarios. These timeseries can then be used for performing a climate stress test.

This package is under heavy development and should not be used until its official release.

## Installation

This package must be installed from github using devtools. Note that the vignettes will fail at the moment because they rely on external data.

```R
library(devtools)
install_github('walkerjeffd/r-weathergen', build_vignettes=FALSE)
library(weatherGen)
```

## Building

To build from source, run:

```R
library(devtools)
document()
build()
```

To build only vignettes, run:

```R
library(devtools)
build_vignettes()
```