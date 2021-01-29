# nemTar
A method based on graphical model for cancer regulatory network inference and prioritization of potential therapeutic targets
## Quick Installation

### This package is available under R(>=3.5.0, 3.6.x is preferred).

If you have installed `devtools` package, you only need to call `install_github` function in `devtools` to install `nemTar`.

```
Please refer to the following commands to install **nemTar**

library("devtools")

devtools::install_github("CityUHK-CompBio/nemTar", dependencies=TRUE)

```

## Dependency

**PAGnet** requires the following R/Bioconductor packages for its normal function:

- nem
- dplyr

## Quick Start to Apply **nemTar**

A toy example on simulated network and two case studies involved in GC and CRC could be found in our [vignette] (https://github.com/CityUHK-CompBio/PAGnet/blob/master/vignettes/PAGnet.pdf)

**PAGnet** also provides a local shiny interface to perform MRA in PAGnet with signature gene sets in Gene Ontology (GO) and KEGG databases obtained from Pseudomonas Genome DB.
Call local shiny interface:

### For more details, pelease refer to our [vignette](https://github.com/CityUHK-CompBio/PAGnet/blob/master/vignettes/PAGnet.pdf)

## Getting help

Should you have any questions about this package, you can either email to the developers listed in the DESCRIPTION part of this package or create an issue in the [issue part](https://github.com/CityUHK-CompBio/PAGnet/issues).
