# Support for Albers Equal Area projection in MET/METPlus

The current (September 2025) version of MET and METPlus does not support the
Albers Equal Area projection (hereafter simply referred to as the Albers
projection). The Albers projection is a standard map projection for Australia,
and is used in the Bureau of Meteorology's IMPROVER software. There is a need
within the Bureau to have MET and METPlus support data on the Albers
projection.

This document is structured into following sections:

* Installation on a Bureau Linux computer
* Data format requirements for data on the Albers projection
* An example run of METPlus with data on the Albers projection

## Installation

To support Albers projection data in METPlus, the underlying MET code (which is
written in C++) needs to be modified and recompiled. No modifications for
METPlus are required.

### Preparing the environment for compiling MET/METPlus

Before compiling MET, create a conda environment with the various dependencies
needed by MET and METPlus. The following environment should suffice:

```
name: metenv
channels:
    - conda-forge
dependencies:
    - automake=1.16.5
    - python
    - netcdf4
    - netcdf-cxx4
    - gcc
    - gxx
    - gsl
    - proj
    - hdf4
    - hdf5
    - g2clib
    - nceplibs-bufr
    - hdfeos2
    - cairo
    - freetype
    - metplus
```

If the environment is saved in a file named `metenv.yml`, the following
commands will build and activate the environment. The `conda deactivate` is
simply in case another environment is currently loaded.

```
conda deactivate
conda env create -f metenv.yml
conda activate metenv
```

### Checkout the new MET code and compile it
