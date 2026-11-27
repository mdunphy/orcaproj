# orcaproj

## Overview

`orcaproj` is a C library that implements the ORCA05 tripolar-grid coordinate
transformations [1]. A Python wrapper exposes the library through
`ctypes`.

This is an initial release of the library. Packaging, cross-platform support, optimization, and other improvements are still in progress. The API should not be considered stable at this stage.

## Requirements

For the C library:

- GCC or a compatible C compiler
- GNU Make

For the Python wrapper:

- Python 3.9 or newer
- NumPy

## Build and install

From the repository root:

```bash
make
python -m pip install -e .
```

This builds the native library and installs the Python wrapper in editable
mode.

## Usage

```python
import orcaproj

# Forward transform
i, j, status = orcaproj.lonlat_to_ij(300.0, 45.0)

# Inverse transform
lon, lat, status = orcaproj.ij_to_lonlat(75.0, 495.0)
```

Forward arguments `lon` and `lat` are specified in decimal degrees.
The inverse transform accepts fractional ORCA05 grid coordinates i in [0, 720] and j in [1, 510.5].

## References

[1] Dunphy, M. (2027). A projection for the ORCA05 tripolar coordinate
system. *Ocean Modelling*, 205, 102862.
<https://doi.org/10.1016/j.ocemod.2026.102862>
