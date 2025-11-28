# VSTPOL pipeline

Note: this software is currently very much alpha software!

Polaris is a prototype for the polarization reduction step for VSTPol.

## VSTPOL

VSTPOL is the polarization variant of the VLT (Very Large Telescope)
Survey Telescope, a 2.6 meter survey telescope located at the Paranal
VLT site in Chile. Attached to it is the 1-square degree OmegaCAM
instrument, consistent of 32 CCDs with a total of 268 million (16K
times 16K) pixels. A rotatable polarization filter is placed in a
spare slot in the light pathway, and the combination of different
rotations allows for measuring the various I, Q, U and V Stokes
parameters of the polarization of imaged sources.

## Pipeline

The pipeline handles the input images for different rotations,
combining them such that the Stokes parameters are retrieved, through
demodulation of the input frames. Since these frames are not obtained
simultaneously, but in sequence after each other, corrections are then
applied, to handle shifts and differences in the point spread
functions between the frames. This is largely done through use of
stars in the field, which for a large fraction are taken to be
unpolarized. Selection of these stars, from the GAIA catalog, can
strengthen this assumption, for example by using stars relatively in
the Solar neighbourhood, which avoids polarization due to gas and
dust; and by selecting stellar types that are known to be little to
un-polarized.

## Installation

Installation is straightforward. Create a Python virtual environment
and activate it. Then, from the downloaded or cloned project, run

```
pip install .
```

or install it from its GitHub repository directly 

```
pip install git+https://github.com/evertrol/polaris.git
```

to install Polaris in that virtual environment.

To install a development environment, use

```
pip install -e '.[dev]' .
```

This will install some of the packages used for development
(formatters, linters etc).


## Documentation

Documentation can be find inside the [docs](docs/index.md) directory.

## License

Polaris is copyright 2025 University of Amsterdam, and licensed under
the MIT license. See LICENSE for the details.
