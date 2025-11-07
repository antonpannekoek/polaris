# VSTPOL pipeline

Note: this software is currently very much pre-alpha software!

## VSTPOL

VSTPOL is the polarization variant of the Very Large Telescope (VLT)
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

Installation is straightforward. Create a Python virtual environment and activate it. Then, from the downloaded or cloned project, run

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

This will install some of the packages used for development (formatters, linters etc).


## Some example usage

### Simulate a set of observations

Use the configuration file `sim.cfg`, in TOML format, to set any parameters.

The default will create a 4 by 2 mosaic of 2K by 4K pixel chips, with random stars 25,000 random stars, plus some stars added from a selected set of GAIA stars (the `gaia.duckdb` is added to this repository for convenience). A small shift of 0.2 pixels is given for the x-position of the 60-degree-polarizer frame (this applies to all chips in the mosaic simultaneously).

Running

```
python -m polaris.simulate sim.cfg
```

will create three FITS files, one for each measurement angle. The FITS files will contain several extensions, each extension matching 1 chip; the primary extension is empty.

### Reduce a set of frames

If the default settings are used, the simulation will create a set of `polaris-deg*_0.fits` files. These can then be reduced to their I, Q and U (mosaic) frames, as well as the relative frames, and the corrections for shifts are estimated, by using

```
python -m polaris.frames polaris-deg0_0.fits polaris-deg60_0.fits polaris-deg120_0.fits --config polaris.cfg
```

Note that both the simulation and reduction will print various simple debugging output; this can be ignored. The last print outputs from the reduction will show the estimated scale factors of the Q & U polarization with respect to the I gradient, alpha.
