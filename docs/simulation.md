# Simulating a mosaic of polarization images

The `polaris.simulation` module can create a set of multi-extension FITS files, one file per polarizer angle, and each file containing a set of FITS files mapping to a mosaic on individual chips.

The simulation is fairly simple, using `photutils` to simulate a field of stars within a given magnitude range. A first pass adds unpolarized stars. A second pass adds a small set of polarized stars, with random values for U and Q. Before these are added to the polarization images, the Stokes flux values are modulated so they map to a polarizer angle, which is what is actually observed.

The mosaic is optional: it is configurable and a 1 x 1 mosaic consisting of just one image per angle is entirely possible. A shift between images (x and y separately) can also be added, in units of pixels.

A database with existing stars can be used to add selected stars to the image (on top of random simulated stars). Alternatively, for testing purposes, the database entries can be created instead from a random selection of simulated stars, which can then be used in the pipeline as reference stars for correcting instrumental polarization.

The database should contain at least a table `catalog`, with the following three columns minimum: `ra`, `dec` and `mag`. Each column should be a (single or double precision) floating point number.


Below is a list of configuration options. These are given in as a file in TOML format, subdivided by sections.

## Configuration options

- `[logging]`
```
level = "info"
```

Set the logging output to one of `"warning"`, `"info"` or `"debug"`. The default is `"info"`. Logging is currently rather minimal, but it's easy to add further logging output.


- `[simulation]`

```
random_seed = -1
```

Set the seed for NumPy's random number generator. A negative integer means unset, and NumPy will pick a random seed.

```
t_exp = 300.0
```

The exposure time in seconds


```
nstars = 25_000
```

The number of stars to simulate. This is the number of stars across the entire mosaic, so the field will be relatively more dense if fewer chips are used.


```
magrange = [8.0, 22.0]
```

The magnitude range of the simulated stars. The distribution follows a simple logarithmic formula, with more stars at the fainter magnitudes.


```
zeropoint = 25.0
```

The zeropoint of the instrument.

```
skylevel = 2
```

The skylevel in ADU (analog to digital unit) or counts per second.


```
poisson_noise = true
```

Set this to False is no (Poisson) noise is wanted (probably only for testing purposes).


- Section `[simulation.catalog]`

```
action = "use"
```

One of `"use"` or `"create"`.

With "use", this will read an existing catalog and add stars found in the catalog to the simulated stars.

With "create", only randomly simulated stars are used, and a selected fraction is stored in the database.


```
dbname = "sim.duckdb"
```

Currently, the file name of a DuckDB database, that either contains the stars to be used (action = "use") or where to write a selection of stars to (action = "create"). This file will be automatically overwritten when the action is "create".


```
fraction = 0.1
```
The fraction of (unpolarized) simulated stars to save to the database, in case of action = "create".


- Section `[simulation.psf]`

```
fwhm = 4.0
```
The full-width half-maximum of the stellar point-spread function (PSF).
The simulation uses the Moffat function for the PSF.



```
fwhm_std = 0.2
```

A standard deviation of the above FWHM of the PSF.


```
beta = 3.
```

The `beta` ("wings") factor of the Moffat PSF function. `alpha` is derived from the `fwhm` value above.


- Section `[simulation.ccd]`

This is the geometry for a single chip.

```
shape = [2048, 4096]
```

The x by y size of the CCD in pixels


```
pixelsize = 0.2
```

The pixel size in arcseconds.

```
ron = 5.0
```

The read-out noise of the CCD.

- Section `[simulation.mosaic]`

The layout of the chips in the full mosaic

```
layout = [1, 2]
```

The grid of individual CCDs. n by m CCDs (x and y direction, respectively).

```
gaps = [30, 40]
```

The gaps between the CCDs in the n and m (x and y) directions, in pixels


- Setion `[simulation.sky]`

```
ra = 123.0
dec = 0.45
```

The central pointing of the mosaic. Use degrees, not sexagesimal notation.


- Section `[simulation.polarisation]`

```
angles = [0, 60, 120]
```

The polarizer angles for each individual observation. These are saved with the keyword `polangle` in the FITS header.


```
nstars = 10
magrange = [15, 20]
```

The number and magnitude range of the polarized stars.


```
qrange = [0.05, 0.3]
urange = [0.05, 0.3]
```

The degree of polarization for Q and U. This is a range (here, between 5 and 30 percent) with a random degree used inside that range for each of the polarized stars.


- Section `[simulation.shifts]`


# shift in x or y for a frame at _<deg> angle, in pixels
```
x_60 = 0.0
x_60 = 0.0
x_120 = 0.2
y_120 = 0.0
```

The pointing shift in (fractional) pixels in x and y for each of the angles, relative to the 0-angle frame.


- Section `[simulation.region]`

This saves the polarized stars to the given region file, which can be used as a reference to find the polarized stars in the Q and U images.
```
regfile = "ds9.reg"
```
