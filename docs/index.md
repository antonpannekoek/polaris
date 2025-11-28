# Polaris - VSTPol polarization reduction package

Polaris is a prototype reduction pipeline for VSTPol data, processing data after the last stage of the general AstroWISE pipeline for OmegaCam, to produce Stokes-vector images of polarization observations.

VSTPol[^1][^2] is an instrument consistenting of the VST, the 2.6 meter VLT (Very Large Telescope) Survey Telescope, equipped with a rotatable polarizer plate; it uses OmegaCam as the imaging backend, a wide-field imager (about 1 degree diameter) consistenting of an 8 by 4 mosaic of 2K by 4K CCDs, totalling 256 Megapixels.

This pipeline demodulates a set of images taken at different polarizer angles and turns them into images corresponding to a set of Stokes-vectors I, Q and U; a step called demodulation. The input images are taken at the same sky position, three or four times with differerent polarizer angles: 0, 60 and 120 degrees, or 0, 45, 90 and 135 degrees. A set of two images taken at angles theta and theta + 90 are also possible, but this will not disentangle the Stokes Q and U parameters, and are only meant to obtain the polarization degree at a certain angle on the sky.

The pipeline contains an optional correction step for instrumental (induced) polarization because of potential shifts between the individual set of images; the pointing of the VST is very accurate, but small, subpixel shifts may occur between consecutive images, which can result in induced polarization in the demodulated Q and U images. There are a number of ways to correct for this, and several of these options have been implemented in the pipeline.

The details of the reduction are described in Van Vorstenbosch, A.T.P., 2019, master thesis[^3].

## Specific documents

Details on the pipeline can be found in the [pipeline](pipeline.md) document.

In addition, a very simple simulation module also exists. This can create multi-HDU (header-data unit) FITS images with random stars in a magnitude range, and additionally use a stellar catalog as input, or create such a catalog as output. This catalog can then be used in the correction step in the reduction pipeline. See the [simulation](simulation.md) document for more details.

For a quick overview of practical use of the package, see the [usage](usage.md) document.

For details on the software itself, see the [implementation document](implementation.md).

### Dependencies

The software has various existing packages as dependencies: numpy, scipy, pandas, astropy, sep, astroalign, duckdb. Several of these may be replaced:

- `duckdb` is for convenience, but could be replaced by the built-in `sqlite` module and in a more final version by another SQL library or SQLAlchemy, to connect directly to e.g. the GAIA database.

- AstroAlign may be replaced by Scikit-learn, which is a dependency of AstroAlign, or similar functionality from SciPy (providing an interpolated shift between 2D images). The AstroAlign functionality used here is mainly a linear shift correction: matching between images is already done through matching source lists, and it is not expected that transformations other than a shift (e.g., rotation or shear) are necessary.

## References

[^1]: Smette, A., et al, 2018, " VST: The First Large Survey Telescope for Optical Polarimetry ", https://zenodo.org/records/1304780
[^2]: Schipani, P., et al, 2025, "VSTPOL: making the VST a large survey telescope for optical polarimetry", https://arxiv.org/abs/2505.07113
[^3]: Van Vorstenbosch, A.T.P., 2019, "Theoretical performance of a wide field polarimetry survey usinga polarimetric upgrade of OmegaPOL on the VLT Survey Telescope", master thesis, University of Leiden.
