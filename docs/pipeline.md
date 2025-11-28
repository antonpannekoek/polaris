# Overview of the pipeline module

## Terminology

Practical terms used in this section:

- polarizer or polarization frames: the actual observed frames (2D images), containing the observations for different polarizer angles

- Stokes frames: the frames derived from the polarization frames, that contain the values of the I, Q and U Stokes parameter values (as ADUs / fluxes) for objects in the frames. Note that circular polarization, V, is not used.

- demodulation: the process of adding the polarization frames with various weights, to obtain the stokes frames. Since the polarization frames are often taken at 3 or 4 angles, and the Stokes frames come in sets of 3, the relevant weights are often given in a demodulation matrix (3x3 or 3x4). The circular polarization can't be obtained this way, only linear polarization.

- induced polarization: instrumental polarization that is caused by small offsets between the individual polarization frames. The result is that the addition during demodulation is not perfect, resulting in induced polarization in the Stokes frames. This can be corrected for in various ways.


## Procedure

The pipeline `run` function is the entry point of the actual pipeline, and takes a set of polarization frames as input, as well as a configuration dictionary. It returns a set of Stokes frames. (The input and output frame sets are are Python dictionaries, with the keys being the polarizer angles and the Stokes parameters, respectively, and the values the actual frames.)

Each frame contains the data as a two-dimensional NumPy array in an `image` attribute, header information in the form of a FITS header as a `header` attribute, and one or more additional attributes such as the polarizer angle or Stokes parameter.


The function tests if the correction method given in the configuration dict is valid, then proceeds as follows:

- if the correction method is `"interpolation"`, it will register the polarization frames and align them to the 0-angle polarization frame, to remove any shifts between the individual frames. See the [image registration](#image-registration) section below for the details.

- the input frames are demodulated using the relevant matrix, yielding a set of Stokes frames. See the demodulation section below for the details.

- if the correction method is `"induced"`, the Stokes frames are corrected for induced polarization using a least-squares fit between the gradient of the Stokes I frame and the Q and U frames. See the [induced polarization correction](#induced-polarization-correction) section below for details.

- if the correction method is `"shift-induced"`, The Stokes frames are corrected for induced polarization, by calculating the shift between the input polarization frames and applying the correction from the gradient of the Stokes I frame. See the [induced polarization correction with shift](#induced-polarization-correction-with-shift) section below for details.

Note that all the correction methods are exclusive of each other: only one can be applied. `"none"` is also a valid option, to apply no correction.


### Image registration

For the image registration, the set of frames is looped over to correct each frame against the frame with polarizer angle 0.

First, for both the 0-angle frame and the current frame, sources in the image are detected using `sep`, a Python package that uses the SExtractor functionality through a Python interface (it uses the routines as a C-extension; it does not rely on SExtractor being actually available). This detects all sources in both images, retaining only those with a SExtractor flag of 0 (non-zero flags may not be problematic, but the majority of objects likely have a flag of 0, keeping enough objects to match the images).

The positions of the detected objects are converted to RA and Dec using the WCS supplied by the frames header information, but also keeping the (x, y) pixel coordinates. The RA and Dec are then used to match sources between the two images using `astropy.coordinates.match_coordinates_sky`.

The shift between the two images is then estimated using AstroAlign's `estimate_transform`, which is a re-import of `skimage`'s `estimate_transform. It uses a Euclidean transform, which is a rigid transform: only a rotation followed by a shift. Alternatively, a simpler and perhaps better approach is to estimate an x and y pixel shift by subtracting the two coordinate sets and calculate the medians of the x and y shifts, as it is expected that the offsets between the frames are minimal and the rotation between frames is neglible.

The calculated transformation is then applied to the second image, to align it with the first (0-angle) image, using AstroAlign's `apply_transform`, which uses `skimage`'s `warp` function under the hood, with bi-cubic interpolation.

The header of resulting transformed image is updated with the WCS of the 0-angle frame, since these should now be the same.


### Induced polarization correction

The correction for induced polarization follows the procedure described in Van Vorstenbosch, A.T.P., 2019.

To first order, the polarization caused by small misalignments between the individual polarization frames is proportional to the gradient of the Stokes I image. In particular, this is for the Q and U polarization frames. The misalignments are assumed to be shifts only, relatively small, and independent in x and y (as these are orthogonal), for this approximation to work. This gives

Q(x, y) = alpha * dI(x, y)/dx + beta * dI(x, y)/dy

and

U(x, y) = gamma * dI(x, y)/dx + delta * dI(x, y)/dy

with I(x, y), Q(x, y) and U(x, y) the Stokes vector images; dI/dx and dI/dy the derivatives in the x and y direction, respectively; and alpha, beta, gamma and delta the proportionality factors to be fitted for.

Note that there is a correspondence between alpha and gamma, and between beta and delta, which the pipeline currently does not take into account, allowing for two (simpler) independent least-squares fits.

The derivates (gradients) are calculated using `numpy.gradient`, and the least-squares solution is obtained using `numpy.linalg.lstsq`.


#### Using selected stars only

Before the least-squares fit is applied, optional masking of the images (including the gradient image) is performed, as follows:

If a catalogue is provided (in the form of a database file name), stars are selected within the field-of-view of the frame(s). These stars are assumed to be a selection of unpolarized stars, to be used for calibrating the induced polarization. Sources are then detected for the Stokes I frame only, in the same way as for the image registration. This source lists is then filtered to retain only sources that correspond to stars from the catalog. A square area around each of these stars is then selected in the frame, and only these pixels are kept (overlapping squares only keep a pixel once), by creating a boolean mask. This mask is then used to index the various Stokes images above: the gradient-I image, the Q and U image. (Source detection naturally has to be done in the I frame, since the Q and U images may not show any images and the gradient image will not show actual stellar sources.)

This way, the induced polarization correction is calculated only for the combined area around these stars; with the assumption that these stars are unpolarized, the Q and U polarization detected in the combined area is then only the result of induced polarization, avoiding the inclusion of actually polarized sources in the least-squares fit.

Note that in the pipeline, providing a (valid) catalog name sets this option; if the catalog name is empty, this masking step is not performed and full images are used for the least-squares fit.



### Induced polarization correction with shift

A more direct method of correcting for induced polarization relies on calculating the shifts between the polarization frames, then multiplying these shifts with a known constant factor to the gradient of the Stokes I image, to obtain the induced Stokes Q and U images.

The calculation of the shift between the two images is done similarly to that in the registration correction step, except that the transform is replaced by calculating the median of the difference in x and y positions between the two source lists.

The x and y shifts are then applied (with their constant factor) to the x and y gradients of the Stokes I frame, to obtain an induced Q and U image. Since there are at least two sets of frames to compare (0, 60 and 120 degree polarizer angles), or possibly three (0, 45, 90 and 135 degrees), there are multiple Q and U images, that need to be added together with appropriate factors.

Finally, the final induced Q and U images are subtracted from the demodulated Q and U images, to remove the induced polarization.

Note: the correction for 4 different angles has not been implemented as of this writing.


## Wrapper

The modules `main` function wraps this procedure so that it accepts a set of files of polarization images, which can be multi-extension FITS files, and outputs a set of (multi-extension) FITS files of Stokes images.

In principle, frames are processed independently, and can thus be processed in parallel. This is not implemented, for two reasons:

- the final pipeline where this code may run may already take care of this

- implementing this directly in Python, for example in this `main` function, poses the additional issue that it requires use of shared memory. While perhaps not that difficult to implement, doing so may distract from the actual implementation.



## Configuration


All settings are given through a configuration file that is provided as the only argument to the module. For example,

```
python -m polaris.pipeline polaris.cfg
```

The configuration file is in TOML format; its sections are detailed inside the file itself, as comments.


[^1]: Van Vorstenbosch, A.T.P., 2019, "Theoretical performance of a wide field polarimetry survey using a polarimetric upgrade of OmegaPOL on the VLT Survey Telescope", master thesis, University of Leiden.
