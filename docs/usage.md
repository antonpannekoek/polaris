For simple usage, start with a simulation of a set of files. Edit the `sim.cfg` file (make a copy first if wanted), and adjust the settings to your liking: the size and number of CCDs, a possible shift, the number of angles etc.

Use `action = "create"` to create a database with selected (unpolarized) field stars; make sure the output file (`dbname`) can be overwritten or is a non-existing file.


Then, simulate some data; this will be written as FITS file. The file name is fixed at the moment, and is `polaris-deg<angle>.fits`, with `<angle>` one of the polarizer angles. To run a simulation:

```
python -m polaris.simulation sim.cfg
```


To then run the pipeline on the output data, edit the `polaris.cfg` file. Adjust the input files to match the files created in the simulation step. Set the correction step for induced polarization to your liking, and adjust the output file name accordingly. You can leave the debugging on for some intermediate files (e.g. the gradient of Stokes I as a FITS file), or turn it off.

Run the pipeline with

```
python -m polaris.pipeline polaris.cfg
```

and view the output FITS files. Note that there may be multiple extensions, depending on the input files; the primary HDU of the FITS file is always empty.

There will also a region file that contains the positions of the polarized stars, for ease of viewing in e.g. DS9, or for other analysis.
