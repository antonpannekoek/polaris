# VSTPOL pipeline

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
