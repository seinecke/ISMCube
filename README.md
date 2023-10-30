# ISMCube
Collection of scripts to work with ISM cubes (in FITS format), such as from the Mopra CO survey or the HI SGPS.

## Installation
The following packages are required: `matplotlib, numpy, astropy, scipy, spectral-cube, photutils, click`.

If using *anaconda*, you can create a new environment with the provided `environment.yml`. Download the file, go to the folder with the file, and execute:

```
conda env create -f environment.yml
```
Every time you want to use this environment, activate it via:
```
conda activate ism
```


## Usage
Every script comes with a help functionality, simply type `--help` after the function call. For example:
```
python crop_cube.py --help
```
Output:
```
Usage: crop_cube.py [OPTIONS] INPUT_FILE GLON_LO GLON_HI

  Simple script to crop a Mopra cube in Galactic longitude. Optional: Crop
  also in Galactic latitude.

  Example: python crop_cube.py G340-350-12CO.fits 342d 344d

Options:
  --glat_lo ANGLE
  --glat_hi ANGLE
  -o, --output_dir PATH
  --overwrite
  --help                 Show this message and exit.
```

## Typical Workflow
1. Crop cube (here from Galactic longitudes 342-344 deg):
```
python crop_cube.py G340-350-12CO.fits 342d 344d
```
2. Integrate brightness (here from velocities 0-10 km/s, masking out emission below 3 x the noise and using the cube cropped in the previous step):
```
python integrate_brightness.py G340-350-12CO_cropped.fits --vel_lo 0 --vel_hi 10 -m 3
```
3. Calculate physical parameters, including average brightness temperature, column density, number density, near and far distance (here the distance is calculated for a velocity of -10 km/s and using the integrated brightness from the previous step):
```
python calc_physical_parameters_roi.py G340-350-12CO_cropped_integrated.fits --vel -10
```

## Additional Information
### crop_cube.py
The cube needs to come in Galactic coordinates. For cropping, the longitudes (`glon_lo, glon_hi`) need to be provided, but the latitudes (`glat_lo, glat_hi`) are optional. The coordinates need to be followed by a `d` to indicate degree.

### integrate_brightness.py
By default the integration is performed over the entire velocity range in the cube. The minimum and maximum velocities (`vel_lo, vel_hi`) for the integration can be specified (in units of km/s). There is an option to remove noise before integrating. The noise is determined via the 68% percentile of negative temperatures and is removed by setting brightness temperatures to zero if they are below a multiple `mask_level` (also `m`) of this noise value. 

### calc_physical_parameters.py
The input here needs to be a map of integrated brightness. 

By default, the average brightness temperature, column density and number density are calculated over the entire map. There is an option to calculate these parameters for a smaller, circular region, defined by the circle's centre (`glon`, `glat`) and the circle's radius `radius`. By default, an X-factor (`xfactor`) of 1.8e20 s/cm2/K/km is used to calculate the column density - this factor should be adapted especially when not using 12CO. 

For calculating the distance, the Galactic Rotation Model by Brand & Blitz (1993) is used. The velocity `vel` (which should be within the velocity range that was integrated over) needs to be specified in km/s and is used in the calculation of the near and far distance.
