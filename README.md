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
python ../crop_cube.py --help
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
2. Integrate brightness (here from velocities 0-10 km/s and masking out emission below 3 x the noise):
```
python integrate_brightness.py G340-350-12CO.fits --vel_lo 0 --vel_hi 10 -m 3
```
3. Calculate physical parameters, including average brightness temperature, column density, number density, near and far distance (here the distance is calculated for a velocity of -10 km/s):
```
python calc_physical_parameters_roi.py G340-350-12CO_integrated.fits --vel -10
```
