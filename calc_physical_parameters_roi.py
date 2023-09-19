import click
import os

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from photutils.aperture import SkyCircularAperture, CircularAperture
import numpy as np
from scipy.optimize import fsolve

import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)

@click.command()
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False))
@click.option('--vel', type=u.Quantity, default=None)
@click.option('--glon', type=Angle, default=None)
@click.option('--glat', type=Angle, default=None)
@click.option('--radius', type=Angle, default='0.5d')
@click.option('--xfactor', type=u.Quantity, default=1.8e20)
@click.option('--output_dir', '-o', default='./', type=click.Path(exists=True, 
                dir_okay=True))
@click.option('--overwrite', is_flag=True, default=False)
def main(input_file, vel, glon, glat, radius, xfactor, output_dir, overwrite):
    """
    Simple script to calculate basic physical parameters, 
    such as average brightness temperature, column density, 
    number density, near and far distance.

    Example: python calc_physical_parameters_roi.py G340-350-12CO_integrated.fits --vel -10
    """
    if vel is None:
        raise Exception('Please specify the velocity in km/s')

    xfactor *= u.s / u.cm**2 / u.K / u.km
    vel *= u.km/u.s


    f = fits.open(input_file)
    header = f[0].header
    data = u.Quantity(f[0].data, unit=header['BUNIT'])
    wcs = WCS(header)

    if glon is None or glat is None:
        R = vlsr_to_R(vel, l=wcs.wcs.crval[0]*wcs.wcs.cunit[0], b=wcs.wcs.crval[1]*wcs.wcs.cunit[1])
        dist_near, dist_far = R_to_dist(R, l=wcs.wcs.crval[0]*wcs.wcs.cunit[0], b=wcs.wcs.crval[1]*wcs.wcs.cunit[1])
        print(f'Near distance: {dist_near}')
        print(f'Far distance: {dist_far}')

    else: 
        R = vlsr_to_R(vel, l=glon, b=glat)
        dist_near, dist_far = R_to_dist(R, l=glon, b=glat)
        print(f'Near distance: {dist_near}')
        print(f'Far distance: {dist_far}')

    if glon is None or glat is None:
        W = np.nanmean(data)

    else:
        position = SkyCoord(l=glon, b=glat, frame='galactic')
        aperture = SkyCircularAperture(position, radius).to_pixel(wcs)
        mask = aperture.to_mask(method='center') 
        W = np.nanmean(mask.multiply(data))

    print(f'Average integrated brightness temperature : {W}')


    N = xfactor * W
    print(f'Corresponding column density: {N.to(u.cm**-2)}')


    n = N / (2 * dist_near * np.tan(radius))
    print(f'Corresponding number density (near distance): {n.to(u.cm**-3)}')

    n = N / (2 * dist_far * np.tan(radius))
    print(f'Corresponding number density (far distance): {n.to(u.cm**-3)}')


# Galactic Rotation Curve Model
# Brand and Blitz 1993
R0 = 8.5 * u.kpc
theta0 = 220 * u.km / u.s

# Table 4
a1 = 1.00767
a2 = 0.0394
a3 = 0.00712

def vlsr_to_R(v, l, b):
    
    def func(R):
        x = R / R0.to(u.kpc).value
        f = theta0 * ( (a1 * x**a2 + a3) / x - 1) \
            * np.sin(l) * np.cos(b) - v
        return f
    
    return fsolve(func, 1)[0] * u.kpc


def R_to_dist(R, l, b):
    dist_near = R0 * np.cos(l) / np.cos(b) \
            - (1/np.cos(b)**2) * np.sqrt(R**2 - (R0 * np.sin(l))**2)
    dist_far = R0 * np.cos(l) / np.cos(b) \
            + (1/np.cos(b)**2) * np.sqrt(R**2 - (R0 * np.sin(l))**2)
    return dist_near, dist_far

if __name__ == '__main__':
    main()

