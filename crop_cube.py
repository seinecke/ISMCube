import click
import os

from astropy import units as u
from astropy.coordinates import Angle

from spectral_cube import SpectralCube

import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)

@click.command()
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False))
@click.argument('glon_lo', type=Angle)
@click.argument('glon_hi', type=Angle)
@click.option('--glat_lo', type=Angle, default=None)
@click.option('--glat_hi', type=Angle, default=None)
@click.option('--output_dir', '-o', default='./', type=click.Path(exists=True, 
                dir_okay=True))
@click.option('--overwrite', is_flag=True, default=False)
def main(input_file, glon_lo, glon_hi, glat_lo, glat_hi, output_dir, overwrite):
    """
    Simple script to crop a Mopra cube in Galactic longitude.
    Optional: Crop also in Galactic latitude.

    Example: python crop_cube.py G340-350-12CO.fits 342d 344d
    """
    cube = SpectralCube.read(input_file).with_spectral_unit(u.km/u.s)


    if glat_lo is not None and glat_hi is not None:
        if glat_lo > glat_hi:
            raise Exception('glat_lo needs to be larger than glat_hi')

    if glat_lo is None:
        glat_lo = 'min'
    elif glat_lo < cube.latitude_extrema[0]:
        glat_lo = 'min'
        warnings.warn('Requested minimum latitude is smaller than extension of cube. \
            Minimum latitude of cube will be used.', stacklevel=2)
        
    if glat_hi is None:
        glat_hi = 'max'
    elif glat_hi > cube.latitude_extrema[1]:
        glat_hi = 'max'
        warnings.warn('Requested maximum latitude is larger than extension of cube. \
            Maximum latitude of cube will be used.', stacklevel=2)
    

    if glon_lo > glon_hi:
        raise Exception('glon_lo needs to be larger than glon_hi')

    if glon_lo < cube.longitude_extrema[0]:
        glon_lo = 'min'
        warnings.warn('Requested minimum longitude is smaller than extension of cube. Minimum longitude of cube will be used.', stacklevel=2)
    
    if glon_hi > cube.longitude_extrema[1]:
        glon_hi = 'max'
        warnings.warn('Requested maximum longitude is larger than extension of cube. Maximum longitude of cube will be used.', stacklevel=2)


    cube = cube.subcube(
        glon_lo, glon_hi,
        glat_lo, glat_hi,
        )
    
    print('\n', cube)

    basename = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f'{output_dir}{basename}_cropped.fits'
    cube.write(output_file, format='fits', overwrite=overwrite)
    print(f'Cropped cube written to {output_file}')


if __name__ == '__main__':
    main()