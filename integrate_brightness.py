import click
import os

from astropy import units as u
import numpy as np

from spectral_cube import SpectralCube

import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)

@click.command()
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False))
@click.option('--vel_lo', type=u.Quantity)
@click.option('--vel_hi', type=u.Quantity)
@click.option('--mask_level', '-m', default=None, type=float)
@click.option('--output_dir', '-o', default='./', type=click.Path(exists=True, 
                dir_okay=True))
@click.option('--overwrite', is_flag=True, default=False)
def main(input_file, vel_lo, vel_hi, mask_level, output_dir, overwrite):
    """
    Simple script to integrate cube over specified velocity range [in km/s].
    Optional: Mask out noise voxels, specified via mask_level x sigma.

    Example: python integrate_brightness.py G340-350-12CO.fits --vel_lo 0 --vel_hi 10 -m 3
    """
    cube = SpectralCube.read(input_file).with_spectral_unit(u.km/u.s)
    cube = cube.spectral_slab(vel_lo*u.km/u.s, vel_hi*u.km/u.s)

    if mask_level is not None:
        data = cube.unmasked_data[:,:,:]
        sigma = np.percentile(np.abs(data[data<0]), 68)
        mask = (cube > mask_level * sigma) 
        cube = cube.with_fill_value(0.0)
        cube = cube.with_mask(mask)
        cube.write('temp.fits', format='fits', overwrite=True)
        cube = SpectralCube.read('temp.fits')
        os.remove('temp.fits')

    mmap = cube.moment(order=0)

    basename = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f'{output_dir}{basename}_integrated.fits'
    mmap.write(output_file, format='fits', overwrite=overwrite)
    print(f'Integrated cube written to {output_file}')


if __name__ == '__main__':
    main()