from astropy.io import fits
import numpy as np
from lib_fits import flatten


def process_fits(input_filename, output_filename):
    with fits.open(input_filename) as hdul:

        head, data = flatten(input_filename)

        fits.writeto(output_filename, data, header = head)
        print(f'File salvato come {output_filename}')

# Esempio di utilizzo
input_filename = 'A1213_HBA-MFS-image.fits'  # Sostituisci con il percorso del tuo file FITS di input
output_filename = 'output_file.fits'  # Sostituisci con il percorso del tuo file FITS di output

process_fits(input_filename, output_filename)
