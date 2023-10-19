#!/usr/local/miniconda3/envs/science/bin/python
# -*- coding: utf-8 -*-

import os
import pyregion
from astropy.io import fits
#import pyfits
import logging
import numpy

def add_coloring_to_emit_ansi(fn):

	def new(*args):
		levelno = args[0].levelno
		if(levelno>=50):
			color = '\x1b[31m' # red
			pass
		elif(levelno>=40):
			color = '\x1b[31m' # red
			pass
		elif(levelno>=30):
			color = '\x1b[33m' # yellow
			pass
		elif(levelno>=20):
			color = '\x1b[32m' # green
			pass
		elif(levelno>=10):
			color = '\x1b[35m' # pink
			pass
		else:
			color = '\x1b[0m' # normal
			pass
		args[0].msg = color + args[0].msg +  '\x1b[0m'  # normal
		return fn(*args)
		pass
	return new
	pass

###############################################################################
def main(fits_name = None, z = 0.0, region_name = None, output_name = None):

    ## check fits image if present
    if not fits_name:
        logging.error('No FITS image specified.')
        return(1)
    if not os.path.isfile(fits_name):
        logging.error(str(fits_name) + ' could not be found.')
        return(1)

    ## check regions file if present
    if not region_name:
        logging.error('No region file specified.')
        return(1)
    if not os.path.isfile(region_name):
        logging.error(str(region_name) + ' could not be found.')
        return(1)

    ## reading image file
    logging.info('Reading FITS image: ' + str(fits_name))
    image  = fits.open(fits_name)[0]
    data   = fits.open(fits_name)['PRIMARY'].data
    header = image.header

    ## reading region file
    logging.info('Reading region file: ' + str(region_name))
    regions = pyregion.open(region_name).as_imagecoord(header)

    ## writing solutions to file
    if output_name:
        logging.info('Output file will be written to: ' + str(output_name))
        output = open(output_name, 'w')

    ## iterate over region files and calculate the sum, mean, standard deviation
    logging.info('region sum      pixel     mean      stdev')
    region_length = '{:>' + str(len(str(len(regions)))) + '}'
    for i, item in enumerate(regions):
        region = pyregion.ShapeList([item])
        mask = region.get_mask(hdu = image)
        image_field = numpy.copy(data)
        
        image_field[~mask] = 0
        
        filtered_area = image_field[image_field != 0]
        pixel_sum  = '{:9.2e}'.format(numpy.sum(filtered_area))
        pixel_mean = '{:9.2e}'.format(numpy.mean(filtered_area))
        pixel_std  = '{:9.2e}'.format(numpy.std(filtered_area))
        npix       = '{:9.2e}'.format(len(filtered_area))
        index      = region_length.format(i)
        
        values_to_print = index + ' ' + pixel_sum + ' ' + npix + ' ' + pixel_mean + ' ' + pixel_std
        logging.info(values_to_print)
        
        if output_name:
            output.write(values_to_print + '\n')

    if output_name:
        output.close()

    return(0)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Extracts pixel statistics from a FITS image.')

    parser.add_argument('fits',       type = str,   default = None, help = 'FITS image to extract.')
    parser.add_argument('--region',   type = str,   default = None, help = 'Region file to be used')
    parser.add_argument('--output',   type = str,   default = None, help = 'Output file to write results to')

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s","%Y-%m-%d %H:%M:%S")
    format_file   = logging.Formatter("%(asctime)s %(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log      = logging.StreamHandler()
    log.setFormatter(format_stream)
    log.emit = add_coloring_to_emit_ansi(log.emit)
    logging.root.addHandler(log)

    main(args.fits, region_name = args.region, output_name = args.output)
