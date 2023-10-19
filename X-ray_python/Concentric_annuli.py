#/Users/thomas_1/Desktop/A1446/4975/repro/acisf04975_broad_thresh.img

import numpy as np
import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

if os.path.exists('annuli.reg'):
    os.remove('annuli.reg')

# Read the FITS file
hdulist = fits.open('/Users/thomas_1/Desktop/A1446/4975/repro/flux/broad_thresh.img')
data = hdulist[0].data
wcs = WCS(hdulist[0].header)

# Define the center of the annulus in sky coordinates
center = SkyCoord('12h02m04.1530s', '+58d02m00.806s', frame='fk5')

# Define the inner and outer radii of the annulus in arcsec
inner_radius = 0.5
outer_radius = 1

counts_threshold = 3500
num_annuli = 10

for i in range(1,num_annuli+1):

    photon_counts = 0

    while photon_counts < counts_threshold:

        # Convert the radii to pixel units
        center_pix = wcs.world_to_pixel(center)
        inner_radius_pix = (inner_radius * u.arcsec).to(u.deg) / wcs.pixel_scale_matrix[1, 1]
        outer_radius_pix = (outer_radius * u.arcsec).to(u.deg) / wcs.pixel_scale_matrix[1, 1]

        # Create a boolean mask for the annulus
        y, x = np.indices(data.shape)
        r = np.sqrt((x - center_pix[0])**2 + (y - center_pix[1])**2)
        annulus_mask = (r >= inner_radius_pix.value) & (r <= outer_radius_pix.value)

        # Extract the annulus from the data
        annulus_data = np.where(annulus_mask, data, 0)

        # Calculate the total number of photon counts in the annulus
        photon_counts = np.sum(annulus_data)
        outer_radius += 1

    # Print the result
    print(f"Number of photon counts in annulus {i}: {photon_counts}")

    # Write the annulus to a .reg file
    with open('annuli.reg', 'a') as f:
        f.write(f"annulus({center.to_string('hmsdms', precision=2)}, {inner_radius}\", {outer_radius}\", 0) # color=green\n")

    with open(f'annulus{i}_toconvert.reg', 'w') as f:
        f.write(f"annulus({center.to_string('hmsdms', precision=2)}, {inner_radius}\", {outer_radius}\", 0) # color=green\n")

    inner_radius = outer_radius