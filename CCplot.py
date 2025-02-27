from astropy.io import fits
from astropy.wcs import WCS
from lib_fits import flatten
import sys, os, shutil, re
import pyregion
import numpy as np
import lib_img
from astropy.wcs import wcs
import warnings
from astropy.utils.exceptions import AstropyWarning
from termcolor import colored
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

warnings.filterwarnings("ignore", category=AstropyWarning)

if not os.path.exists('boxregions'):
    os.mkdir('boxregions')
else:
    shutil.rmtree('boxregions')
    os.mkdir('boxregions')

if not os.path.exists('excluded'):
    os.mkdir('excluded')
else:
    shutil.rmtree('excluded')
    os.mkdir('excluded')

def approx(number, precision=4):
    return round(np.ceil(number * 10 ** precision) / 10 ** precision, precision)

def get_data(fname,colname):
    data=fits.open(fname)
    data=data[1].data
    return data[colname]

def create_region(ra, dec, width, height, angle, color='yellow', name='region', dist_param=None):

    if dist_param is not None:
        regtext = [f'# Region file format: DS9 version 4.1 distance={dist_param}']
    else:
        regtext = [f'# Region file format: DS9 version 4.1']
    regtext.append(
        f'global color={color} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')
    regtext.append('fk5')
    regtext.append('box(' + str(ra) + ',' + str(dec) + f',{width}",{height}",{angle})')
    nline = '\n'
    target = f"{nline}{nline.join(regtext)}"
    if os.path.exists(f'{name}.reg'):
        os.remove(f'{name}.reg')
    with open(f'{name}.reg', 'w') as f:
        f.write(target)

    return target

def pixel_to_deg(fits_file, x, y):
    head, data = flatten(fits_file)
    wcs = WCS(head)
    ra, dec = wcs.all_pix2world(x, y, 1)
    return ra, dec

def create_grid_regions(fits_file, box_region_file, output_prefix, regsize, rainj, decinj, redshift):

    distance_list = []
    regions = pyregion.open(box_region_file)
    initial_box = regions[0].coord_list
    ra, dec, width, height, angle = float(initial_box[0]),float(initial_box[1]),float(initial_box[2]),float(initial_box[3]),float(initial_box[4])
    width = approx(width, 4)
    height = approx(height, 4)

    head, data = flatten(fits_file)
    beam_maj = head['BMAJ']

    box_size = beam_maj*regsize
    box_size = round(box_size,5)

    for i in range(int(width/box_size)):
        for j in range(int(height/box_size)):
            x1 = ra + ((width / 2)+box_size) - (i * box_size*1.2) #There is some approximation issue here probably, trying to compensate...
            y1 = dec - ((height / 2)-box_size/2) + (j * box_size)
            x2 = box_size*3600
            y2 = box_size*3600

            coord_1 = SkyCoord(rainj * u.deg, decinj * u.deg, frame='fk5')
            coord_2 = SkyCoord(x1 * u.deg, y1 * u.deg)
            sep = coord_1.separation(coord_2)*60 *u.arcmin/u.deg
            conv = cosmo.kpc_proper_per_arcmin(redshift)
            dist_kpc = sep*conv
            distance_list.append(dist_kpc.value)

            if dist_kpc.value < 100:
                create_region(x1, y1, x2, y2, angle, name=f'boxregions/{output_prefix}_{i}_{j}', color='red', dist_param=dist_kpc.value)
            elif dist_kpc.value >= 100 and dist_kpc.value < 200:
                create_region(x1, y1, x2, y2, angle, name=f'boxregions/{output_prefix}_{i}_{j}', color='yellow', dist_param=dist_kpc.value)
            elif dist_kpc.value >= 200 and dist_kpc.value < 300:
                create_region(x1, y1, x2, y2, angle, name=f'boxregions/{output_prefix}_{i}_{j}', color='green', dist_param=dist_kpc.value)
            elif dist_kpc.value >= 300 and dist_kpc.value < 400:
                create_region(x1, y1, x2, y2, angle, name=f'boxregions/{output_prefix}_{i}_{j}', color='cyan', dist_param=dist_kpc.value)
            elif dist_kpc.value >= 400 and dist_kpc.value < 500:
                create_region(x1, y1, x2, y2, angle, name=f'boxregions/{output_prefix}_{i}_{j}', color='blue', dist_param=dist_kpc.value)
            else:
                create_region(x1, y1, x2, y2, angle, name=f'boxregions/{output_prefix}_{i}_{j}', color='purple', dist_param=dist_kpc.value)


    return distance_list

def calc_flux(fits_file, region, threshold, main=True):

    roi = region

    with fits.open(fits_file) as phdu:
        head, lhdu = flatten(fits_file)
        gfactor = 2.0 * np.sqrt(2.0 * np.log(2.0))
        f = phdu[0]
        prhd = phdu[0].header
        units = prhd.get('BUNIT')
        if units is None:
            units = prhd.get('UNIT')
        if units != 'JY/BEAM' and units != 'Jy/beam':
            print('Warning: units are', units, 'but code expects JY/BEAM')
        bmaj = prhd.get('BMAJ')
        bmin = prhd.get('BMIN')

        bmaj = np.abs(bmaj)
        bmin = np.abs(bmin)

        w = wcs.WCS(prhd)
        cd1 = -w.wcs.cdelt[0]
        cd2 = w.wcs.cdelt[1]
        if ((cd1 - cd2) / cd1) > 1.0001 and ((bmaj - bmin) / bmin) > 1.0001:
            print('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

        bmaj /= cd1
        bmin /= cd2
        area = 2.0 * np.pi * (bmaj * bmin) / (gfactor * gfactor)

        d = [lhdu]

        region = pyregion.open(roi).as_imagecoord(prhd)

        for i, n in enumerate(d):
            mask = region.get_mask(hdu=f, shape=np.shape(n))
            data = np.extract(mask, d)
            nndata = data[~np.isnan(data)]
            flux = np.sum(nndata) / area

    if main:
        if flux < threshold:
            halfstring = roi.split("/")
            halfstring = halfstring[-1]
            os.rename(f'{roi}', f'excluded/{halfstring}')
            keep = False
            dist = None
        else:
            dist = read_distance(roi)
            keep = True

        return flux, keep, dist

    else:
        dist = read_distance(roi)
        keep = True

        return flux, keep, dist


def read_distance(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('# Region file format'):
            match = re.search(r'distance=(\d+)', line)
            if match:
                return int(match.group(1))
    return None


parser = argparse.ArgumentParser(description='Make color-color plot of 3 radio images. Creates the regions, extracts the flux density and plots the spectral index. Usage: CCplot_radio.py -mainimage *.fits -otherimages *.fits *.fits -region *.reg.')
parser.add_argument('-mainimage', default=[], help='Fits image on which the region grid is designed.')
parser.add_argument('-otherimages', nargs=2, help='List of images at other frequencies on which to compute the color-color.')
parser.add_argument('-region', help='DS9 region file in which the grid is designed.')
parser.add_argument('-radec', type=float, nargs=2, help='RA and DEC in degree of injection point. Will be used to calculate the distance of the regions.')
parser.add_argument('-redshift', type=float, help='Redshift of the source.')
parser.add_argument('--multiplier', type=int, default=5, help='Defines the N number by which the sigma is multiplied to exclude unwanted pixels.')
parser.add_argument('--regsize', type=float, default=2, help='Defines the extent of the regions. This number gets multiplied by the beam size to derive the final extent.')
parser.add_argument('--noise', type=float, help='Manually input the noise level of the main image.')
parser.add_argument('--plotgrid', action='store_true', help='Use makeimage.py to plot the grid above the main image.')
parser.add_argument('--plotgrid_excluded', action='store_true', help='Use makeimage.py to plot the excluded regions above the main image.')

args = parser.parse_args()

fits_file = args.mainimage
box_region_file = args.region
other_images = args.otherimages
regsize = args.regsize
multiplier = args.multiplier
noise = args.noise
ra = args.radec[0]
dec = args.radec[1]
redshift = args.redshift
plotgrid = args.plotgrid
plotgrid_excluded = args.plotgrid_excluded

output_prefix = 'region'

print(colored('Creating grid...', 'green'))
dist_list = create_grid_regions(fits_file, box_region_file, output_prefix, regsize, ra, dec, redshift)

nimage = lib_img.Image(fits_file)
if not args.noise:
    noise = nimage.getNoise()
threshold = multiplier*noise

fluxes_main = []
dist_final = []
print(colored(f'Calculating flux density for regions of {fits_file}...', 'green'))
for region in os.listdir('boxregions'):
    main, keep, distance = calc_flux(fits_file, f'boxregions/{region}', threshold, main=True)
    if keep == True:
        fluxes_main.append(main)
        dist_final.append(distance)

if plotgrid:
    print(colored(f'Plotting the grid overlay...', 'green'))
    try:
        os.system(f'python3 makeimage.py --z {redshift} --region boxregions/* --radec {ra} {dec} -o grid_plot --show_beam {fits_file}')
    except:
        print(colored(f'Grid plot gave an error, maybe check whether makeimage.py is in the working directory...', 'red'))

if plotgrid_excluded:
    print(colored(f'Plotting the excluded regions...', 'green'))
    try:
        os.system(f'python3 makeimage.py --z {redshift} --region excluded/* --radec {ra} {dec} -o excluded_plot --show_beam {fits_file}')
    except:
        print(colored(f'Grid plot gave an error, maybe check whether makeimage.py is in the working directory...', 'red'))

fluxes_lists = [[] for _ in range(len(other_images))]
dist_lists = [[] for _ in range(len(other_images))]
for n, image in enumerate(other_images):
    print(colored(f'Calculating flux density for image {image}...', 'green'))
    fluxes_list = fluxes_lists[n]
    dist_list = dist_lists[n]
    nimage = lib_img.Image(image)
    try:
        noise = nimage.getNoise()
    except:
        noise = input(colored(f'Noise calculation for {image} did not converge, please input manually (Jy/beam):', 'red'))
    threshold = multiplier * noise
    for region in os.listdir('boxregions'):
        intflux, keep, dist_final2 = calc_flux(image, f'boxregions/{region}', threshold, main=False)
        fluxes_list.append(intflux)
        dist_list.append(dist_final2)

all_images = [fits_file, other_images[0], other_images[1]]
all_lists = [fluxes_main, fluxes_lists[0], fluxes_lists[1]]
all_dists = dist_final

freqs = []
for image in all_images:
    head, data = flatten(image)
    freq = head['FREQ']
    freqs.append(freq)
sorted_freqs, sorted_all_images, sorted_flux = zip(*sorted(zip(freqs, all_images, all_lists)))

with open('fluxes.txt', "w") as f:
    f.write(f"#F_{round(sorted_freqs[0]/1e6)}MHz\tF_{round(sorted_freqs[1]/1e6)}MHz\tF_{round(sorted_freqs[2]/1e6)}MHz\terrF_{round(sorted_freqs[0]/1e6)}MHz"
            f"\terrF_{round(sorted_freqs[1]/1e6)}MHz\terrF_{round(sorted_freqs[2]/1e6)}MHz\talpha_low\talpha_high\terralpha_low\terralpha_high\n")
    f.write('\n')
    i=0
    for val_1, val_2, val_3 in zip(*sorted_flux):
        if not val_1 < 0 and not val_2 < 0 and not val_3 < 0:
            f.write(f"{val_1}\t{val_2}\t{val_3}\t")
            err_1 = val_1 * 0.1
            err_2 = val_2 * 0.1
            err_3 = val_3 * 0.06
            f.write(f"{err_1}\t{err_2}\t{err_3}\t")
            f.write(f'{dist_final[i]}\n')
            i += 1

datafile = np.loadtxt('fluxes.txt')

alphas_low = []
alphas_high = []
erralphas_low = []
erralphas_high = []

for i, val in enumerate(datafile):
    flux_low = datafile[i, 0]
    flux_mid = datafile[i, 1]
    flux_high = datafile[i, 2]
    err_low = datafile[i, 3]
    err_mid = datafile[i, 4]
    err_high= datafile[i, 5]
    dist = datafile[i, 6]

    alpha_low = (np.log10(flux_low / flux_mid)) / (np.log10((round(sorted_freqs[0]/1e6)) / (round(sorted_freqs[1]/1e6))))
    erralpha_low = (1 / (np.log(round(sorted_freqs[1]/1e6)) - np.log(round(sorted_freqs[0]/1e6)))) * (np.sqrt(((err_low / flux_low) ** 2) + ((err_mid / flux_mid) ** 2)))
    #erralpha = (1 / (np.log(freq2) - np.log(freq1))) * (np.sqrt(((err1 / flux1) ** 2) + ((err2 / flux2) ** 2)))

    if not np.isnan(alpha_low):
        alphas_low.append(alpha_low)
        erralphas_low.append(erralpha_low)

    alpha_high = (np.log10(flux_mid / flux_high)) / (np.log10((round(sorted_freqs[1] / 1e6)) / (round(sorted_freqs[2] / 1e6))))
    erralpha_high = (1 / (np.log(round(sorted_freqs[2]/1e6)) - np.log(round(sorted_freqs[1]/1e6)))) * (np.sqrt(((err_mid / flux_mid) ** 2) + ((err_high / flux_high) ** 2)))

    if not np.isnan(alpha_high):
        alphas_high.append(alpha_high)
        erralphas_high.append(erralpha_high)

if os.path.exists('Results.fits'):
    os.remove('Results.fits')

a0 = np.array(datafile[:,6])
a1 = np.array(datafile[:,0])
a2 = np.array(datafile[:,1])
a3 = np.array(datafile[:,2])
a4 = np.array(datafile[:,3])
a5 = np.array(datafile[:,4])
a6 = np.array(datafile[:,5])
a7 = np.array(alphas_low)
a8 = np.array(alphas_high)
a9 = np.array(erralphas_low)
a10 = np.array(erralphas_high)
col0 = fits.Column(name=f'Distance [kpc]', array=a0, format='F', unit = 'kpc')
col1 = fits.Column(name=f'F_{round(sorted_freqs[0]/1e6)}MHz', array=a1, format='F', unit = 'Jy')
col2 = fits.Column(name=f'F_{round(sorted_freqs[1]/1e6)}MHz', array=a2, format='F', unit = 'Jy')
col3 = fits.Column(name=f'F_{round(sorted_freqs[2]/1e6)}MHz', array=a3, format='F', unit = 'Jy')
col4 = fits.Column(name=f'errF_{round(sorted_freqs[0]/1e6)}MHz', array=a4, format='F', unit = 'Jy')
col5 = fits.Column(name=f'errF_{round(sorted_freqs[1]/1e6)}MHz', array=a5, format='F', unit = 'Jy')
col6 = fits.Column(name=f'errF_{round(sorted_freqs[2]/1e6)}MHz', array=a6, format='F', unit = 'Jy')
col7 = fits.Column(name='alpha_low', array=a7, format='F')
col8 = fits.Column(name='alpha_high', array=a8, format='F')
col9 = fits.Column(name='erralpha_low', array=a9, format='F')
col10 = fits.Column(name='erralpha_high', array=a10, format='F')
cols = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto('Results.fits')

fdata="Results.fits"

alpha_low = get_data(fdata,'alpha_low')
erralpha_low = get_data(fdata,'erralpha_low')
alpha_high = get_data(fdata,'alpha_high')
erralpha_high = get_data(fdata,'erralpha_high')
distance = get_data(fdata,'Distance [kpc]')

sc = plt.scatter(alpha_low, alpha_high, c=distance, cmap='rainbow_r', marker='s', alpha=1, s=80)
plt.errorbar(alpha_low, alpha_high, xerr=erralpha_low, yerr=erralpha_high, fmt='none', ecolor='black', capsize=2, alpha=1, elinewidth=1)
cbar = plt.colorbar(sc, label='Distance (kpc)', alpha=1)
plt.xlabel(r'$\alpha_{54 \rm MHz}^{144 \rm MHz}$', fontsize=20)
plt.ylabel(r'$\alpha_{144 \rm MHz}^{383 \rm MHz}$', fontsize=20)
plt.tick_params(labelsize=19)
# plt.xlim(-5,0)
# plt.ylim(-5,0)
bisector = np.arange(-5,0,0.1)
plt.errorbar(x=bisector, y=bisector,fmt='--k',capsize=3)
plt.gcf().set_size_inches(10, 10)
plt.savefig('CCP.pdf')

os.remove('fluxes.txt')














