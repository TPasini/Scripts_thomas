import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.colors as colors
from lib_fits import flatten
from astropy.wcs import WCS


# Funzione per salvare la maschera come file .mask.fits
def salva_maschera(maschera, nome_file):
    header = fits.Header()
    header['BITPIX'] = 8
    header['NAXIS'] = 2
    header['NAXIS1'] = maschera.shape[1]
    header['NAXIS2'] = maschera.shape[0]

    hdu = fits.PrimaryHDU(maschera.astype(np.uint8), header)
    hdu.writeto(nome_file, overwrite=True)


# Funzione per leggere un file .reg e ottenere le coordinate del cerchio
def leggi_file_reg(file_reg):
    with open(file_reg, 'r') as f:
        linee = f.readlines()

    x_centro, y_centro, raggio = None, None, None
    for linea in linee:
        if 'circle' in linea:
            valori = linea.split('(')[1].split(')')[0].split(',')
            x_centro = float(valori[0])
            y_centro = float(valori[1])
            raggio = float(valori[2])

    return x_centro, y_centro, raggio


# Funzione per creare una maschera basata sul valore di soglia
def crea_maschera(image, rms_noise, soglia, x_centro, y_centro, raggio):
    maschera = np.zeros(image.shape, dtype=bool)

    y_size, x_size = image.shape

    for y in range(y_size):
        for x in range(x_size):
            distanza = np.sqrt((x - x_centro) ** 2 + (y - y_centro) ** 2)
            if distanza <= raggio:
                maschera[image > soglia * rms_noise] = True

    return maschera


def calc_noise(data, niter=100, eps=1e-6):
    """
    Return the rms
    """

    rms = 1.;
    oldrms = 1.
    for i in range(niter):
        rms = np.nanstd(data)
        if np.abs(oldrms - rms) / rms < eps:
            return rms

        oldrms = rms

    raise Exception('Noise estimation failed to converge.')


# Leggi i file FITS
file_fits1 = 'images/PSZ2G099.86+58.45_briggs_HBA-MFS-image.fits'
file_fits2 = 'images/PSZ2G099.86+58.45_briggs_LBA-MFS-image.fits'
hdulist1 = fits.open(file_fits1)
hdulist2 = fits.open(file_fits2)

# Leggi i dati delle immagini
image_data1 = hdulist1[0].data[0, 0]
image_data2 = hdulist2[0].data[0, 0]

# Leggi la regione circolare dal file .reg
file_reg = 're/r_e_PSZ2G099.86+58.45.reg'
x_centro, y_centro, raggio = leggi_file_reg(file_reg)

# Calcola il raggio in pixel per ciascuna immagine
pixel_scale1 = hdulist1[0].header['CDELT2']  # scala dei pixel per immagine 1
pixel_scale2 = hdulist2[0].header['CDELT2']  # scala dei pixel per immagine 2
raggio_pixel1 = raggio / pixel_scale1
raggio_pixel2 = raggio / pixel_scale2

# Calcola l'RMS noise per ciascuna immagine
head, datafits = flatten(file_fits1)
w = WCS(head)

rapix, decpix = w.wcs_world2pix(x_centro, y_centro, 0)

sizemask=400 #half size of the mask

min_x = int(rapix) - sizemask
max_x = int(rapix) - int(sizemask/2)
min_y = int(decpix) - sizemask
max_y = int(decpix) - int(sizemask/2)

restrict1 = image_data1[min_y:max_y, min_x:max_x]
rms_noise1 = calc_noise(restrict1)

restrict2 = image_data2[min_y:max_y, min_x:max_x]
rms_noise2 = calc_noise(restrict2)

# Crea le maschere
maschera1 = crea_maschera(image_data1, rms_noise1, 1, x_centro, y_centro, 0.5)
maschera2 = crea_maschera(image_data2, rms_noise2, 1, x_centro, y_centro, 0.5)

salva_maschera(maschera1, 'maschera1.mask.fits')
salva_maschera(maschera2, 'maschera2.mask.fits')

# Crea una maschera comune
maschera_comune = np.logical_and(maschera1, maschera2)

# Salva la maschera comune in un file .reg
with open('maschera_comune.reg', 'w') as f:
    f.write(f'circle({x_centro},{y_centro},{raggio})\n')
    for i in range(maschera_comune.shape[0]):
        for j in range(maschera_comune.shape[1]):
            if maschera_comune[i, j].any():
                f.write(f'point({i},{j})\n')

# Visualizza le immagini e le maschere
fig, axs = plt.subplots(2, 2)

# Plot dell'immagine 1
axs[0, 0].imshow(image_data1, cmap='gray', norm=colors.LogNorm())
axs[0, 0].set_title('Immagine 1')
axs[0, 0].axis('off')

# Zoom verso la regione di interesse per l'immagine 1
x_start, x_end = x_centro - raggio, x_centro + raggio
y_start, y_end = y_centro - raggio, y_centro + raggio
axs[0, 0].set_xlim(x_start, x_end)
axs[0, 0].set_ylim(y_start, y_end)

# Plot della maschera 1
axs[0, 1].imshow(maschera1, cmap='gray')
axs[0, 1].set_title('Maschera 1')
axs[0, 1].axis('off')

# Plot dell'immagine 2
axs[1, 0].imshow(image_data2, cmap='gray', norm=colors.LogNorm())
axs[1, 0].set_title('Immagine 2')
axs[1, 0].axis('off')

# Zoom verso la regione di interesse per l'immagine 2
axs[1, 0].set_xlim(x_start, x_end)
axs[1, 0].set_ylim(y_start, y_end)

# Plot della maschera 2
axs[1, 1].imshow(maschera2, cmap='gray')
axs[1, 1].set_title('Maschera 2')
axs[1, 1].axis('off')

plt.tight_layout()
plt.show()
