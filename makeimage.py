#!/usr/bin/env python
# coding: utf-8

import aplpy
from astropy.io import fits
from astropy import units as u
import os, sys
import tempfile
import numpy as np
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from lib_fits import flatten
import argparse
import warnings
import lib_log
from termcolor import colored
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

logger = lib_log.logger

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='Make pdf image of fits file. Add contours, regions etc. \n    makeimage.py <fits image> *args')
parser.add_argument('image', default=[], help='Fits image to plot.')
parser.add_argument('--type', default='radio', help='radio/X')
parser.add_argument('--cmap', default='gnuplot2', help='Colormap to use.')
parser.add_argument('-z', '--z', type=float, help='Source redshift. Defaults to A1550.')
parser.add_argument('--radec', dest='radec', nargs=2, type=float, help='RA/DEC where to center image in deg (if not given, get full extent.)')
parser.add_argument('-s', '--size', default=600, type=float, help='Size of the final image in arcsec. Default to 600 arcsec')
parser.add_argument('-n', '--noise', type=float, help='Manually input noise (Jy/beam).')
parser.add_argument('--show_beam', action='store_true', help='Show image beam.')
parser.add_argument('--sbar', help='Length of size bar in kpc.')
parser.add_argument('--show_scale', action='store_true', help='Show arcsec to kpc scale.')
parser.add_argument('--show_contours', action='store_true', help='Show image contours.')
parser.add_argument('--interval', nargs='+', help='Manually input vmin, vmax, [vmid for log].')
parser.add_argument('--stretch', default='linear', help='Stretch to apply to image. Default: linear.')
parser.add_argument('--smooth', default=7, help='Smoothing to apply to the X-ray image. Default: 7')
parser.add_argument('--region', nargs='+', help='DS9 region file to plot.')
parser.add_argument('-o', '--outfile', default='plot', help='Prefix of output image.')
parser.add_argument('--fix_figure', action='store_true', help='Fix fits file if too many axes.')
parser.add_argument('--contours_color', default='white', help='Contours color.')
parser.add_argument('--cmap_type', default='sequential', choices=['sequential', 'qualitative'], help='Type of colormap to use')
parser.add_argument('--hide_nan', action='store_true', help='If set, NaN pixels are rendered with the same color as the minimum of the colorscale (invisible background).')
parser.add_argument('--fill_nan', action='store_true', help='If set, NaN pixels are filled by interpolation from neighboring pixels (only for display, the input FITS is not modified).')
parser.add_argument('--fill_nan_kernel', default=3, type=float, help='Std dev (pixels) of the Gaussian kernel used for NaN interpolation. Default: 3.')
parser.add_argument('--show_colorbar', action='store_true', help='Show colorbar on the image.')
parser.add_argument('--colorbar_unit', default=None, help='Label for the colorbar. If not given, uses the BUNIT keyword from the FITS header.')
parser.add_argument('--colorbar_scale', default=1.0, type=float, help='Multiplicative factor applied to the image data before display (e.g. 1000 to convert Jy/beam → mJy/beam). Default: 1.0.')
parser.add_argument('--bvec', default=None,
                    help='FITS file with polarization angle map (EVPA). '
                         'Vectors are rotated by 90° to show the B-field orientation.')
parser.add_argument('--bvec_pi', default=None,
                    help='FITS file with polarized intensity map. '
                         'If given, vector length is proportional to PI; otherwise uniform.')
parser.add_argument('--bvec_step', default=10, type=int,
                    help='Sampling step in pixels between adjacent B-field vectors. Default: 10.')
parser.add_argument('--bvec_scale', default=None, type=float,
                    help='Length scale for PI-proportional vectors: pixels per flux unit. '
                         'If omitted and --bvec_pi is given, auto-scaled to the 90th percentile of PI.')
parser.add_argument('--bvec_threshold', default=None, type=float,
                    help='Minimum polarized intensity (same units as --bvec_pi) to draw a vector. '
                         'Requires --bvec_pi. Overridden by --bvec_nsigma.')
parser.add_argument('--bvec_nsigma', default=None, type=float,
                    help='Set the vector threshold automatically as N×σ, where σ is the noise '
                         'estimated from the PI map (if --bvec_pi is given) or from the main image. '
                         'Overrides --bvec_threshold. Example: --bvec_nsigma 3')
parser.add_argument('--bvec_color', default='white',
                    help='Color of the B-field bars. Default: white.')
parser.add_argument('--bvec_linewidth', default=1.0, type=float,
                    help='Line width of the B-field bars. Default: 1.0.')
parser.add_argument('--bvec_pa_deg', action='store_true',
                    help='PA map is in degrees. If not set, the unit is auto-detected '
                         'from the BUNIT header keyword or from the value range of the map.')
parser.add_argument('--bvec_avg', default=1, type=int,
                    help='Block-averaging width in pixels before computing vector direction '
                         '(circular mean over avg×avg blocks, like CARTA\'s "Averaging width"). '
                         'Default: 1 (no averaging). CARTA default is 8.')
parser.add_argument('--bvec_evpa', action='store_true',
                    help='Draw EVPA vectors (electric field direction) instead of B-field vectors. '
                         'By default vectors are rotated 90° to show the magnetic field direction.')


def fix_aplpy_fits(aplpy_obj, dropaxis=2):
    temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
    temp_wcs = temp_wcs.dropaxis(dropaxis)
    aplpy_obj._wcs = temp_wcs


def calc_noise(data, niter=100, eps=1e-6):
    """Estimate RMS noise via iterative sigma-clipping (ignores NaN)."""
    rms = oldrms = 1.0
    for _ in range(niter):
        rms = np.nanstd(data)
        if np.abs(oldrms - rms) / rms < eps:
            return rms
        oldrms = rms
    raise RuntimeError('Noise estimation failed to converge.')


args = parser.parse_args()

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

filename = args.image
regions = args.region
z = args.z
plottype = args.type
plotcmap = args.cmap
noise = args.noise
size = args.size
beam = args.show_beam
sbar = args.sbar
scale = args.show_scale
contours = args.show_contours
output = args.outfile
interval = args.interval
imstretch = args.stretch
imsmooth = args.smooth
fix_figure = args.fix_figure
ctr_color = args.contours_color

if plotcmap in plt.colormaps():
    base_cmap = plt.get_cmap(plotcmap)
else:
    from palettable import cubehelix
    from palettable.colorbrewer import qualitative, sequential

    if hasattr(cubehelix, plotcmap):
        obj = getattr(cubehelix, plotcmap)
        base_cmap = obj.mpl_colormap
    elif hasattr(sequential, plotcmap):
        obj = getattr(sequential, plotcmap)
        base_cmap = obj.mpl_colormap
    elif hasattr(qualitative, plotcmap):
        obj = getattr(qualitative, plotcmap)
        base_cmap = ListedColormap(obj.mpl_colors)
    else:
        logger.warning(f"Cmap '{plotcmap}' not found anywhere")
        base_cmap = plt.get_cmap('viridis')

if args.cmap_type == 'qualitative':
    colors = base_cmap(np.linspace(0, 1, 8))
    plotcmap = ListedColormap(colors)
else:
    plotcmap = base_cmap


# --- Pre-processing dei dati: fill_nan e/o colorbar_scale (solo per il plot) ---
tmp_filename = None
needs_tmpfile = args.fill_nan or (args.colorbar_scale != 1.0)

if needs_tmpfile:
    with fits.open(filename) as hdul:
        hdu_orig = hdul[0].copy()

    data_orig = hdu_orig.data.copy()
    data_2d   = data_orig.squeeze()

    if args.fill_nan:
        logger.info('Interpolating NaN pixels for display...')
        kernel   = Gaussian2DKernel(x_stddev=args.fill_nan_kernel)
        data_2d  = interpolate_replace_nans(data_2d.astype(float), kernel)

    if args.colorbar_scale != 1.0:
        logger.info(f'Applying colorbar scale factor: ×{args.colorbar_scale}')
        data_2d = data_2d * args.colorbar_scale

    hdu_orig.data = data_2d.reshape(data_orig.shape)

    tmp_file     = tempfile.NamedTemporaryFile(suffix='.fits', delete=False)
    tmp_filename = tmp_file.name
    tmp_file.close()
    hdu_orig.writeto(tmp_filename, overwrite=True)
    filename = tmp_filename
    logger.info(f'Pre-processed data written to temp file: {tmp_filename}')

f = aplpy.FITSFigure(filename)

# Determina l'unità della colorbar: --colorbar_unit ha priorità, poi BUNIT dall'header
if args.colorbar_unit:
    colorbar_unit = args.colorbar_unit
else:
    with fits.open(args.image) as hdul:
        colorbar_unit = hdul[0].header.get('BUNIT', '')

if args.radec is None:
    logger.info('No RA/DEC provided: using FITS center')
    with fits.open(args.image) as hdul:
        wcs = WCS(hdul[0].header)
        ny, nx = hdul[0].data.shape[-2:]
        ra_c, dec_c = wcs.wcs_pix2world(nx/2, ny/2, 0)
    center = [ra_c, dec_c]
else:
    center = args.radec

# Se è stato applicato un fattore di scala, scala anche i valori di --interval
scaled_interval = None
if interval and args.colorbar_scale != 1.0:
    scaled_interval = [str(float(v) * args.colorbar_scale) for v in interval]
elif interval:
    scaled_interval = interval

if plottype=='radio':
    if fix_figure:
        fix_aplpy_fits(f)
    f.recenter(center[0], center[1], radius=(size/3600.))
    if scaled_interval:
        if imstretch != 'log':
            f.show_colorscale(cmap=plotcmap, stretch=imstretch, vmin=float(scaled_interval[0]), vmax=float(scaled_interval[1]))
        else:
            f.show_colorscale(cmap=plotcmap, stretch=imstretch, vmin=float(scaled_interval[0]), vmax=float(scaled_interval[1]), vmid=float(scaled_interval[2]))
    else:
        f.show_colorscale(cmap=plotcmap, stretch=imstretch)
    # f.add_colorbar()
    # f.colorbar.set_width(0.15)
    # f.colorbar.set_axis_label_text('[Jy beam$^{-1}$]')
    # #f.colorbar.set_axis_label_text('Polarisation degree')
    # f.colorbar.set_font(size=18)
    # f.colorbar.set_axis_label_font(size=18)
    # #f.colorbar.set_location('top')

    if args.show_colorbar:
        f.add_colorbar()
        f.colorbar.set_width(0.15)
        f.colorbar.set_font(size=18)
        f.colorbar.set_axis_label_font(size=18)
        if colorbar_unit:
            f.colorbar.set_axis_label_text(colorbar_unit)

    if contours:
        if not noise:

            head, datafits = flatten(filename)
            w = WCS(head)

            # CONVERT COORDINATES INTO PHYSICAL
            rapix, decpix = w.wcs_world2pix(center[0], center[1], 0)

            sizemask=400 #half size of the mask

            min_x = int(rapix) - sizemask
            max_x = int(rapix) - int(sizemask/2)
            min_y = int(decpix) - sizemask
            max_y = int(decpix) - int(sizemask/2)

            restrict = datafits[min_y:max_y, min_x:max_x]

            rms = calc_noise(restrict)

        else:
            rms = noise

        levels = np.array([3,6,12,24,48,96,192,384,768,1536,3072])
        levels = levels[:]*rms
        logger.info('Plotting contours..')
        f.show_contour(levels=levels, colors=f'{ctr_color}', linewidths=0.6)

    if beam:
        f.add_beam()
        f.beam.set_color('black')
        f.beam.set_frame(True)

elif plottype=='X':
    f.recenter(center[0], center[1], radius=(size/3600.))
    if scaled_interval:
        f.show_colorscale(cmap=plotcmap, stretch=imstretch, vmin=float(scaled_interval[0]), vmax=float(scaled_interval[1]), smooth=float(imsmooth))
    else:
        f.show_colorscale(cmap=plotcmap, stretch=imstretch, smooth=float(imsmooth))
    f.add_colorbar()
    f.colorbar.set_width(0.15)
    f.colorbar.set_font(size=12)
    f.colorbar.set_axis_label_font(size=12)
    unit_x = colorbar_unit if colorbar_unit else '[erg s$^{-1}$ cm$^{-2}$]'
    f.colorbar.set_axis_label_text(unit_x)

else:
    print(colored('ERROR:', 'red'), 'Plottype unrecognized.')
    sys.exit()

if args.hide_nan:
    current_cmap = f.image.get_cmap()
    cmap_fixed = plt.get_cmap(current_cmap.name).copy() if hasattr(current_cmap, 'name') else current_cmap.__copy__()
    cmap_fixed.set_bad(color=cmap_fixed(0.0))
    f.image.set_cmap(cmap_fixed)

# --- Vettori di campo magnetico da mappa di angolo di polarizzazione ---
if args.bvec:
    logger.info('Plotting B-field vectors from PA map...')

    # ── Lettura mappa PA ────────────────────────────────────────────────────
    with fits.open(args.bvec) as hdul:
        pa_data  = hdul[0].data.squeeze().astype(float)
        pa_bunit = hdul[0].header.get('BUNIT', '').strip().lower()

    # ── Auto-rilevamento unità (deg vs rad) ─────────────────────────────────
    # Priorità: 1) --bvec_pa_deg flag  2) BUNIT header  3) range dei valori
    finite_vals = pa_data[np.isfinite(pa_data)]
    if args.bvec_pa_deg:
        pa_is_deg = True
    elif 'deg' in pa_bunit:
        pa_is_deg = True
        logger.info('PA map: degrees detected from BUNIT header.')
    elif 'rad' in pa_bunit:
        pa_is_deg = False
        logger.info('PA map: radians detected from BUNIT header.')
    elif len(finite_vals) > 0 and np.max(np.abs(finite_vals)) > 2 * np.pi:
        pa_is_deg = True
        logger.info('PA map: degrees inferred from value range '
                    f'(|max|={np.max(np.abs(finite_vals)):.2f}).')
    else:
        pa_is_deg = False
        logger.info('PA map: assuming radians (values within ±2π).')

    pa_rad = np.radians(pa_data) if pa_is_deg else pa_data.copy()

    # ── B-field = EVPA + 90°  (oppure solo EVPA con --bvec_evpa) ────────────
    vec_rad = pa_rad if args.bvec_evpa else pa_rad + np.pi / 2.0

    # ── Lettura mappa PI (opzionale) ────────────────────────────────────────
    pi_data = None
    if args.bvec_pi:
        with fits.open(args.bvec_pi) as hdul:
            pi_data = hdul[0].data.squeeze().astype(float)
        if args.colorbar_scale != 1.0:
            pi_data = pi_data * args.colorbar_scale

    # ── Soglia automatica N×σ ────────────────────────────────────────────────
    bvec_threshold = args.bvec_threshold   # valore assoluto (può essere None)
    if args.bvec_nsigma is not None:
        # Stima il noise dalla mappa PI (se disponibile) o dall'immagine principale
        if pi_data is not None:
            noise_map = pi_data
            noise_src = '--bvec_pi map'
        else:
            with fits.open(args.image) as hdul:
                noise_map = hdul[0].data.squeeze().astype(float)
                if args.colorbar_scale != 1.0:
                    noise_map = noise_map * args.colorbar_scale
            noise_src = 'main image'
        sigma_est     = calc_noise(noise_map)
        bvec_threshold = args.bvec_nsigma * sigma_est
        logger.info(f'bvec threshold: {args.bvec_nsigma}×σ = '
                    f'{args.bvec_nsigma} × {sigma_est:.4e} = {bvec_threshold:.4e} '
                    f'(noise from {noise_src})')

    # ── Parametri di lunghezza ───────────────────────────────────────────────
    step             = args.bvec_step
    avg              = max(1, args.bvec_avg)
    half_len_uniform = step * 0.45

    if pi_data is not None:
        if args.bvec_scale is not None:
            pi_to_pix = args.bvec_scale
        else:
            valid_pi  = pi_data[np.isfinite(pi_data) & (pi_data > 0)]
            pi_ref    = np.percentile(valid_pi, 90) if len(valid_pi) > 0 else 1.0
            pi_to_pix = half_len_uniform / pi_ref

    # ── Griglia di campionamento ─────────────────────────────────────────────
    ny, nx = pa_data.shape
    ax     = f.ax
    pix_tr = ax.get_transform('pixel')

    xs = np.arange(step // 2, nx, step)
    ys = np.arange(step // 2, ny, step)

    segments = []
    n_drawn  = 0

    for iy in ys:
        for ix in xs:
            # ── Block averaging (media circolare) ────────────────────────────
            if avg > 1:
                r0 = max(0, iy - avg // 2); r1 = min(ny, r0 + avg)
                c0 = max(0, ix - avg // 2); c1 = min(nx, c0 + avg)
                block = vec_rad[r0:r1, c0:c1].ravel()
                block = block[np.isfinite(block)]
                if len(block) == 0:
                    continue
                # Circular mean: evita l'artefatto dello zero-crossing a ±π
                angle = np.arctan2(np.mean(np.sin(block)), np.mean(np.cos(block)))

                if pi_data is not None:
                    pi_block = pi_data[r0:r1, c0:c1].ravel()
                    pi_block = pi_block[np.isfinite(pi_block)]
                    pi_val   = float(np.mean(pi_block)) if len(pi_block) > 0 else np.nan
                else:
                    pi_val = None
            else:
                if iy >= ny or ix >= nx:
                    continue
                angle  = vec_rad[iy, ix]
                pi_val = float(pi_data[iy, ix]) if pi_data is not None else None

            if not np.isfinite(angle):
                continue

            # ── Soglia su PI ─────────────────────────────────────────────────
            if pi_data is not None and bvec_threshold is not None:
                if pi_val is None or not np.isfinite(pi_val) or pi_val < bvec_threshold:
                    continue

            # ── Direzione in pixel ────────────────────────────────────────────
            # Angolo α (N→E) in pixel (Nord=+y, Est=−x per immagini standard):
            #   dx = −sin(α),  dy = cos(α)
            dx = -np.sin(angle)
            dy =  np.cos(angle)

            # ── Lunghezza ────────────────────────────────────────────────────
            if pi_data is not None and pi_val is not None:
                if not np.isfinite(pi_val) or pi_val <= 0:
                    continue
                hl = min(pi_val * pi_to_pix, step * 0.45)
            else:
                hl = half_len_uniform

            segments.append([[ix - dx * hl, iy - dy * hl],
                              [ix + dx * hl, iy + dy * hl]])
            n_drawn += 1

    if segments:
        from matplotlib.collections import LineCollection
        lc = LineCollection(segments,
                            colors=args.bvec_color,
                            linewidths=args.bvec_linewidth,
                            capstyle='round',
                            transform=pix_tr)
        ax.add_collection(lc)

    vtype = 'EVPA' if args.bvec_evpa else 'B-field'
    logger.info(f'{vtype} vectors plotted: {n_drawn} bars '
                f'(step={step}px, avg={avg}px, units={"deg" if pa_is_deg else "rad"}).')

if scale:
    if not z:
        print(colored('ERROR:', 'red'), 'Redshift needed to plot scale.')
        sys.exit()
    arcsec_to_kpc = cosmo.kpc_proper_per_arcmin(z).value / 60
    text = f"{arcsec_to_kpc:.2f}"
    f.add_label(0.5, 0.97, "1'' / " + text + " kpc", relative=True, size=20, bbox=dict(facecolor=f'{ctr_color}', edgecolor='black', pad=5.0))

if sbar:
    if not z:
        print(colored('ERROR:', 'red'), 'Redshift needed to plot scale.')
        sys.exit()
    bar_length = sbar
    arcsec_to_kpc = cosmo.kpc_proper_per_arcmin(z).value/60
    if float(bar_length) >= 1000:
        f.add_scalebar((float(bar_length) / arcsec_to_kpc) / 3600, str(float(bar_length)/1000) + ' Mpc', color=f'{ctr_color}')
    else:
        f.add_scalebar((float(bar_length)/arcsec_to_kpc)/3600, bar_length+' kpc', color=f'{ctr_color}')
    f.scalebar.set_font(size='xx-large')


if regions is not None:
    if isinstance(regions, str):
        regions = [regions]
    for region in regions:
        logger.info('Plotting regions..')
        f.show_regions(region)

f.ticks.set_color('white')
f.ticks.set_linewidth(2)
f.tick_labels.set_font(size='xx-large')
f.axis_labels.set_font(size='xx-large')
f.axis_labels.set_ytext('Declination',)
f.axis_labels.set_xtext('Right Ascension')

if output:
    f.save(output+'.pdf', dpi=300)
else:
    f.save('plot.pdf', dpi=300)

# Rimuovi il file temporaneo creato da --fill_nan
if tmp_filename is not None and os.path.exists(tmp_filename):
    os.remove(tmp_filename)
    logger.info(f'Removed temp file: {tmp_filename}')

