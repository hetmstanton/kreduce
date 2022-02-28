"""
Adapted from: https://github.com/charlottenosam/kmos_tools/blob/master/kmos_tools/sky_clean.py
"""

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

import datetime

import os
import shutil
import glob

import warnings
warnings.filterwarnings("ignore")

def make_sky_residual_spectrum(exposure, show=True):
    """Calculate all the sky residual corrections for one DIT ~ a la Trevor Mendel/Charlotte Mason
    
    Generate median 1D sky spectra for each detector,
    only from galaxy cubes (skip star IFUs), and save them to use for subtraction.

    Note this only works because if objects are generally small with respect to the detector size.
    """

    print(f'\n*** [kreduce-sky]: calculating sky residuals in {exposure.filename}\n')

    # need to know which ifus contain stars so that we skip them:
    exposure.get_star_ifus()

    # Create stacks of 'empty' cubes:
    exposure.filename_skyspec = exposure.filename.replace('.fits', '_SKYSPEC.fits')

    detector1, detector2, detector3 = [], [], []

    for ifu in range(1, 25):

        if ifu in exposure.star_ifus:
            print(f'*** [kreduce-sky]: skipping IFU {ifu} (contains ref star)')
            continue

        ext = exposure.hdulist[f'IFU.{ifu}.DATA']

        # make sure to skip empty ifus
        if len(ext.shape) > 0:

            ifu_comment = f'HIERARCH ESO OCS ARM {ifu} NAME'
            ifu_header  = ext.header
            ifu_cube = ext.data

            if 1 <= ifu <= 8:
                detector1.append(ifu_cube)
            elif 9 <= ifu <= 16:
                detector2.append(ifu_cube)
            else:
                detector3.append(ifu_cube)

    len_for_stack = len(detector1) + len(detector2) + len(detector3)

    detector1, detector2, detector3 = np.array(detector1), np.array(detector2), np.array(detector3)

    detector_all = np.concatenate((detector1, detector2, detector3), axis=0)

    skyspec_1D_all = np.nanmedian(detector_all, axis=(0, 2, 3))
    skyspec_1D_all[np.isnan(skyspec_1D_all)] = 0.0

    skyspec_1D = {}
    detectors = [detector1, detector2, detector3]
    for i in range(len(detectors)):
        detector_sky = np.nanmedian(detectors[i], axis=(0, 2, 3))
        detector_sky[np.isnan(detector_sky)] = 0
        skyspec_1D[i] = detector_sky

    if show:

        plt.figure(figsize=(10, 5))

        plt.plot(skyspec_1D_all, lw=1, alpha=0.8, 
            label=f'All detectors ({detector_all.shape[0]} IFUs)', zorder=10)

        for i in range(len(skyspec_1D)):
            if skyspec_1D[i] is not None:
                plt.plot(skyspec_1D[i], lw=1, alpha=0.8, 
                    label=f'Detector {i} ({detectors[i].shape[0]} IFUs)')

        ymin, ymax = np.nanpercentile(skyspec_1D_all, 1), np.nanpercentile(skyspec_1D_all, 99)
        
        if ymin > 0.: ymin = -1.e-18

        plt.ylim(ymin, ymax)
        plt.xlabel('Wavelength [pix]')
        plt.ylabel('Flux')
        plt.title('Sky Subtraction Residuals')
        plt.legend()
        plt.tight_layout()
        plot_filename = exposure.filename_skyspec.replace('.fits', '.pdf')
        plt.savefig(f'{plot_filename}')

     # Save spectra to fits file
    # Headers
    prihdr = exposure.hdr.copy()
    prihdr.add_comment('Sky Spectrum on each detector, and median sky spectrum')
    hdr1D = fits.Header()
    hdr1D['SIMPLE'] = 'T'
    hdr1D['BITPIX'] = -32
    hdr1D['NAXIS'] = 1
    hdr1D['NAXIS1'] = 2048
    hdr1D['PCOUNT'] = 0
    hdr1D['GCOUNT'] = 1
    hdr1D['CUNIT1'] = ifu_header['CUNIT3']
    hdr1D['CRPIX1'] = ifu_header['CRPIX3']
    hdr1D['CRVAL1'] = ifu_header['CRVAL3']
    hdr1D['CDELT1'] = ifu_header['CDELT3']
    hdr1D['BUNIT'] = 'cgs'

    hdr1D_1, hdr1D_2, hdr1D_3 = hdr1D.copy(), hdr1D.copy(), hdr1D.copy()
    hdr1D['EXTNAME']   = 'ALL'
    hdr1D_1['EXTNAME'] = 'DETECTOR1'
    hdr1D_2['EXTNAME'] = 'DETECTOR2'
    hdr1D_3['EXTNAME'] = 'DETECTOR3'

    # Extensions
    hdu     = fits.PrimaryHDU(header=prihdr)
    hdu_all = fits.ImageHDU(skyspec_1D_all, header=hdr1D)
    hdu_1   = fits.ImageHDU(skyspec_1D[0], header=hdr1D_1)
    hdu_2   = fits.ImageHDU(skyspec_1D[1], header=hdr1D_2)
    hdu_3   = fits.ImageHDU(skyspec_1D[2], header=hdr1D_3)

    # Create hdu list and write
    hdulist = fits.HDUList([hdu, hdu_all, hdu_1, hdu_2, hdu_3])
    hdulist.writeto(exposure.filename_skyspec, overwrite=True)
    print('*** [kreduce-sky]: Saved sky fits file to {:s}'.format(exposure.filename_skyspec))


def sky_subtract_residual_sky(exposure, show=True):

    print(f'\n*** [kreduce-sky]: subtracting sky residuals in {exposure.filename}\n')

    exposure.filename_skycorr = exposure.filename.replace('.fits', '_SKYCORR.fits')
    skyspec_all = fits.open(exposure.filename_skyspec)

    # need to know which ifus contain stars so that we skip them:
    exposure.get_star_ifus()

    if show:
        if not os.path.exists(f'{exposure.filedir}/sky_correct_plots'):
            os.mkdir(f'{exposure.filedir}/sky_correct_plots')

    for ifu in range(1, 25):

        if ifu in exposure.star_ifus:
            continue

        ext = exposure.hdulist[f'IFU.{ifu}.DATA']

        # make sure to skip empty ifus
        if len(ext.shape) > 0:

            # Finding detector
            if 1 <= ifu <= 8:
                detector = 1
            elif 9 <= ifu <= 16:
                detector = 2
            elif 17 <= ifu <= 24:
                detector = 3

            sky  = skyspec_all[detector+1].data

            # Estimate 1d error from std of flux
            ext_copy = ext.data.copy()
            error1d = np.nanstd(ext_copy, axis=(1, 2))/np.sqrt(ext_copy.shape[1] * ext_copy.shape[2])

            if show:
                
                figure_mosaic = """
                ACC
                BCC
                """
                fig, axes = plt.subplot_mosaic(mosaic=figure_mosaic, figsize=(11,5))

                collapsed_original = np.nanmedian(ext_copy, axis=(1,2))
                axes["A"].plot(collapsed_original, lw=1.,
                    alpha=0.8, label='Original full spectrum', color='k')
                axes["A"].plot(sky, lw=1.,
                    alpha=0.8, label='Sky', color='C0', ls=':')
                ymin_sky, ymax_sky = np.nanpercentile(collapsed_original, 1), np.nanpercentile(collapsed_original, 99)
        
                if ymin_sky > 0.: ymin_sky = -1.e-18

                axes["A"].set_ylim(ymin_sky, ymax_sky)
                axes["A"].set_ylabel('Flux')
                axes["A"].legend()

                nanmask = np.isfinite(collapsed_original)
                xmin_hist, xmax_hist = np.nanpercentile(collapsed_original[nanmask], 1), np.nanpercentile(collapsed_original[nanmask], 99)
                limmask = (collapsed_original[nanmask] >= xmin_hist) & (collapsed_original[nanmask] <= xmax_hist) 
                collapsed_to_plot = collapsed_original[nanmask][limmask]

                axes["C"].hist(collapsed_to_plot, histtype='step', color='k', lw=1.,
                    bins=100)
                axes["C"].set_xlabel('Flux')

            # make a copy of the original data cube and and 
            # empty array that will hold the sky-subtracted copy:
            data_orig = np.copy(ext.data)
            data_corr = np.empty_like(data_orig)

            for i in range(data_orig.shape[1]):
                for j in range(data_orig.shape[2]):

                    pixel_spec = data_orig[:,i,j]

                    if np.isnan(pixel_spec).all():
                        data_corr[:,i,j] = np.zeros_like(data_orig[:,i,j]) * np.nan
                    else:
                        sky_norm =  get_sky_normalization(spec=pixel_spec, 
                            spec_err=error1d, sky=sky)
                        data_corr[:,i,j] = data_orig[:,i,j] - sky_norm*sky

            if show:
                collapsed_skycorr = np.nanmedian(data_corr, axis=(1,2))
                axes["B"].plot(collapsed_skycorr, lw=1.,
                    alpha=0.8, label='Sky corrected spectrum', color='C3')

                axes["B"].set_ylabel('Flux')
                axes["B"].set_xlabel('Wavelength [pix]')

                axes["B"].legend()
                axes["B"].set_ylim(ymin_sky, ymax_sky)

                nanmask = np.isfinite(collapsed_skycorr)
                xmin_hist, xmax_hist = np.nanpercentile(collapsed_skycorr[nanmask], 1), np.nanpercentile(collapsed_skycorr[nanmask], 99)
                limmask = (collapsed_skycorr[nanmask] >= xmin_hist) & (collapsed_skycorr[nanmask] <= xmax_hist) 
                collapsed_to_plot = collapsed_skycorr[nanmask][limmask]
                axes["C"].hist(collapsed_to_plot, histtype='step', color='C3', lw=1.,
                    bins=100)

                axes["C"].set_xlabel('Flux')

                plt.tight_layout()
                
                plot_filename = exposure.filename_skycorr.replace('.fits', f'_IFU{ifu}.pdf')
                fig.savefig(f'{plot_filename}')

                # remove the plot if it already exists:
                root_pfile = f'SCI_' + plot_filename.split('SCI_')[1]
                if os.path.exists(f"{exposure.filedir}/sky_correct_plots/{root_pfile}"):
                    os.remove(f"{exposure.filedir}/sky_correct_plots/{root_pfile}")

                shutil.move(src=os.path.abspath(plot_filename), 
                    dst=os.path.abspath(f'{exposure.filedir}/sky_correct_plots'))

            ext.data = data_corr

        exposure.hdulist.writeto(exposure.filename_skycorr, overwrite=True)


def get_sky_normalization(spec, spec_err, sky):

    # mask an nan/non-finte values:
    ok = np.isfinite(spec) & np.isfinite(spec_err) & np.isfinite(sky)

    num = np.sum((sky[ok] * spec[ok]) / spec_err[ok]**2)
    denom = np.sum((sky[ok] / spec_err[ok])**2)

    return num / denom


























