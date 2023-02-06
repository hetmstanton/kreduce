"""
Adapted from: https://github.com/charlottenosam/kmos_tools/blob/master/kmos_tools/sky_clean.py
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import datetime
import os, sys
import shutil
import glob
import warnings
warnings.filterwarnings("ignore")

def make_sky_residual_spectrum(exposure, show=True, masking=True):
    """Calculate all the sky residual corrections for one DIT ~ a la Trevor Mendel/Charlotte Mason
    
    Generate median 1D sky spectra for each detector (excluding strong emission lines),
    only from galaxy cubes (skip star IFUs), and save them to use for subtraction.

    Note this only works if objects are generally small with respect to the detector size.
    """

    print(f'\n*** [kreduce-sky]: calculating sky residuals in {exposure.filename}\n')

    # need to know which ifus contain stars so that we skip them:
    exposure.get_star_ifus()

    # Create stacks of 'empty' cubes:
    exposure.filename_skyspec = exposure.filename.replace('.fits', '_SKYSPEC.fits')
    detector1, detector2, detector3 = [], [], []

    # determine redshift data 
    if masking: obj_ids, obj_redshifts = systemic_z('./analysis/lineresults/summary/redshift_catalogue.txt') # find_z(filepath=f'{os.getcwd()}/data/ob-data/kvs_final_target_catalogue.txt')

    # Add non-empty, non reference star IFUS to their relevant detector arrays
    for ifu in range(1, 25):

        # Skup IFUs containing stars
        if ifu in exposure.star_ifus:
            print(f'*** [kreduce-sky]: skipping IFU {ifu} (contains ref star)')
            continue

        # Add relevant IFU data to the ext
        ext = exposure.hdulist[f'IFU.{ifu}.DATA']

        # make sure to skip empty ifus
        if len(ext.shape) > 0:

            ifu_comment = f'HIERARCH ESO OCS ARM {ifu} NAME'
            ifu_header  = ext.header
            ifu_cube = ext.data

            if masking:

                # get redshifts
                redshift = obj_redshifts[obj_ids==exposure.hdr[f'HIERARCH ESO OCS ARM{ifu} NAME']][0]

                # generate wavelengths
                wl0 = ifu_header['CRVAL3']
                dwl = ifu_header['CDELT3']
                naxis = ifu_cube.shape[0]
                wlf = wl0 + (dwl*naxis)
                wavelengths = np.arange(wl0, wlf, dwl) * 1e4 #μm -> Å

                # convert from micrometers to anstroms
                wavelengths *= 1e4

                # make and apply mask
                cube_mask = mask_ifu_cube(ifu_cube, redshift.astype(float), wavelengths)
                m_ifu_cube = np.copy(ifu_cube)
                m_ifu_cube[cube_mask] = np.nan

            if 1 <= ifu <= 8:
                detector1.append(m_ifu_cube)
            elif 9 <= ifu <= 16:
                detector2.append(m_ifu_cube)
            else:
                detector3.append(m_ifu_cube)

    # Set values as arrays and work out cumulative length
    len_for_stack = len(detector1) + len(detector2) + len(detector3)

    # Set Detectors as arrays
    detector1, detector2, detector3 = np.array(detector1), np.array(detector2), np.array(detector3)

    # Combine all detectors into one overall array
    detector_all = np.concatenate((detector1, detector2, detector3), axis=0)
    # Handle nans, take median of all detectors, setting any nans to zero
    skyspec_1D_all = np.nanmedian(detector_all, axis=(0, 2, 3))
    skyspec_1D_all[np.isnan(skyspec_1D_all)] = 0.0

    # Work out individual medians (with nans set to 0) for each individual detector
    skyspec_1D = {}
    detectors = [detector1, detector2, detector3]
    for i in range(len(detectors)):
        # Take median 
        detector_sky = np.nanmedian(detectors[i], axis=(0, 2, 3))
        detector_sky[np.isnan(detector_sky)] = 0
        skyspec_1D[i] = detector_sky

    # Plotting
    if show:

        plt.figure(figsize=(10, 5))

        # Plot overall sky spectrum of all detectors
        plt.plot(skyspec_1D_all, lw=1, alpha=0.8, 
            label=f'All detectors ({detector_all.shape[0]} IFUs)', zorder=10)

        # Plot individual Detector Sky Spectra
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


    # Write data for sky spectra out to its own fits file, for each detector and the overall median
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

    """ subtracts the residual sky from each ifu to yield the sky-subtracted cubes. 
        + also produces any relevant plots. """

    print(f'\n*** [kreduce-sky]: subtracting sky residuals in {exposure.filename}\n')

    # Open spectrum corresponding to sky residual spectra
    exposure.filename_skycorr = exposure.filename.replace('.fits', '_SKYCORR.fits')
    skyspec_all = fits.open(exposure.filename_skyspec)

    # need to know which ifus contain stars so that we skip them:
    exposure.get_star_ifus()

    # Make corresponding directory to store sky_corrected plots if applicable.
    if show:
        if not os.path.exists(f'{exposure.filedir}/sky_correct_plots'):
            os.mkdir(f'{exposure.filedir}/sky_correct_plots')

    # Iterate over IFUs
    for ifu in range(1, 25):

        # skip ifus containing stars
        if ifu in exposure.star_ifus:
            continue

        # generate ext object based on the ifu fits data
        ext = exposure.hdulist[f'IFU.{ifu}.DATA']

        # make sure to skip empty ifus
        if len(ext.shape) > 0:

            # select corresponding detector
            if 1 <= ifu <= 8:
                detector = 1
            elif 9 <= ifu <= 16:
                detector = 2
            elif 17 <= ifu <= 24:
                detector = 3

            # collect sky data
            sky  = skyspec_all[detector+1].data

            # Estimate 1d error from std of flux
            ext_copy = ext.data.copy()
            # calculate corresponding error spectrum
            error1d = np.nanstd(ext_copy, axis=(1, 2))/np.sqrt(ext_copy.shape[1] * ext_copy.shape[2])

            if show:
                
                # Simply a fancy way of producing subplots of the correct format
                figure_mosaic = """
                ACC
                BCC
                """
                fig, axes = plt.subplot_mosaic(mosaic=figure_mosaic, figsize=(11,5))

                collapsed_original = np.nanmedian(ext_copy, axis=(1,2))
                # Plot original spectrum with sky overplotted
                axes["A"].plot(collapsed_original, lw=1.,
                    alpha=0.8, label='Original full spectrum', color='k')
                axes["A"].plot(sky, lw=1.,
                    alpha=0.8, label='Sky', color='C0', ls=':')
                ymin_sky, ymax_sky = np.nanpercentile(collapsed_original, 1), np.nanpercentile(collapsed_original, 99)
        
                if ymin_sky > 0.: ymin_sky = -1.e-18

                axes["A"].set_ylim(ymin_sky, ymax_sky)
                axes["A"].set_ylabel('Flux')
                axes["A"].legend()

                # Plot black histogram of flux values pre-subtraction
                nanmask = np.isfinite(collapsed_original)
                xmin_hist, xmax_hist = np.nanpercentile(collapsed_original[nanmask], 1), np.nanpercentile(collapsed_original[nanmask], 99)
                limmask = (collapsed_original[nanmask] >= xmin_hist) & (collapsed_original[nanmask] <= xmax_hist) 
                collapsed_to_plot = collapsed_original[nanmask][limmask]

                # TODO 
                label = f'nocorr'

                axes["C"].hist(collapsed_to_plot, histtype='step', color='k', lw=1.,
                    bins=100, label=label)
                axes["C"].set_xlabel('Flux')

            # make a copy of the original data cube and an 
            # empty array that will hold the sky-subtracted copy:
            data_orig = np.copy(ext.data)
            data_corr = np.empty_like(data_orig)

            # Iterate over every pixel
            for i in range(data_orig.shape[1]):
                for j in range(data_orig.shape[2]):

                    # Collect data from the given pixel
                    pixel_spec = data_orig[:,i,j]

                    if np.isnan(pixel_spec).all():
                        # If all the pixels are nan, set the data correction to an array of nans
                        data_corr[:,i,j] = np.zeros_like(data_orig[:,i,j]) * np.nan
                    else:
                        # get sky normalisation given the pixel spectrum, error spectrum and sky.
                        sky_norm =  get_sky_normalization(spec=pixel_spec, 
                            spec_err=error1d, sky=sky)
                        # get the data correction by subtracting the product of the sky normalisation and the sky
                        data_corr[:,i,j] = data_orig[:,i,j] - sky_norm*sky

            # Plotting
            if show:
                # plot the sky corrected spectrum
                collapsed_skycorr = np.nanmedian(data_corr, axis=(1,2))
                axes["B"].plot(collapsed_skycorr, lw=1.,
                    alpha=0.8, label='Sky corrected spectrum', color='C3')

                axes["B"].set_ylabel('Flux')
                axes["B"].set_xlabel('Wavelength [pix]')

                axes["B"].legend()
                axes["B"].set_ylim(ymin_sky, ymax_sky)

                # plot the sky subtracted flux histogram
                nanmask = np.isfinite(collapsed_skycorr)
                xmin_hist, xmax_hist = np.nanpercentile(collapsed_skycorr[nanmask], 1), np.nanpercentile(collapsed_skycorr[nanmask], 99)
                limmask = (collapsed_skycorr[nanmask] >= xmin_hist) & (collapsed_skycorr[nanmask] <= xmax_hist) # get mask between the limits 
                collapsed_to_plot = collapsed_skycorr[nanmask][limmask] # mask out nans and values outside the limits
                # TODO Label
                label = f'skycorr'
                axes["C"].hist(collapsed_to_plot, histtype='step', color='C3', lw=1.,
                    bins=100, label=label)

                axes["C"].set_xlabel('Flux')
                axes["C"].legend()
                if True: axes["C"].set_title('Masked')


                plt.tight_layout()
                
                plot_filename = exposure.filename_skycorr.replace('.fits', f'_IFU{ifu}.pdf')
                fig.savefig(f'{plot_filename}')

                # remove the plot if it already exists:
                root_pfile = f'SCI_' + plot_filename.split('SCI_')[1]
                if os.path.exists(f"{exposure.filedir}/sky_correct_plots/{root_pfile}"):
                    os.remove(f"{exposure.filedir}/sky_correct_plots/{root_pfile}")

                shutil.move(src=os.path.abspath(plot_filename), 
                    dst=os.path.abspath(f'{exposure.filedir}/sky_correct_plots'))

            # set the exposure data to the corrected sky spectrum 
            ext.data = data_corr

        # Write new sky corrected spectrum out
        exposure.hdulist.writeto(exposure.filename_skycorr, overwrite=True)


def get_sky_normalization(spec, spec_err, sky):

    """ calculates the normalisation of the spectrum given the error spectrum and the sky spectrum """
    # mask an nan/non-finte values:
    ok = np.isfinite(spec) & np.isfinite(spec_err) & np.isfinite(sky)

    # Multiply sky values by spec values and divide by error
    num = np.sum((sky[ok] * spec[ok]) / spec_err[ok]**2)
    denom = np.sum((sky[ok] / spec_err[ok])**2)
    # Return normalistion
    return num / denom

def mask_ifu_cube(cube:object, redshift:float, wavelengths:object) -> object:
    """
    this method takes an ifu cube, and determines its relevant redshift from a given data file. It then calcualtes the positions
    of specific stellar lines [OII 3727, Hβ , OIII 4969, OIII 5007] within the spectral indices of the data cube and masks them out.
    The masked cube is then returned.
    Arguments:  cube (object)
    Returns:    masked_cube (object)
    """

    # read in line data
    linedata = np.genfromtxt(f'{os.getcwd()}/data/ob-data/linedata.txt', dtype=str, unpack=True)
    _, lines = linedata[0], linedata[1].astype(float)

    # mask arrays
    maskrange = np.ones(shape=(2, lines.shape[0]))

    for i, linewl in enumerate(lines):

        # scale line to redshift
        scaled_wl = (1+redshift)*linewl

        # determine line width under assumption lines are around 200kms-1 wide
        halfwidth = find_line_width(scaled_wl, v=200) / 2

        # generate masking values
        maskrange[0,i], maskrange[1,i] = scaled_wl-halfwidth, scaled_wl+halfwidth

    # generate mask
    mask = np.array([np.any((maskrange[0] <= wl) & (wl <= maskrange[1])) for wl in wavelengths])

    return mask

def find_line_width(wl:float, v:int) -> float:
    """ return delta lambda based on  """
    return (v/3e5)*wl

def find_z(filepath:str) -> tuple:
    """
    Open file, find corresponding value, return z
    """
    # extract object ids and redshifts
    data = np.genfromtxt(filepath, dtype=str, usecols=[0,2], unpack=True)
    return (data[0], data[1])

def systemic_z(path:str) -> tuple:
    redshift_tabl = np.genfromtxt(path, dtype=str)
    obj_names = redshift_tabl.T[0]
    sys_z = redshift_tabl.T[3].astype(float)
    return (obj_names, sys_z)

























