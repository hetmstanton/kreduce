"""
This file calculates star offsets from good exposures and combines into the final datacubes
"""
from . import exposure

import numpy as np
import subprocess
import astropy.io.fits as fits

import datetime

import glob

import os

def kmos_combine(frames, detector):

    # make the star prameter file:
    starparams_filename, combinefiles_filename = write_star_data_to_file(frames=frames,
      psf_cut=0.8, edge_x=2., edge_y=2., detector=detector)

    # now make the shifts file:
    usershifts_filename = make_user_shifts_file(starparams_filename=starparams_filename,
        detector=detector)

    # get the names of the objects on the detector from one of the frames:
    exp = exposure.Exposure(sci_reconstructed_file=frames[0])
    exp.get_object_names_for_each_detector()
    object_names = exp.objects_per_detector["{:d}".format(detector)]

    for object_name in object_names:
        kmo_comb = 'esorex kmos_combine --method=user --filename={:s} --name={:s} {:s}'.format(usershifts_filename,
            object_name, combinefiles_filename)
        subprocess.call(kmo_comb, shell=True)


def make_user_shifts_file(starparams_filename, detector):
    """C
    reate file of user shifts to use with ``kmos_combine --method="user"``
    
    Note:
        To run ``kmos_combine`` you *must* use the frames in the order they are used here
        i.e. so the offset matches the right frame
    Args:
        starparams_filename (str): filepath of table with the star positions
        usershifts_filename (str, optional): filename to save usershift table to
    """

    # Get info from star table
    star_table = np.genfromtxt(starparams_filename, 
        names=True, skip_header=1)
    n_frames = len(star_table)

    # Calculate shifts and save to file
    print('\n*** [kreduce-combine]: Calculating shifts')
    shiftx, shifty = [0]*(n_frames-1), [0]*(n_frames-1)
    for i in range(1, n_frames):
        shiftx[i-1] = star_table['XCEN_pix'][i] - star_table['XCEN_pix'][0]
        shifty[i-1] = star_table['YCEN_pix'][0] - star_table['YCEN_pix'][i]

    usershifts_filename = starparams_filename.replace('.txt', '_det{:d}_usershifts.txt'.format(detector))

    np.savetxt(usershifts_filename, np.array([shiftx, shifty]).T, fmt='%f', delimiter='    ')
    print('*** [kreduce-combine]: Saved shifts file to', usershifts_filename)

    return usershifts_filename


def write_star_data_to_file(frames, detector, psf_cut=0.8, edge_x=2., edge_y=2.):
    """
    Given a directory containin exposures to combine, will find stars and measure PSFs.
    Then creates a list of frames with PSF FWHM below a given value, saves
    the list of star positions (to calculate user shifts) and creates a .sof
    file of the 'good' frames.

    Note this must be done one detector at a time since the shifts may vary detector
    to detector.
    """

    frames_dir = os.path.abspath(frames[0].split('/SCI_RECON')[0])

    # move into the frames directory to run these functions:
    root_dir = os.getcwd()
    os.chdir(frames_dir)

    combine_files = []
    star_table, star_table_bad = [], []

    datenow = str(datetime.date.today())
    starparams_filename = f'star_psf_{datenow}_det{detector}.txt'
    combinefiles_filename = f'combine_{datenow}_det{detector}.sof'

    for i, frame_full_path in enumerate(frames):

        frame = 'SCI_' + frame_full_path.split('SCI_')[1]

        print('\n*** [kreduce-combine]: Calculating shifts for {:s}'.format(frame))

        # load an exposure object:
        sci_reconstructed = exposure.Exposure(sci_reconstructed_file=frame)

        psf_center_x, psf_center_y, psf_fwhm, psf_ba,\
            psf_pa, invert_comment = star_psf(exposure=sci_reconstructed, detector=detector)

        if psf_fwhm < psf_cut and psf_center_x > edge_x and psf_center_y > edge_y:
            combine_files.append([frame, 'SCI_RECONSTRUCTED'])
            star_table.append([frame, sci_reconstructed.frame_time, psf_center_x, 
                psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])
        elif psf_center_x == np.nan:
            star_table_bad.append([frame, sci_reconstructed.frame_time, psf_center_x, 
                psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])
        else:
            print('WARNING: Something wrong with PSF in {:s}'.format(frame))
            star_table_bad.append([frame, sci_reconstructed.frame_time, psf_center_x, 
                psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])

    # Save star parameter output
    star_table = np.array(star_table)

    FWHM  = np.array([float(x) for x in star_table[:, 4]])
    med_FWHM  = np.nanmedian(FWHM[FWHM > 0.])
    ba  = np.array([float(x) for x in star_table[:, 5]])
    med_ba  = np.nanmedian(ba[ba > 0.])
    pa  = np.array([float(x) for x in star_table[:, 6]])
    med_pa  = np.nanmedian(pa[pa > 0.])

    np.savetxt(starparams_filename, star_table, fmt='%s', delimiter='    ', header='Star PSFs. Median FWHM = %.3f Median BA = %.3f Median PA = %.3f \nFrame    OB_time     XCEN_pix     YCEN_pix    FWHM_arcsec    BA    PA_deg' % (med_FWHM, med_ba, med_pa), comments='# ')
    np.savetxt(starparams_filename.replace('.txt', '_bad.txt'), star_table_bad, fmt='%s', delimiter='    ', header='Bad Star PSFs. \nFrame    OB_time     XCEN_pix     YCEN_pix    FWHM_arcsec    BA    PA_deg', comments='# ')
    print(' - Saved PSF info file to {:s} '.format(starparams_filename))

    # Save .sof file with all the successful frames
    combine_files = np.array(combine_files)
    np.savetxt(combinefiles_filename, combine_files, fmt='%s')
    print(' - Saved combine.sof file to {:s}'.format(combinefiles_filename))

    os.chdir(root_dir)

    return starparams_filename, combinefiles_filename

def star_fit_profile(exposure, detector):
    """
    Fit gaussian profiles to star to find PSFs
    Using ESO command line pipeline tools (``esorex``) extract the star from the exposure,
    fit a Gaussian profile to it, and save the fitted star to a new fits file ``star_file_name``.

    Args:
        exposure (object): exposure object
    Returns:
        star_file_name (str): filepath of FITS file with star with PSF fitted to it
        invert (bool): was the star flux weird and inverted?
    """

    # get the star ifus:
    exposure.get_star_ifus()

    if detector == 1:
        star_ifu = exposure.star_ifu_detector1
    elif detector == 2:
        star_ifu = exposure.star_ifu_detector2
    else:
        star_ifu = exposure.star_ifu_detector3

    # copy the IFU into its own file:
    kmo_copy = 'esorex kmo_copy -x=1 -y=1 -z=1 -xsize=14 -ysize=14 -zsize=2048 -ifu={:s} {:s}'.format(str(star_ifu), 
        exposure.filename)
    subprocess.call(kmo_copy, shell=True)
    status, copyfile = subprocess.getstatusoutput("find . -maxdepth 1 -iname COPY.fits")

    # collapse the IFU to make an image:
    kmo_make_image = 'esorex kmo_make_image %s' % copyfile
    subprocess.call(kmo_make_image, shell=True)
    status, makeimagefile = subprocess.getstatusoutput("find . -maxdepth 1 -iname MAKE_IMAGE.fits")

    # Check the image, if weird, invert (stolen from Charlottes scripts)
    image_hdu = fits.open(makeimagefile)
    image     = image_hdu[1].data
    test_star = np.nansum(image[3:-3, 3:-3] - np.nanmedian(image[3:-3, 3:-3]))
    invert = False
    
    if test_star < 0.:
        print('WARNING: Weird star in {:s}, multiplying image by -1'.format(exposure.filename))
        image_hdu[1].data = -1. * image
        image_hdu.writeto(makeimagefile, clobber=True)
        invert = True

    # Rename star file
    exposure.star_image_file = exposure.filename.strip('.fits') + '_star_image.fits'
    subprocess.call('cp {:s} {:s}'.format(makeimagefile, exposure.star_image_file), 
        shell=True)

    # Fit profile to star image
    kmo_fit_profile = 'esorex kmo_fit_profile {:s}'.format(exposure.star_image_file)
    subprocess.call(kmo_fit_profile, shell=True)
    status, fitprofilefile = subprocess.getstatusoutput("find . -maxdepth 1 -iname FIT_PROFILE.fits")

    # Tidy up
    star_file_name = exposure.filename.strip('.fits')+'_star_psf.fits'
    rename_fit_profile = 'mv {:s} {:s}'.format(fitprofilefile, star_file_name)
    delete_temps = 'rm {:s} {:s}'.format(copyfile, makeimagefile)

    subprocess.call(rename_fit_profile, shell=True)
    subprocess.call(delete_temps, shell=True)

    print('\n*** [kreduce-combine]: Saved star psf profile to {:s}'.format(star_file_name))

    return star_file_name, invert


def star_psf(exposure, detector):
    """
    Get the PSF profiles of the star for calculating user shifts
    Args:
        exposure (object): exposure object
    Returns:
        psf_center_x (float): PSF centroid x in pixels
        psf_center_y (float): PSF centroid y in pixels
        psf_fwhm (float): PSF FWHM in arcsec
        psf_ba (float): PSF BA # TODO I don't remember what this is!!
        psf_pa (float): PSF position angle in degrees
        exposure.invert (bool): is the flux weird?
    """

    exposure.invert = False
    
    exposure.starfile, exposure.invert = star_fit_profile(exposure, detector) 

    if exposure.invert:
        invert_comment = '# weird star, inverted flux'
    else:
        invert_comment = ''

    star_hdulist  = fits.open(exposure.starfile)
    star_hdr      = star_hdulist[1].header
    # Load fit parameters
    psf_center_x = star_hdr['HIERARCH ESO PRO FIT CENTROID X']
    psf_center_y = star_hdr['HIERARCH ESO PRO FIT CENTROID Y']
    psf_r_x      = star_hdr['HIERARCH ESO PRO FIT RADIUS X']
    psf_r_y      = star_hdr['HIERARCH ESO PRO FIT RADIUS Y']

    # Get FWHM, BA, and PA
    psf_fwhm = 2.3548 * 0.5 * (psf_r_x + psf_r_y) * 0.2  # fwhm in arcsec
    psf_ba   = psf_r_y / psf_r_x
    psf_pa   = star_hdr['HIERARCH ESO PRO FIT ROT']
    if psf_pa > 0.:
        psf_pa = psf_pa % 360.
    
    print('\n*** [kreduce-combine]: XPIX={:.2f}, YPIX={:.2f}, FWHM={:.3f} arcsec, BA={:.3f}, PA={:.3f}'.format(psf_center_x, 
        psf_center_y, psf_fwhm, psf_ba, psf_pa))

    return psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment








