""" This file calculates star offsets from good exposures and combines into the final datacubes """
from . import exposure

from astropy.table import Table
import numpy as np
import subprocess
import astropy.io.fits as fits
import datetime
import time
import glob
import os, sys
import matplotlib.pyplot as plt
from matplotlib import patches as pat

def kmos_combine(frames, detector, combinefiles_filename, usershifts_filename):
    """ combines all the exposures of a given object into a combined final cube """

    # get the names of the objects on the detector from one of the frames and set.
    exp = exposure.Exposure(sci_reconstructed_file=frames[0])
    exp.get_object_names_for_each_detector()
    object_names = exp.objects_per_detector["{:d}".format(detector)]

    # Run Esorex kmos_combine method
    for object_name in object_names:

        kmo_comb = 'esorex kmos_combine --method=user --filename={:s} --name={:s} {:s}'.format(usershifts_filename,
            object_name, combinefiles_filename)
        subprocess.call(kmo_comb, shell=True)

        # check combination for zero cubes 
        while True:

            # check final cube to make sure its not zero
            cube_path = f'{os.getcwd()}/COMBINE_SCI_RECONSTRUCTED_{object_name}.fits'
            hdul = fits.open(cube_path)
            data = hdul[1].data

            checksum = np.sum(np.nansum(data, axis=0))
            if checksum != 0:
                break
            print(f'-> [kreduce-combine]: {object_name} failed to combine properly. Re-running... ')
            with open('./combination_log.txt', 'a') as file:
                file.write(f'Re-combining {object_name} as combination bugged...')
            subprocess.call(kmo_comb, shell=True)

    # Process produces a COMBINED_CUBE, COMBINED_IMAGE
    # and EXP_MASK - exposure time frame, every spaxel indicating how many input frames taken into account when combining.


def kmos_calculate_shifts(frames:list, detector:int, edge_x:float=2., edge_y:float=2., shifts_only:bool=False):
    """ creates the files containing the user shifts, first by creating a starparams_file and combinefiles using 
    the internal routine write_star_data_to_file.
    """
    # make the star parameter file:
    starparams_filename, combinefiles_filename = write_star_data_to_file(frames=frames,
      psf_cut=0.8, edge_x=edge_x, edge_y=edge_y, detector=detector)

    # now make the user shifts file that is fed into kmos_combine:
    usershifts_filename = make_user_shifts_file(starparams_filename=starparams_filename,
        detector=detector, frames=frames)

def make_user_shifts_file(frames:list, starparams_filename:str, detector:int):
    """
    create file of user shifts to use with ``kmos_combine --method="user"``
    
    Note:
        To run ``kmos_combine`` you *must* use the frames in the order they are used here
        i.e. so the offset matches the right frame
    Args:
        starparams_filename (str): filepath of table with the star positions
        usershifts_filename (str, optional): filename to save usershift table to
    """

    # Get info from star table and determine the number of frames
    star_table = np.genfromtxt(starparams_filename, dtype=None,
        names=True, skip_header=1)
    framenames = star_table['Frame']
    star_table = np.genfromtxt(starparams_filename,
        names=True, skip_header=1)
    n_frames = len(star_table)

    # Calculate shifts and save to file
    print('\n*** [kreduce-combine]: Calculating shifts')
    shiftx, shifty = [0]*(n_frames-1), [0]*(n_frames-1)
    for i in range(1, n_frames):
        shiftx[i-1] = star_table['XCEN_pix'][i] - star_table['XCEN_pix'][0]
        shifty[i-1] = star_table['YCEN_pix'][0] - star_table['YCEN_pix'][i]

    # Shifts calculated here so this is what
    all_shiftx = np.append([0], shiftx)
    all_shifty = np.append([0], shifty)
    seeing = np.zeros(all_shiftx.shape[0])
    alldecdata = np.array([framenames, star_table['FWHM_arcsec'], seeing, star_table['XCEN_pix'],
        star_table['YCEN_pix'], all_shiftx, all_shifty], dtype='str') # changed str from U40

    # generate relevant filenaes
    usershifts_filename = starparams_filename.replace('.txt', '_usershifts.txt'.format(detector))
    detectorshifts_filename = starparams_filename.replace('.txt', '_detectorshifts.txt'.format(detector))

    np.savetxt(usershifts_filename, np.array([shiftx, shifty]).T, fmt='%f', delimiter='    ')
    print('*** [kreduce-combine]: Saved shifts file to', usershifts_filename)

    np.savetxt(detectorshifts_filename, alldecdata.T, fmt='%s', delimiter='    ',
        header='Frame    FWHM_arcsec     Seeing     XCEN_pix     YCEN_pix    XSHIFT_pix    YSHIFT_pix')
    print('*** [kreduce-combine]: Saved full detector shift file to', detectorshifts_filename)

    return usershifts_filename


def write_star_data_to_file(frames:list, detector:int, psf_cut:float=0.8, edge_x:float=6., edge_y:float=6.):
    """
    Given a directory containin exposures to combine, will find stars and measure Point spread functions (PSF). 
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

    # Set up empty lists
    combine_files = []
    star_table, star_table_bad = [], []

    datenow = str(datetime.date.today()) # Names files - useful for testing i guess
    # Make directories for these? is this useful to have or should I just yeet them into the aether
    if not os.path.exists(f'./star_psf'):
        os.mkdir(f'./star_psf')
    if not os.path.exists(f'./combine'):
        os.mkdir(f'./combine')

    # Set file names for accessing
    starparams_filename = f'./star_psf/star_psf_det{detector}.txt'
    combinefiles_filename = f'./combine/combine_det{detector}.sof'

    for i, frame_full_path in enumerate(frames):

        # Get OB Name
        ob = int(frame_full_path.split('.')[2][2:])        

        # Split SCI_ part of name out to get the frame ID
        frame = 'SCI_' + frame_full_path.split('SCI_')[1]

        print('\n*** [kreduce-combine]: Calculating shifts for {:s}'.format(frame))

        # load an exposure object:
        sci_reconstructed = exposure.Exposure(sci_reconstructed_file=frame)

        # Get the parameters from the star_psf: (see star_psf for descriptions of the parameters)
        psf_center_x, psf_center_y, psf_fwhm, psf_ba,\
            psf_pa, invert_comment = star_psf(exposure=sci_reconstructed, ob=ob, detector=detector, dither=(i%4)) 

        # How to handle star point spread function based on parameters
        # What does this top one do? #
        if psf_ba == 1.0 and psf_pa == 0.0:
            star_table_bad.append([frame, sci_reconstructed.frame_time, psf_center_x, 
                psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])
            #with open('../textfiles/Oops.txt', 'a') as file:
            #    file.write(f'{frame} got psf_ba = 1.0, psf_pa = 0.0.\n{psf_center_x} {psf_center_y} {psf_fwhm} {psf_ba} {psf_pa}\n\n')
        # If FWHM is less than the cut off and the center is bigger than the edge cut off
        elif psf_fwhm < psf_cut and psf_center_x > edge_x and psf_center_y > edge_y:
            combine_files.append([frame, 'SCI_RECONSTRUCTED'])
            star_table.append([frame, sci_reconstructed.frame_time, psf_center_x, 
                psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])
        # If Nan, append to the bad file
        elif psf_center_x == np.nan:
            star_table_bad.append([frame, sci_reconstructed.frame_time, psf_center_x, 
                psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])
            #with open('../textfiles/Oops.txt', 'a') as file:
            #    file.write(f'{frame} got psf_center_x=nan.\n{psf_center_x} {psf_center_y} {psf_fwhm} {psf_ba} {psf_pa}\n\n')
        # Handle other bad frames
        else:
            print('WARNING: Something wrong with PSF in {:s}'.format(frame))
            #with open('../textfiles/Oops.txt', 'a') as file:
            #    file.write(f'{frame} had something else wrong.\n{psf_center_x} {psf_center_y} {psf_fwhm} {psf_ba} {psf_pa}\n\n')
            star_table_bad.append([frame, sci_reconstructed.frame_time, psf_center_x, 
                psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment])

    # Save star parameter output
    star_table = np.array(star_table)

    # Work out medians of the FWHM, ba and pa
    FWHM  = np.array([float(x) for x in star_table[:, 4]])
    med_FWHM  = np.nanmedian(FWHM[FWHM > 0.])
    ba  = np.array([float(x) for x in star_table[:, 5]])
    med_ba  = np.nanmedian(ba[ba > 0.])
    pa  = np.array([float(x) for x in star_table[:, 6]])
    med_pa  = np.nanmedian(pa[pa > 0.])

    # Save psf data out
    np.savetxt(starparams_filename, star_table, fmt='%s', delimiter='    ', header='Star PSFs. Median FWHM = %.3f Median BA = %.3f Median PA = %.3f \nFrame    OB_time     XCEN_pix     YCEN_pix    FWHM_arcsec    BA    PA_deg' % (med_FWHM, med_ba, med_pa), comments='# ')
    np.savetxt(starparams_filename.replace('.txt', '_bad.txt'), star_table_bad, fmt='%s', delimiter='    ', header='Bad Star PSFs. \nFrame    OB_time     XCEN_pix     YCEN_pix    FWHM_arcsec    BA    PA_deg', comments='# ')
    print(' - Saved PSF info file to {:s} '.format(starparams_filename))

    # Save .sof file with all the successful frames
    combine_files = np.array(combine_files)
    np.savetxt(combinefiles_filename, combine_files, fmt='%s')
    print(' - Saved combine.sof file to {:s}'.format(combinefiles_filename))

    os.chdir(root_dir)

    return starparams_filename, combinefiles_filename

def star_fit_profile(exposure:object, detector:int, edge_cut:int=3):
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

    # Set IFU corresponding to given detector
    if detector == 1:
        star_ifu = exposure.star_ifu_detector1
    elif detector == 2:
        star_ifu = exposure.star_ifu_detector2
    else:
        star_ifu = exposure.star_ifu_detector3

    # copy the IFU into its own file:
    kmo_copy = 'esorex kmo_copy -x=1 -y=1 -z=1 -xsize=28 -ysize=28 -zsize=2048 -ifu={:s} {:s}'.format(str(star_ifu), 
        exposure.filename)
    subprocess.call(kmo_copy, shell=True)
    status, copyfile = subprocess.getstatusoutput("find . -maxdepth 1 -iname COPY.fits")

    # collapse the IFU to make an image:
    kmo_make_image = 'esorex kmo_make_image %s' % copyfile
    subprocess.call(kmo_make_image, shell=True)
    status, makeimagefile = subprocess.getstatusoutput("find . -maxdepth 1 -iname MAKE_IMAGE.fits")

    # Check the image, if weird, invert (See Charlotte Mason's scripts)
    image_hdu = fits.open(makeimagefile)
    image     = image_hdu[1].data
    test_star = np.nansum(image[3:-3, 3:-3] - np.nanmedian(image[3:-3, 3:-3]))
    invert = False
    
    if test_star < 0.:
        print('WARNING: Weird star in {:s}, multiplying image by -1'.format(exposure.filename))
        image_hdu[1].data = -1. * image
        image_hdu.writeto(makeimagefile, clobber=True)
        invert = True

    # Remove edge pixels
    edgeless_image = np.copy(image)
    edgeless_image[edge_cut:-edge_cut, edge_cut:-edge_cut] = 0.0
    image_hdu[1].data -= edgeless_image
    image_hdu.writeto(makeimagefile, clobber=True)

    # Rename star file
    exposure.star_image_file = exposure.filename.strip('.fits') + '_star_image.fits'
    subprocess.call('cp {:s} {:s}'.format(makeimagefile, exposure.star_image_file), 
        shell=True)

    # Fit profile to star image
    kmo_fit_profile = 'esorex kmo_fit_profile {:s}'.format(exposure.star_image_file)
    subprocess.call(kmo_fit_profile, shell=True)
    status, fitprofilefile = subprocess.getstatusoutput("find . -maxdepth 1 -iname FIT_PROFILE.fits")

    # Tidy up files (set star file name to correct format, remove image file and previous fit)
    star_file_name = exposure.filename.strip('.fits')+'_star_psf.fits'
    rename_fit_profile = 'mv {:s} {:s}'.format(fitprofilefile, star_file_name)
    delete_temps = 'rm {:s} {:s}'.format(copyfile, makeimagefile)

    subprocess.call(rename_fit_profile, shell=True)
    subprocess.call(delete_temps, shell=True)

    print('\n*** [kreduce-combine]: Saved star psf profile to {:s}'.format(star_file_name))

    return star_file_name, invert, image


def star_psf(exposure:object, detector:int, ob:int, dither:int, edge_cut:int=2.):
    """
    Get the PSF profiles of the star for calculating user shifts
    Args:
        exposure (object): exposure object
    Returns:
        psf_center_x (float): PSF centroid x [pixels]
        psf_center_y (float): PSF centroid y [pixels]
        psf_fwhm (float): PSF FWHM [arcsec]
        psf_ba (float): PSF BA Ellipticity [dimensionless]
        psf_pa (float): PSF position angle [degrees]
        exposure.invert (bool): is the flux weird?
    """

    exposure.invert = False
    # get the star file name and the invert value by fitting the profile to the star
    exposure.starfile, exposure.invert, image = star_fit_profile(exposure, detector) 

    # set relevant comment 
    invert_comment = '# weird star, inverted flux' if exposure.invert else ''

    # open star data and collect relevant information
    star_hdulist  = fits.open(exposure.starfile)
    star_hdr      = star_hdulist[1].header
    # Load fit parameters (star center and radius)
    psf_center_x = star_hdr['HIERARCH ESO PRO FIT CENTROID X']
    psf_center_y = star_hdr['HIERARCH ESO PRO FIT CENTROID Y']
    psf_r_x      = star_hdr['HIERARCH ESO PRO FIT RADIUS X']
    psf_r_y      = star_hdr['HIERARCH ESO PRO FIT RADIUS Y']

    # Insert Plotting Method Here
    show = True
    if show: plot_star_psf(exposure, detector, ob, dither, image)

    # Get FWHM, BA, and PA
    arcsec_per_pix = 0.1 # only the case for dithering, set some parameter to check 
    psf_fwhm = 2.3548 * 0.5 * (psf_r_x + psf_r_y) * arcsec_per_pix  # fwhm in arcsec [essentially the seeing, star width ~ seeing of obs]
    psf_ba   = psf_r_y / psf_r_x # ratio of y radius to x radius? [ellipticity]
    psf_pa   = star_hdr['HIERARCH ESO PRO FIT ROT'] # Rotation Angle
    # If angle is greater than zero, calculate the modulus by 360
    if psf_pa > 0.: # TODO: What about if the value is negative? what then? still modulo?
        psf_pa %= 360.
    elif psf_pa < 0:
        psf_pa %= (-360)
    
    print('\n*** [kreduce-combine]: XPIX={:.2f}, YPIX={:.2f}, FWHM={:.3f} arcsec, BA={:.3f}, PA={:.3f}'.format(psf_center_x, 
        psf_center_y, psf_fwhm, psf_ba, psf_pa))

    # return parameters
    return psf_center_x, psf_center_y, psf_fwhm, psf_ba, psf_pa, invert_comment

def plot_star_psf(exposure: object, detector: int, ob: int, dither: int, image_data):
    
    # Get Star Data
    star_hdulist  = fits.open(exposure.starfile)
    star_hdr      = star_hdulist[1].header
    # Load fit parameters (star center and radius)
    psf_center_x = star_hdr['HIERARCH ESO PRO FIT CENTROID X']
    psf_center_y = star_hdr['HIERARCH ESO PRO FIT CENTROID Y']
    psf_r_x      = star_hdr['HIERARCH ESO PRO FIT RADIUS X']
    psf_r_y      = star_hdr['HIERARCH ESO PRO FIT RADIUS Y']
    psf_pa   = star_hdr['HIERARCH ESO PRO FIT ROT'] # Rotation Angle
    # If angle is greater than zero, calculate the modulus by 360
    if psf_pa > 0.:
        psf_pa = psf_pa % 360.

    # Produce Star Image
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.imshow(image_data, cmap='plasma', origin='lower', vmin=np.percentile(image_data, 2), vmax=np.percentile(image_data, 98))
    plt.colorbar()
    # Overplot Parameters
    plt.plot(psf_center_x-1, psf_center_y-1, 'x', color='black', label='center')
    ellipse = pat.Ellipse((psf_center_x-1,psf_center_y-1), 2*psf_r_x, 2*psf_r_y, angle=psf_pa, color='black', label='ellipse', fc=(0,0,0,0))
    ax.add_patch(ellipse)
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')

    # Get Star IFU ID
    if detector == 1:
        ifu = exposure.star_ifu_detector1
    if detector == 2:
        ifu = exposure.star_ifu_detector2
    if detector == 3:
        ifu = exposure.star_ifu_detector3

    FileID = int(exposure.filename.split('.')[2].split('-')[0])
    ObjName = exposure.hdr[f'HIERARCH ESO OCS ARM{ifu} NAME']
    bsc = 'SKYCORR' if 'SKYCORR' in exposure.starfile else 'BASIC'

    # Save Figure
    savename = f'../../../star_psf/{bsc.lower()}/d{detector}/EC_{ObjName}_OB{ob}_{FileID:04}_psfplot.png'
    title = f'{ObjName} [{bsc}] [D:{detector} OB:{ob:02} FID:{FileID:04}]'
    plt.title(title)
    plt.savefig(savename)

    # Comment if weird
    with open(f'../../../star_psf/{bsc.lower()}/d{detector}/WeirdStars.txt', 'a') as file:
        file.write(f'{exposure.filename} | Obj: {ObjName} | OB: {ob:02} | Weird: {exposure.invert}\n')




