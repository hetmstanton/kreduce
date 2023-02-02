"""
Adapted from: https://github.com/charlottenosam/kmos_tools/blob/master/kmos_tools/exposure.py
"""

import astropy.io.fits as fits
import os
class Exposure:

    def __init__(self, sci_reconstructed_file):

        self.filename = sci_reconstructed_file
        self.filedir = os.path.abspath(self.filename.split('/SCI_')[0])

        # Read in fits data 
        self.hdulist  = fits.open(self.filename)
        self.hdr      = self.hdulist[0].header
        self.filter   = self.hdr['HIERARCH ESO INS FILT1 ID']

        self.frame_time = self.hdr['DATE-OBS']

    def get_star_ifus(self):
        """
        All objects in the catalogues have a name: KVS_{number}
        If number >= 430 this is not a science target and therefore
        will be one of the stars placed in each dectector.

        Note we skip ifu 14 because it was not used on our
        program.
        """

        self.star_ifus = []

        for ifu in range(1, 25, 1): # Run through the 24 arms 

            ext = self.hdulist[f'IFU.{ifu}.DATA']
            
            # make sure to skip empty IFUs (Integral Field Units?)
            if len(ext.shape) > 0:

                ifu_hdr = self.hdulist[f'IFU.{ifu}.DATA'].header
                name = ifu_hdr[f'HIERARCH ESO OCS ARM{ifu} NAME']

                number = int(name.split('_')[1])

                # all numbers > 430 are not science targets:
                if number >= 430:

                    self.star_ifus.append(ifu)

                    # useful to also separate by detector:
                    if 1 <= ifu <= 8:
                        self.star_ifu_detector1 = ifu
                    elif 9 <= ifu <= 16:
                        self.star_ifu_detector2 = ifu
                    else:
                        self.star_ifu_detector3 = ifu

    def get_object_names_for_each_detector(self):
        """
        Returns a dictonary with a list for each detector
        specifiing the unique indentifying names of the objects
        on each detector.

        This is useful for kmos_combine, when we want to combine
        only objects on one detector to ensure accurate shifts.
        """

        self.objects_per_detector = {'1': [],
                                     '2': [],
                                     '3': []}

        self.valid_ifus_per_detector = {'1': [],
                                        '2': [],
                                        '3': []}

        for ifu in range(1, 25, 1): # Iterate through all IFUs

            ext = self.hdulist[f'IFU.{ifu}.DATA']

            if len(ext.shape) > 0:

                ifu_hdr = self.hdulist[f'IFU.{ifu}.DATA'].header # Get Header
                name = ifu_hdr[f'HIERARCH ESO OCS ARM{ifu} NAME']# Get name of object

                # Append object name to the appropriate detector list.
                if 1 <= ifu <= 8:
                    self.objects_per_detector['1'].append(name)
                    self.valid_ifus_per_detector['1'].append(ifu)
                elif 9 <= ifu <= 16:
                    self.objects_per_detector['2'].append(name)
                    self.valid_ifus_per_detector['2'].append(ifu)
                else:
                    self.objects_per_detector['3'].append(name)
                    self.valid_ifus_per_detector['3'].append(ifu)
  


