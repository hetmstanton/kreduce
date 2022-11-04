# Handles imports of specific files when kreduce is imported into a script
from . import calibration
from . import utils
from . import reconstruct
from . import exposure
from . import sky_correct
from . import combine

import os
kmos_calib_path = os.getenv("KMOS_CALIB_PATH") # "/Users/s1621489/Documents/Package_Installs/KMOS/Calibrations/share/esopipes/datastatic/kmos-4.0.4" 
