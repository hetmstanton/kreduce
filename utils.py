import numpy as np

import glob
import shutil

def add_file_to_sof(sof, file, file_type):
	sof.write(f"{file}\t{file_type}\n")
	return True


def move_reduced_ob_files(ob_dir, destination, skycorr=True):

	if skycorr is True:
		dfiles = glob.glob("{:}/sci/SCI_RECON*-sci_SKYCORR.fits".format(ob_dir))
		for dfile in dfiles:
			print(dfile)
			shutil.copy(src=dfile, dst=destination)
	else:
		dfiles = glob.glob("{:}/sci/SCI_RECON*-sci.fits".format(ob_dir))
		for dfile in dfiles:
			print(dfile)
			shutil.copy(src=dfile, dst=destination)



