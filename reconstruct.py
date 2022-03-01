import kreduce
from . import utils

from astropy.io import fits

import glob

import sys

import os
import shutil

import subprocess

import argparse

def kmos_sci_reduce(ob_dir, task='all'):

	abs_ob_dir = os.path.abspath(ob_dir)

	root_dir = os.getcwd()

	print(f"\n*** [kreduce-sci-red]: running {task.upper()} calibrations for {abs_ob_dir}")

	if task == "all" or task == "sci":

		print(f"\n*** [kreduce-sci-red]: processing SCIENCE frames")

		os.chdir(f"{abs_ob_dir}/sci")

		# if there are already SCI_CONSTRUCTED files in directory, delete them

		sci_files = glob.glob(f"*-sci.fits")

		sof = open(f"sci.sof", 'w')
		for sci in sci_files:
			utils.add_file_to_sof(sof=sof, file=sci, file_type='SCIENCE')

		sof.write(f"{abs_ob_dir}/calib/XCAL_HKHKHK.fits\tXCAL\n")
		sof.write(f"{abs_ob_dir}/calib/YCAL_HKHKHK.fits\tYCAL\n")
		sof.write(f"{abs_ob_dir}/calib/LCAL_HKHKHK.fits\tLCAL\n")
		sof.write(f"{abs_ob_dir}/calib/MASTER_FLAT_HKHKHK.fits\tMASTER_FLAT\n")
		sof.write(f"{abs_ob_dir}/calib/ILLUM_CORR_HKHKHK.fits\tILLUM_CORR\n")
		sof.write(f"{abs_ob_dir}/calib-std/TELLURIC_HKHKHK.fits\tTELLURIC\n")

		sof.write(f"{kreduce.kmos_calib_path}/kmos_wave_band.fits\tWAVE_BAND\n")
		sof.write(f"{kreduce.kmos_calib_path}/kmos_oh_spec_hk.fits\tOH_SPEC\n")
		sof.close()

		process = subprocess.run(["esorex", "kmos_sci_red", 
			"-no_combine", "-background", "-sky_tweak", "-pix_scale=0.1", "sci.sof"], 
	       	stdout=subprocess.PIPE, 
	        universal_newlines=True)

		shutil.move(src="esorex.log", dst="esorex-sci.log")
		os.chdir(root_dir)

		sof.close()

	if task == "all" or task == "acq":

		print(f"\n*** [kreduce-sci-red]: processing ACQUISITION frames")

		os.chdir(f"{abs_ob_dir}/sci")

		# if there are already SCI_CONSTRUCTED files in directory, delete them

		sci_files = glob.glob(f"*-acq.fits")

		sof = open(f"acq.sof", 'w')
		for sci in sci_files:
			utils.add_file_to_sof(sof=sof, file=sci, file_type='SCIENCE')

		sof.write(f"{abs_ob_dir}/calib/XCAL_HKHKHK.fits\tXCAL\n")
		sof.write(f"{abs_ob_dir}/calib/YCAL_HKHKHK.fits\tYCAL\n")
		sof.write(f"{abs_ob_dir}/calib/LCAL_HKHKHK.fits\tLCAL\n")
		sof.write(f"{abs_ob_dir}/calib/MASTER_FLAT_HKHKHK.fits\tMASTER_FLAT\n")
		sof.write(f"{abs_ob_dir}/calib/ILLUM_CORR_HKHKHK.fits\tILLUM_CORR\n")
		sof.write(f"{abs_ob_dir}/calib-std/TELLURIC_HKHKHK.fits\tTELLURIC\n")

		sof.write(f"{kreduce.kmos_calib_path}/kmos_wave_band.fits\tWAVE_BAND\n")
		sof.write(f"{kreduce.kmos_calib_path}/kmos_oh_spec_hk.fits\tOH_SPEC\n")
		sof.close()

		process = subprocess.run(["esorex", "kmos_sci_red", 
			"-no_combine", "-background", "-pix_scale=0.1", "acq.sof"], 
	       	stdout=subprocess.PIPE, 
	        universal_newlines=True)

		shutil.move(src="esorex.log", dst="esorex-acq.log")
		os.chdir(root_dir)

		sof.close()