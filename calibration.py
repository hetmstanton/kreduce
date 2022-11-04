import kreduce
from . import utils
from astropy.io import fits 
import glob
import sys
import os
import shutil
import subprocess 
import argparse 

def rename_raw_data(ob_dir, verbose=True):
	"""

	This function sets the absolute path of the directory and generates the relevant subdirectories.
	Next any .Z files are unzipped. fits files are read in, and based on the object type in their directory
	are assigned a relevant subdirectory. The files are then moved to their relevant subdirectories and any 
	remaining files are moved to an ancillary directory.

	PARAMETERS:
	ob_dir - directory containing all the raw files from one
			 OB downloaded from the ESO archive.

	verobose - switch on/off certain print statements
	"""

	# the absolute path - basically the entirety of the path, no short path:
	abs_ob_dir = os.path.abspath(ob_dir)

	# make the directories that the files will be sorted into (if they do not already exist):
	for dirname in ['ancil', 'calib', 'calib-std', 'sci']:
		if not os.path.exists(f"{abs_ob_dir}/{dirname}"):
			os.mkdir(f'{abs_ob_dir}/{dirname}')

	# first unzip the .Z files
	dotz_files = glob.glob(f"{abs_ob_dir}/*.Z")

	print("\n*** [kreduce-calibration]: unzipping ESO .Z files")
	# Quits the code if there are no files to unzip.
	if len(dotz_files) == 0:
		print(f"no .Z ESO piplines files found in {abs_ob_dir}")
		sys.exit()

	for dotz in dotz_files:
		if verbose:
			print(f"{dotz} --> {dotz.split('.Z')[0]}")
		# Run subprocesses (using command line inputs) to unzip the files in the dotz list
		process = subprocess.run(['gzip', '-d', f'{dotz}'], # Runs the command line using the given arguments
	                         stdout=subprocess.PIPE, # Handles.. I'm not sure
	                         universal_newlines=True) # Same as above

	# List of files in absolute directory (all fits files)
	dfiles = glob.glob(f"{abs_ob_dir}/*.fits")

	# Rename the files and determine their relevant subdirectories
	print("\n*** [kreduce-calibration]: renaming and moving core data files")
	for dfile in dfiles:

		# Simply takes the file ID from the directory
		file_id = dfile.split(f"{abs_ob_dir}/")[1].split('.fits')[0]

		hdu = fits.open(dfile)
		header = hdu[0].header

		if not "OBJECT" in header:
			print(f"missing OBJECT keyword in {dfile}")
			continue

		obj = header["OBJECT"].lower().replace(",","-").strip()
		
		# Determine relevant subdirectory for a given object type
		if obj == "object":
			obj = "acq"
			subdir = 'sci'
		elif obj == "dark" or obj == "flat-off" or obj == "flat-lamp" or obj == "wave-off" or \
			obj == "wave-lamp" or obj == "flat-sky" or obj == "sky":
			subdir = 'calib'
		elif obj == "object-sky-std-flux":
			subdir = 'calib-std'
		else:
			# must be a science frame, but just make sure:
			if not "HIERARCH ESO TPL ID" in header:
				print(f"missing HIERARCH ESO TPL ID keyword in {dfile}")

			tpl = header["HIERARCH ESO TPL ID"].lower().strip()
			
			if not tpl.startswith("kmos_spec_obs"):
				print("unknown frame type {obj} with template {tpl}")
				continue

			obj = "sci"
			subdir = "sci"

		if verbose:
			print(f"{dfile} --> {abs_ob_dir}/{subdir}/{file_id}-{obj}.fits")
		# Move relevant file into the right location based on determined subdir from object type
		shutil.move(src=dfile, dst=f"{abs_ob_dir}/{subdir}/{file_id}-{obj}.fits")

	print("\n*** [kreduce-calibration]: cleaning up, move remaining files to ancil/")
	# move erverything remaining into the ancilliary directory: (e.g. any other file types or remaining files)
	for ext in ['.fits', '.xml', '.txt']:
		ancil_files = glob.glob(f"{abs_ob_dir}/*{ext}")
		if len(dotz_files) == 0:
			print(f"no ancilliary {ext} files found in {abs_ob_dir}")
		else:
			for ancil_file in ancil_files:
				file_name = ancil_file.split(f"{abs_ob_dir}/")[1]
				print(f"{ancil_file} --> {abs_ob_dir}/ancil/{file_name}")
				shutil.move(src=ancil_file, dst=f"{abs_ob_dir}/ancil/{file_name}")


def process_calibrations(ob_dir, task='all'):
	"""

	Handles all the calirations for the defined tasks.
	Each one changes to the relevant directory, opens the relevant sof file,  and adds the relevant files to the sof file.
	An esorex subprocess is then run on the .sof file to carry out the relevant calibration. 
	The esorex log file is then stored and renamed based on which task was carried out.

	ob_dir - directory containing renamed files from ESO archive
			 note: the above rename_raw_data() function must be run
			 before this function is called.

	task - the vaild inputs are 'all', 'dark', 'flat', 'wave_cal',
	       'illum', 'std_star'. allows the user to select specific
	       parts of the calibration to run, or to run all at once.
	"""

	abs_ob_dir = os.path.abspath(ob_dir) # Absolute path
	root_dir = os.getcwd() # Current working directory

	# Basically ensures the first function has been run so that files are in the correct organisation
	for dirname in ['ancil', 'calib', 'calib-std', 'sci']:
		if not os.path.exists(f"{abs_ob_dir}/{dirname}"):
			print("Incorrect directory structure, please run ")
			print("kreduce.calibration.rename_raw_data() to initialize files for the pipeline")
			sys.exit()

	print(f"\n*** [kreduce-calibration]: running {task.upper()} calibrations for {abs_ob_dir}")

	### DARKS:
	if task == "all" or task == "dark":
		print(f"\n*** [kreduce-calibration]: processing DARKS")

		os.chdir(f"{abs_ob_dir}/calib") # Change to relevant directory
		
		dark_files = glob.glob(f"*-dark.fits") # Get the dark fits
		
		# Raw Frame Input (DARK)
		sof = open(f"dark.sof", 'w')
		for dark in dark_files:
			utils.add_file_to_sof(sof=sof, file=dark, file_type='DARK')
		sof.close()
		
		# Run kmos_dark recipe - Produces Product frames MASTER_DARK and BADPIXEL_DARK
		process = subprocess.run(["esorex", "kmos_dark", "dark.sof"], 
	       	stdout=subprocess.PIPE, 
	        universal_newlines=True)

		# Change the esorex log to be specific to eso-rex dark
		shutil.move(src="esorex.log", dst='esorex-dark.log')
		os.chdir(root_dir) # return to root directory

	### FLATS:
	if task == "all" or task == "flat":
		print(f"\n*** [kreduce-calibration]: processing FLATS")

		os.chdir(f"{abs_ob_dir}/calib")
		sof = open(f"flat.sof", "w")

		# Raw Frame Inputs (FLAT_ON, FLAT_OFF)
		flat_off_files = glob.glob("*-flat-off.fits")
		for flat in flat_off_files:
			utils.add_file_to_sof(sof=sof, file=flat, file_type="FLAT_OFF")

		flat_on_files = glob.glob("*-flat-lamp.fits")
		for flat in flat_on_files:
			utils.add_file_to_sof(sof=sof, file=flat, file_type='FLAT_ON')
		
		# from kmos_dark
		sof.write("BADPIXEL_DARK.fits\tBADPIXEL_DARK\n")
		sof.close()

		# Run kmos_flat recipe - Produces XCAL, YCAL, MASTER_FLAT, FLAT_EDGE, BADPIXEL_FLAT Product Frame
		process = subprocess.run(["esorex", "kmos_flat", "flat.sof"], 
	       	stdout=subprocess.PIPE, 
	        universal_newlines=True)

		# Change esorex.log to a specific flat log
		shutil.move(src="esorex.log", dst="esorex-flat.log")
		os.chdir(root_dir)

	### WAVELENGTH CALIBRATION:
	if task == "all" or task == "wave_cal":
		print(f"\n*** [kreduce-calibration]: processing WAVELENGTH CALIBRATION")

		os.chdir(f"{abs_ob_dir}/calib")
		sof = open("wave_cal.sof", "w")

		# Raw frame inputs (ARC_ON, ARC_OFF)
		wave_off_files = glob.glob('*-wave-off.fits')
		for wave in wave_off_files:
			utils.add_file_to_sof(sof=sof, file=wave, file_type="ARC_OFF")

		wave_on_files = glob.glob('*-wave-lamp.fits')
		for wave in wave_on_files:
			utils.add_file_to_sof(sof=sof, file=wave, file_type="ARC_ON")

		# From kmos_flat
		sof.write("XCAL_HKHKHK.fits\tXCAL\n") # X Calibrations
		sof.write("YCAL_HKHKHK.fits\tYCAL\n") # Y Calibrations
		sof.write("FLAT_EDGE_HKHKHK.fits\tFLAT_EDGE\n") # Flat Edge 

		# Static Calibration Frames (ARC_LIST, REF_LINES, WAVE_BAND)
		sof.write(f"{kreduce.kmos_calib_path}/kmos_ar_ne_list_hk.fits\tARC_LIST\n")
		sof.write(f"{kreduce.kmos_calib_path}/kmos_wave_ref_table.fits\tREF_LINES\n")
		sof.write(f"{kreduce.kmos_calib_path}/kmos_wave_band.fits\tWAVE_BAND\n")
		sof.close()

		# Run kmos_wave_cal recipe - Produces LCAL, DET_IMG_WAVE Product Frames
		process = subprocess.run(["esorex", "kmos_wave_cal", "wave_cal.sof"], 
	       	stdout=subprocess.PIPE, 
	        universal_newlines=True)

		shutil.move(src="esorex.log", dst="esorex-wave.log")
		os.chdir(root_dir)

	### ILLUMINATION CORRECTION:
	if task == "all" or task == "illum":
		print(f"\n*** [kreduce-calibration]: processing ILLUMINATION CORRECTION")

		os.chdir(f"{abs_ob_dir}/calib")
		sof = open("illum.sof", "w")

		# Raw Frame inputs (FLAT_ON)
		flat_on_files = glob.glob("*-flat-lamp.fits")
		for flat in flat_on_files:
			utils.add_file_to_sof(sof=sof, file=flat, file_type="FLAT_ON")

		# From kmos_flat, kmos_wave_cal
		sof.write("XCAL_HKHKHK.fits\tXCAL\n")
		sof.write("YCAL_HKHKHK.fits\tYCAL\n")
		sof.write("LCAL_HKHKHK.fits\tLCAL\n")
		sof.write("FLAT_EDGE_HKHKHK.fits\tFLAT_EDGE\n")

		# Static Calibration Frames (WAVE_BAND)
		sof.write(f"{kreduce.kmos_calib_path}/kmos_wave_band.fits\tWAVE_BAND\n")
		sof.close()

		# Run kmos_illumination recipe - Produces ILLUM_CORR
		process = subprocess.run(["esorex", "kmos_illumination", 
			"-pix_scale=0.1", "illum.sof"], 
	       	stdout=subprocess.PIPE, 
	        universal_newlines=True)

		shutil.move(src="esorex.log", dst="esorex-illum.log")
		os.chdir(root_dir)

	### TELLURIC:
	if task == "all" or task == "std_star":
		print(f"\n*** [kreduce-calibration]: processing STANDARD STAR")

		os.chdir(f"{abs_ob_dir}/calib-std")
		sof = open("std.sof", 'w')

		# Raw Frame Input (STD)
		std_files = glob.glob("*-object-sky-std-flux.fits")
		for std in std_files:
			utils.add_file_to_sof(sof=sof, file=std, file_type="STD")

		# from kmos_flat, kmos_illumination
		sof.write(f"{abs_ob_dir}/calib/XCAL_HKHKHK.fits\tXCAL\n")
		sof.write(f"{abs_ob_dir}/calib/YCAL_HKHKHK.fits\tYCAL\n")
		sof.write(f"{abs_ob_dir}/calib/LCAL_HKHKHK.fits\tLCAL\n")
		sof.write(f"{abs_ob_dir}/calib/MASTER_FLAT_HKHKHK.fits\tMASTER_FLAT\n")
		sof.write(f"{abs_ob_dir}/calib/ILLUM_CORR_HKHKHK.fits\tILLUM_CORR\n")

		# Uses Static Calibration Frames (WAVE_BAND, SOLAR_SPEC, ATMOS_MODEL, SPEC_TYPE_LOOKUP)
		sof.write(f"{kreduce.kmos_calib_path}/kmos_wave_band.fits\tWAVE_BAND\n")
		sof.write(f"{kreduce.kmos_calib_path}/kmos_solar_hk_1100.fits\tSOLAR_SPEC\n")
		sof.write(f"{kreduce.kmos_calib_path}/kmos_atmos_hk.fits\tATMOS_MODEL\n")
		sof.write(f"{kreduce.kmos_calib_path}/kmos_spec_type.fits\tSPEC_TYPE_LOOKUP\n")
		sof.close()

		# Run kmos_std_star recipe
		process = subprocess.run(["esorex", "kmos_std_star", "std.sof"], 
	       	stdout=subprocess.PIPE, 
	        universal_newlines=True)

		# Can produce TELLURIC, STAR_SPEC, STD_IMAGE, STD_MASK - unsure which one is necessarily used going forwarrd but will se.
		shutil.move(src="esorex.log", dst="esorex-std.log")
		os.chdir(root_dir)

	# All tasks completed
	print("\n*** [kreduce-calibration]: Done\n")
















