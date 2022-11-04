# imports
import glob
import shutil 
import subprocess

def add_file_to_sof(sof, file, file_type):
	""" Adds a file to the sof file """
	sof.write(f"{file}\t{file_type}\n")
	return True

def write_ob_p_file(dfile: str, destination: str, pointing: str, ob: int):
	with open(f'{destination}/obj_dict.txt', 'a') as file:
		file.write(f'dfile: {dfile} | p: {pointing} | ob: {ob}\n')

def move_reduced_ob_files(ob_dir, destination, pointing: str, ob: int, skycorr=True):
	""" Moves reduced ob files to a destination based on whether or not they've been put through the sky corrections """
	if skycorr is True:
		dfiles = glob.glob("{:}/sci/SCI_RECON*-sci_SKYCORR.fits".format(ob_dir))
		for dfile in dfiles:
			dfile_new = dfile.replace(dfile.split('.')[2], f'OB{ob:02}') #
			subprocess.call(f'mv {dfile} {dfile_new}', shell=True) #
			dfile = dfile_new #
			shutil.copy(src=dfile, dst=destination)
			write_ob_p_file(dfile, destination, pointing, ob)
	else:
		dfiles = glob.glob("{:}/sci/SCI_RECON*-sci.fits".format(ob_dir))
		for dfile in dfiles:
			dfile_new = dfile.replace(dfile.split('.')[2], f'OB{ob:02}') #
			subprocess.call(f'mv {dfile} {dfile_new}', shell=True) #
			dfile = dfile_new #
			shutil.copy(src=dfile, dst=destination)
			write_ob_p_file(dfile, destination, pointing, ob)



