import sys
import glob 
import os
import subprocess

def getVCFfiles(directory):
	filenames = []
	for root, dirs, files in os.walk(directory):
		for file in files:
			if file.endswith(".vcf"):
				filenames.append(os.path.join(root,file))	
	return filenames

def convertformat(cross_loc,convert_file, ref_file, directory):
	filenames = getVCFfiles(directory)
	if not os.path.exists('converted'):
		os.makedirs('converted')
	for file in filenames:
		filename = file.split('/')[-1]
		output_name = 'converted/'+filename
		header_process = subprocess.Popen(["python",cross_loc,'vcf',convert_file,ref_file,file,output_name], stdout=subprocess.PIPE)
		output, error = header_process.communicate()
		break

if __name__ == "__main__":
	cross_loc = sys.argv[1]
	convert_file = sys.argv[2]
	ref_file = sys.argv[3]
	directory = sys.argv[4]

	convertformat(cross_loc,convert_file, ref_file, directory)















