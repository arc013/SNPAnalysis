import sys

def parseVCF(filename):
	info = {}
	samples = {}
	header = []
	colnames = []
	lastline = ''
	for line in open(filename):
		if '#' in line:
			colnames = line 
			header.append(''.join(lastline))
			lastline = line.split('\t')
		else:			
			sep_line = line.split('\t')
			ID = 

				
