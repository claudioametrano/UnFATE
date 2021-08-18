import os
import re
from os import listdir
from os.path import isfile, join

import click

@click.command()
@click.option('--data_folder','-f', default='./', help='path of the folder containing the trimmed fastq.gz files')

def sample_list(data_folder):
	#script_path = os.path.realpath(__file__)

	#Makes list of file names in current directory, set directory to input data
	onlyfiles = [f for f in listdir(data_folder) if isfile(join(data_folder, f))]
	#print(onlyfiles)
	namelist = []
	#Goes through list and finds the prefix "Name" of the file, must be unique sample name
	#Adds to a list of names
	for i in onlyfiles:
		regex = re.search('(.+)_R[1|2].trimmed_paired.fastq.gz', i)
		if regex != None:
			namelist.append(regex.group(1))
	#print(namelist)

	for i in onlyfiles:
		regex = re.search('(.+)_SE.trimmed.fastq.gz', i)
		if regex != None:
			namelist.append(regex.group(1))
	# make a set from the list and then back to list to get rid of double names
	nodouble_list = list(set(namelist))
	sorted_namelist = sorted(nodouble_list)
	#print(sorted_namelist)
	with open(data_folder + "namelist.txt",'w+') as f:
		f.write('\n'.join(sorted_namelist))

# starts the function					
if __name__ == '__main__':
	sample_list()

