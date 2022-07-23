import click
import os 
import re
from os import path
import subprocess

@click.command()
@click.option('--accession_file','-f', default='accession.txt', help='accession numbers file, one accession per line')
@click.option('--out_file','-f', default='accession_plus_taxonomy.txt', help='output file')
def get_taxonomy(accession_file, out_file):
	""" Uses NCBI e-direct (esearch and other) and bash commands (awk grep etc) to get taxonomy from NCBI using the accession number and then to parse the output"""
	taxonomy_file = open(accession_file, "r")
	taxonomy_file_content = taxonomy_file.readlines()
	#if path.exists(out_file):
	#	print ("OUTPUT FILE already exists!")
	#else:	
	output_file = open(out_file ,"a+")
	for accession in taxonomy_file_content:
		print(" ")
		print("******************************")
		print("analyzying: ",accession)
		print("Searching the NCBI database...")
		output_file.write(accession.rstrip("\n"))
		output_file.write(",")
		NCBI_search = 'esearch -db assembly -query "%s" | elink -target taxonomy | efetch -format native -mode xml | grep ScientificName >temporary_output.txt' %(accession)
		os.system(NCBI_search)
		temp_out_file = open("temporary_output.txt", "r")
		temp_out_file_content = temp_out_file.readlines()
		for l in temp_out_file_content:
			#<ScientificName>Bacteria</ScientificName> output line look like this
			regex_pattern = re.search("<ScientificName>(.+?)</ScientificName>", l)
			print(regex_pattern.group(1))
			output_file.write(regex_pattern.group(1))
			output_file.write(",")
		output_file.write("\n")	
	output_file.close()		

# starts the function					
if __name__ == '__main__':
	get_taxonomy()
