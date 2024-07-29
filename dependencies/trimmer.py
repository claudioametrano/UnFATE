import os
import re
import logging
import click
from glob import glob

@click.command()
@click.option('--data_folder','-f', default='./', help='path of the folder containing the fastq.gz file to trim with trimmomatic')
@click.option('--cpus', '-c', default=20, help="number of cpus for trimmomatic to use")

def trimming(data_folder, cpus):
	logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	script_path = os.path.realpath(__file__)
	os.system("ulimit -n 1024000")
	trimmomatic_path = script_path.rstrip('trimmer_dev_mod.py') + '/' + 'Trimmomatic-0.39'
	samplesR1 = []
	samplesR2 = []
	samplesSE = []

	for name in glob(os.path.join(data_folder, "*R1.fastq*")):
		samplesR1.append(name)
	for name in glob(os.path.join(data_folder, "*R2.fastq*")):
		samplesR2.append(name)
	for name in glob(os.path.join(data_folder, "*SE.fastq*")):
		samplesSE.append(name)

	#print(samplesR1)
	#print(samplesR2)
	#print(samplesSE)
	for path1 in samplesR1:
		regex1 = re.search('(.+)R1.fastq*',path1)
		path1out = path1.split(".")[0]
		#this assumes that there won't be a period before .fastq(.gz)
		for path2 in samplesR2:
			regex2 = re.search('(.+)R2.fastq*',path2)
			path2out = path2.split(".")[0]
			if regex1.group(1) == regex2.group(1):
				logging.info("Trimming sample: " + regex1.group(1) + 'R*.fastq*' )
				TRIMMOMATIC_COMMAND = "java -jar {}/trimmomatic-0.39.jar PE -threads {} -phred33 {} {} {}.trimmed_paired.fastq.gz {}.trimmed_UNpaired.fastq.gz {}.trimmed_paired.fastq.gz {}.trimmed_UNpaired.fastq.gz ILLUMINACLIP:{}/adapters/all_adapters.fas:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:25".format(trimmomatic_path, cpus, path1, path2, path1out, path1out, path2out, path2out, trimmomatic_path)
				os.system(TRIMMOMATIC_COMMAND)
	for path3 in samplesSE:
		path3out = path3.split(".")[0]
		logging.info("Trimming sample: " + path3)
		TRIMMOMATIC_COMMAND = "java -jar {}/trimmomatic-0.39.jar SE -threads {} -phred33 {} {}.trimmed.fastq.gz  ILLUMINACLIP:{}/adapters/all_adapters.fas:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:25".format(trimmomatic_path, cpus, path3, path3out, trimmomatic_path)
		os.system(TRIMMOMATIC_COMMAND)

# starts the function
if __name__ == '__main__':
	trimming()
