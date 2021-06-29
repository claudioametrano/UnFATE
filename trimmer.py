import os 
import re
import logging
import click

@click.command()
@click.option('--data_folder','-f', default='./', help='path of the folder containing the fastq.gz file to trim with trimmomatic')
@click.option('--cpus', '-c', default=20, help="number of cpus for trimmomatic to use")

def trimming(data_folder, cpus):
	logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	script_path = os.path.realpath(__file__)
	os.system("ulimit -n 1024000")
	trimmomatic_path = script_path.rstrip('trimmer.py') + '/' + 'Trimmomatic-0.39'
	samplesR1 = []
	samplesR2 = []
	samplesSE = []
	for root, dirs, files in os.walk(data_folder, topdown=True):
		for name in files:
			if name.endswith("R1.fastq.gz"):
				filepath=os.path.join(root, name)
				samplesR1.append(filepath)
		for name in files:
			if name.endswith("R2.fastq.gz"):
				filepath=os.path.join(root, name)
				samplesR2.append(filepath)
		for name in files:
			if name.endswith("SE.fastq.gz"): 
				filepath=os.path.join(root, name)
				samplesSE.append(filepath)				
	#print(samplesR1)			
	#print(samplesR2)
	#print(samplesSE)
	for path1 in samplesR1:
		regex1 = re.search('(.+)R1.fastq.gz$',path1)
		path1out = path1.rstrip("fastq.gz")
		for path2 in samplesR2:
			regex2 = re.search('(.+)R2.fastq.gz$',path2)
			path2out = path2.rstrip("fastq.gz")
			if regex1.group(1) == regex2.group(1):
				logging.info("Trimming sample: " + regex1.group(1) + 'R*.fastq.gz' )
				TRIMMOMATIC_COMMAND = "java -jar {}/trimmomatic-0.39.jar PE -threads {} -phred33 {} {} {}.trimmed_paired.fastq.gz {}.trimmed_UNpaired.fastq.gz {}.trimmed_paired.fastq.gz {}.trimmed_UNpaired.fastq.gz ILLUMINACLIP:{}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:25".format(trimmomatic_path, cpus, path1, path2, path1out, path1out, path2out, path2out, trimmomatic_path)
				os.system(TRIMMOMATIC_COMMAND)	
	for path3 in samplesSE:
		path3out = path3.rstrip("fastq.gz")
		logging.info("Trimming sample: " + path3)	
		TRIMMOMATIC_COMMAND = "java -jar {}/trimmomatic-0.39.jar SE -threads {} -phred33 {} {}.trimmed.fastq.gz  ILLUMINACLIP:{}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:25".format(trimmomatic_path, cpus, path3, path3out, trimmomatic_path)
		os.system(TRIMMOMATIC_COMMAND)					
			
# starts the function					
if __name__ == '__main__':
	trimming()
