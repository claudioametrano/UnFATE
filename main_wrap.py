#!/usr/bin/python3
'''
Main wrapper for running Hybpiper and subsequent analysis.


'''
import shutil
import re
import os
import sys
import argparse
import shlex
import multiprocessing
import itertools
import subprocess
import logging
from os import path

def run_exonerate_hits(file_, ref_seq_file):
	logging.info("Extracting genes from: " +file_)
	fline=open(file_).readline()
	regex_spades_header =re.search("^>NODE_[0-9]+_length_[0-9]+_cov_[0-9]+",fline)
	#if spades assembly, run exonerate_hits from HybPiper
	if regex_spades_header != None:
		os.system("python3 exonerate_hits.py {} --prefix {} {} ".format(ref_seq_file, file_.rstrip("\.fna"), file_))
	# else use the script version not using coverage information
	else:
		os.system("python3 exonerate_alt.py {} --prefix {} {} ".format(ref_seq_file, file_.rstrip("\.fna"), file_))
	return(file_)

def get_alignment(path_to_data):
	genes_list = []
	for root, dirs, files in os.walk(path_to_data, topdown=True):
		#append any gene to a list, make it a set to eliminate redundancy and then back to a list		
		for f in files:
			if f.endswith(".FAA"):
				genes_list.append(f.rstrip(".FAA"))
	genes_list = set(genes_list)
	genes_list = list(genes_list)
	print("List of genes found: ", genes_list)
	os.chdir(path_to_data)
	for g in genes_list:
		logging.info("Building nucleotide alignment for gene {}".format(g))
		with open("Alignment_" + g + "_nucleotide.FAS",'a+') as alignment:
			for root, dirs, files in os.walk(path_to_data, topdown=True):
				for f in files:
					if f == g + ".FNA":
						os.chdir(root)
						#print("root: " + root + f)
						with open(f, 'r') as gene:
							f_content = gene.read()
							os.chdir(path_to_data)
							#print("path to data:", path_to_data)
							alignment.write(f_content)	
	for g in genes_list:
		logging.info("Building aminoacid alignment for gene {}".format(g))
		with open("Alignment_" + g + "_protein.FAS",'a+') as alignment:
			for root, dirs, files in os.walk(path_to_data, topdown=True):
				for f in files:
					if f == g + ".FAA":
						os.chdir(root)
						print("root" + root + f)
						with open(f, 'r') as gene:
							f_content = gene.read()
							os.chdir(path_to_data)
							#print("path to data: ", path_to_data)
							alignment.write(f_content)
	return()

def merge_alignments(path_to_assemblies, path_to_target_enrichment):
	output_path = path_to_target_enrichment + '../'
	if not os.path.exists('merged_alignments'):
		os.makedirs(output_path + 'alignments')
		os.makedirs(output_path + 'alignments_merged')
	output_path = output_path + 'alignments/'
	output_path_merged = output_path.rstrip("/") + "_merged/"
	#print(output_path)
	for filename in (os.listdir(path_to_assemblies)):
		if filename.endswith(".FAS"):
			os.rename(path_to_assemblies + filename , path_to_assemblies + filename.rstrip("FAS") + "assembly.FA")
	for filename in (os.listdir(path_to_target_enrichment)):
		if filename.endswith(".FAS"):
			os.rename(path_to_target_enrichment + filename , path_to_target_enrichment + filename.rstrip("FAS") + "targenrich.FA")
	for filename in (os.listdir(path_to_assemblies)):
		if filename.endswith("assembly.FA"):
			shutil.move(path_to_assemblies + filename , output_path + filename) 
	for filename in (os.listdir(path_to_target_enrichment)):
		if filename.endswith("targenrich.FA"):
			shutil.move(path_to_target_enrichment + filename , output_path + filename)
	file_list_total = os.listdir(output_path)
	for index, item in  enumerate(file_list_total):
		file_list_total[index] = (file_list_total[index]).replace(".targenrich.FA", "")
		file_list_total[index] = (file_list_total[index]).replace(".assembly.FA", "")	
	#print(file_list_total)
	file_list_total = list(set(file_list_total))
	print(file_list_total)
	for item in file_list_total:
		with open(output_path_merged + item + "_merged.fasta", 'a+') as f:
			for filename in os.listdir(output_path):
				regex = re.search('(^Alignment_[0-9]+at4890_\w+)\.\w+\.FA$', filename )
				if item == regex.group(1):
					my_file = open(output_path + filename, 'r')
					my_file_content = my_file.read()
					f.write(my_file_content)
	return()

def run_gblocks(DNAextension, AAextension, path): 
	"""USE: Launch Gblocks getting the relaxed setting from the alignments characteristics (defaults are relaxed setting from Talavera & Castresana 2007)"""  
	for gene_file in os.listdir(path):
		if gene_file.endswith(DNAextension):
			fraction1=0.56
			fraction2=0.56
			b3=10
			b4=5
			print( 'File being processed: %s' %gene_file)
			count = 0
			my_file = open(gene_file, "r")
			my_file_content = my_file.readlines()
			for line in my_file_content:	
				if line.startswith(">"):
					count = count + 1
			print("The alignment has: ",count," sequences")
			b1 = int(count * fraction1)
			b2 = int(count * fraction2)	
			print("Number of char in a column of the alignment to be considered conserved and flanking regions, respectively: ", b1, b2)		
			start_Gblocks = "./Gblocks %s -t=c -b1=%s -b2=%s -b3=10 -b4=5 -b5=h -e=_Gb.fas" % (gene_file, str(b1), str(b2)) 		
			os.system(start_Gblocks) 
		elif gene_file.endswith(AAextension):
			fraction1=0.56
			fraction2=0.56
			b3=10
			b4=5
			print( 'File being processed: %s' %gene_file)
			count = 0
			my_file = open(gene_file, "r")
			my_file_content = my_file.readlines()
			for line in my_file_content:	
				if line.startswith(">"):
					count = count + 1
			print("The alignment has: ",count," sequences")
			b1 = int(count * fraction1)
			b2 = int(count * fraction2)	
			print("Number of char in a column of the alignment to be considered conserved and flanking regions, respectively: ", b1, b2)		
			start_Gblocks = "./Gblocks %s -t=p -b1=%s -b2=%s -b3=10 -b4=5 -b5=h -e=_Gb.fas" % (gene_file, str(b1), str(b2)) 		
			os.system(start_Gblocks) 

def check_arg(args=None):
	''' 
	Argument input for Wrapper: target enrichment fastq files, assemblies fasta, WGS fasta AND target genes fasta 
	Dependancies: mafft, RAxML, IQTREE, ASTRALIII in your $PATH
	
	'''
	parser = argparse.ArgumentParser(description='Run the whole pipeline to raw data to phylogenetic tree')
	parser.add_argument('-b', '--target_markers', default= '',
						help=' Path to fasta files containg all the sequences used to design the bait set, MUST BE A PROTEIN FASTA, USE  AN ABSOLUTE PATH!'
						)
	parser.add_argument('-ch', '--cpu', default= '4',
						help='CPU number used by Hybpiper or parallel run of Exonerate, MACSE, RAxML etc.' 
						)				
	parser.add_argument('-t', '--target_enrichment_data', default= '',
						help='Path to target enriched data, USE AN ABSOLUTE PATH!',
						)
	parser.add_argument('-w', '--whole_genome_data', default= '',
						help='Input path of de novo whole genome sequence data.'
						)
	parser.add_argument('-a', '--assemblies', default= '',
						help='Path to assemblies, USE AN ABSOLUTE PATH!',
						)	
	parser.add_argument('-f', '--first_use', action= 'store_true', 
						help='Clones Hybpiper and MACSE pipelines to your script directory from Github, use this argument only if is the first time you run the pipeline',
						)		
	parser.add_argument('-g', '--gblocks_relaxed', action= 'store_true', 
						help='applies Gblocks with relaxed parameters (Talavera et al. 2007)',
						)									
	return parser.parse_args(args)
args = check_arg(sys.argv[1:])


def main():
	#print(args)
	os.system("ulimit -n 1024000")
	main_script_dir = os.path.realpath(__file__)
	main_script_dir = main_script_dir.rstrip("main_wrap.py")
	#print(main_script_dir)
	#print(args.target_enrichment_data)
	#print(args.assemblies)
	#Clones hybpiper into current directory
	os.system('ulimit -n 1024000')
	if args.first_use == True:
		clone_hybpiper = 'git clone https://github.com/mossmatters/HybPiper.git'
		os.system(clone_hybpiper)
		logging.info("Hybpiper cloned")
		clone_MACSE = 'git clone https://github.com/ranwez/MACSE_V2_PIPELINES.git'
		os.system(clone_MACSE)
		logging.info("MACSE alignment pipeline cloned")
		# commands to modity OMM_MACSE main script and utilities with the right path to them instead of having to install Singularity to use it
		MACSE_dir = main_script_dir + "MACSE_V2_PIPELINES/OMM_MACSE"
		MACSE_script = MACSE_dir + "S_OMM_MACSE_V10.02.sh"
		MACSE_utils_dir = MACSE_dir + "UTILS"
		for filename in (os.listdir(MACSE_dir)):
			print(filename)
			if filename == "S_OMM_MACSE_V10.02.sh":
					fin = open(MACSE_dir + "/" + filename, "rt")
					data = fin.read()
					data = data.replace('LG_UTILS=${LG_UTILS_PATH}','LG_UTILS=' +MACSE_utils_dir)
					data = data.replace('mafft="${LG_MAFFT} --quiet $ALIGNER_EXTRA_OPTION"', 'mafft="mafft --quiet $ALIGNER_EXTRA_OPTION"')
					data = data.replace('hmmcleaner="perl ${LG_HMMCLEANER}"', 'hmmcleaner="perl ' + MACSE_utils_dir + '/HMMcleanerV1_8_VR2/HMMcleanAA_VR.pl"')
					data = data.replace('macse="java -jar -Xmx${JAVA_MEM} ${LG_MACSE}"', 'macse="java -jar -Xms4g -Xmx8g' + MACSE_script)
					fin.close()
					fin = open(MACSE_dir + "/" + filename, "wt")
					fin.write(data)
					fin.close()
			if filename == "S_fasta1L.sh":
				fin1 = open(MACSE_dir + "UTILS" + "/" + filename, "rt")
				data1 = fin1.read()
				data1 = data1.replace('LG_UTILS=${LG_UTILS_PATH}','LG_UTILS=' + MACSE_utils_dir)
				fin1.close()
				fin1 = open(MACSE_dir + "/" + filename, "wt")
				fin1.write(data1)
				fin1.close()
				
	path_to_sequences = args.target_enrichment_data
	if args.target_enrichment_data:
		logging.info("********************************************************************************************")
		logging.info("* PERFORMING TARGET ENRICHMENT DATA ANALYSIS WITH Hybpiper  *")
		logging.info("********************************************************************************************")
		logging.info('Path to TE data: '+path_to_sequences)
		trimming_cmd = "python3 {}/trimmer.py -f {}".format(main_script_dir, args.target_enrichment_data)
		os.system(trimming_cmd)
		#Get namelist.txt from dataset directory
		namelist_cmd = 'python3 {}/getNameList.py -f {}'.format(main_script_dir, args.target_enrichment_data)
		os.system(namelist_cmd)
		namelist = 'namelist.txt'
		path_to_namelist = os.path.join(path_to_sequences,namelist)
		
		logging.info("Gunzipping paired reads trimmed fastq archives")
		gunzip_fastq =' parallel gunzip ::: {}*_paired.fastq.gz'.format(path_to_sequences) 
		os.system(gunzip_fastq)
		logging.info("Running Hybpiper on target enrichment data provided")
		os.chdir(path_to_sequences)
		with open(path_to_namelist, 'r') as f:
			for line in f:
				logging.info("Processing sample:" + line)
				sample_path = path_to_sequences + '/' + line.rstrip('\n') + '_R*.trimmed_paired.fastq'
				run_Hybpiper =  '{}HybPiper/reads_first.py -b {} -r {}  --prefix {} --cpu {} '.format(main_script_dir, args.target_markers, sample_path, line, args.cpu)
				os.system(run_Hybpiper)
		os.chdir(main_script_dir)	
				
	"""if argument is whole genome input data:
			# SAME AS TARGET ENRICHMENT DATA OR MAYBE ASSEMBLY WITH SPADES THEN EXONERATE?? 
			TEST THIS!!"""
			
	#user input: assemblies
	if args.assemblies:
		path_to_assemblies = args.assemblies
		logging.info("*********************************************************************************")
		logging.info("* PERFORMING ASSEMBLIES DATA ANALYSIS WITH Exonerate  *")
		logging.info("*********************************************************************************")
		logging.info('Path to assemblies '+path_to_assemblies)
		#logging.info('Created hybpiper directory in assembly sequence directory')
		#os.chdir(path_to_assemblies)
		
		for root, dirs, files in os.walk(path_to_assemblies, topdown=True):
			for name in files:
				if name.endswith(".fna.gz"):
					os.system("gunzip "+ path_to_assemblies + name)		
		pezizo_list = []	
		for root, dirs, files in os.walk(path_to_assemblies, topdown=True):
			for name in files:
				if name.endswith(".fna"):
					pezizo_list.append(root + name)
		#print("Samples are: ", pezizo_list)
		empty_list = []
		empty_list.append(args.target_markers)
		#print(empty_list)
		# product function from itertools does the cartesian product (lane * rows), it is like a nested for!!
		list_of_list = list(itertools.product(pezizo_list, empty_list))
		print(list_of_list)
		logging.info("Running exonerate using exonerate_hits.py script from Hybpiper..")	
		args.cpu = int(args.cpu)
		pool = multiprocessing.Pool(processes=args.cpu)
		pool.starmap(run_exonerate_hits, list_of_list)
		
		logging.info("*********************************")
		logging.info("* BUILDING FASTA FILES  *")
		logging.info("*********************************")
		logging.info("Building alignments from assemblies data")
		get_alignment(args.assemblies)
		logging.info("Building alignments from target enrichment data")
		get_alignment(args.target_enrichment_data)
		merge_alignments(args.assemblies, args.target_enrichment_data)
		
		logging.info("***************************************************************")
		logging.info("* PERFORMING ALIGNMENT WITH OMM_MACSE *")
		logging.info("***************************************************************")
		logging.info("**************************************************")
		logging.info(" (mafft version) with HMMcleaner filtering  *")
		logging.info("**************************************************")
		path_to_merged_alignments = args.target_enrichment_data + '../alignments_merged/'
		# maybe better a MACSE2 pipeline using mafft
		#run_mafft_parallel = "find %s -type f -name '*_merged.fasta' | parallel -j %s mafft --maxiterate 1000 --localpair --thread 1 {} > {}_maffted.fas    "(path_to_merged_alignments, args.parallel_mafft, )
		#os.system(run_mafft_parallel)
		MACSE_dir = main_script_dir + "MACSE_V2_PIPELINES/OMM_MACSE/"
		MACSE_script = MACSE_dir + "S_OMM_MACSE_V10.02.sh"
		run_OMM_MACSE = "find %s -type f -name '*_merged.fasta' | parallel -j %s %s --out_dir {}_out --out_file_prefix {} --in_seq_file {} --no_prefiltering --nopostfiltering --aligner_extra_option '--localpair --maxiterate 1000' --min_percent_NT_at_ends 0.5 " %(path_to_merged_alignments, args.cpu, MACSE_script)
		os.system(run_OMM_MACSE)
		logging.info("******************************************************************************************")
		logging.info("* PERFORMING ALIGNMENT FILTERING WITH Gblocks (relaxed param.)  *")
		logging.info("******************************************************************************************")
		make_align_fold = "mkdir {}macsed_alignments".format(path_to_merged_alignments + "../")
		os.system("make_align_fold")
		copy_MACSE_output = "find %s -type f -name '*mask_align_NT.aln | cp -r {} %s '"%(path_to_merged_alignment, path_to_merged_alignments + "../macsed_alignments")
		copy_MACSE_output = "find %s -type f -name '*mask_align_AA.aln | cp -r {} %s '"%(path_to_merged_alignment, path_to_merged_alignments + "../macsed_alignments")
		os.system(copy_MACSE_output)
		os.system(copy_MACSE_output1)
		if args.gblocks_relaxed:
			run_gblocks("_final_mask_align_NT.aln","macsed_final_mask_align_AA.aln", path_to_merged_alignments + "../macsed_alignments")
			os.system(run_gblocks)	
		logging.info("************************************************************************************")
		logging.info("* RECONSTRUCTING SINGLE MARKER TREES WITH RAxML-NG  *")
		logging.info("************************************************************************************")
		#RAxml_parallel = "find %s -type f  -name 'mask_align_NT.aln' | parallel -j %s raxmlHPC-PTHREADS-AVX -s {} -m GTRGAMMA -fa -p %s -x %s -#1000 -n raxml -T 1".format(path_to_merged_alignment + "../macsed_alignments", args.cpu, rand1, rand2 )
		#RAxml_parallel = "find %s -type f  -name 'mask_align_NT.aln_Gb.fas' | parallel -j %s raxmlHPC-PTHREADS-AVX -s {} -m GTRGAMMA -fa -p %s -x %s -#1000 -n raxml -T 1".format(path_to_merged_alignment + "../macsed_alignments", args.cpu, rand1, rand2 )
		raxml_script = main_script_dir
		#RAxML-NG instead of RAxML
		RAxml_parallel = "find %s -type f  -name 'mask_align_NT.aln' | parallel -j %s ./raxml-ng --all --msa {} --model GTR+G --prefix raxml --seed 888 --threads 1 --bs-metric tbe --tree pars{20},rand{20} ".format(path_to_merged_alignment + "../macsed_alignments", args.cpu)
		RAxml_parallel = "find %s -type f  -name 'mask_align_NT.aln_Gb.fas' | parallel -j %s ./raxml-ng --all --msa {} --model GTR+G --prefix raxml --seed 888 --threads 1 --bs-metric tbe --tree pars{20},rand{20} ".format(path_to_merged_alignment + "../macsed_alignments", args.cpu)
		#amino acid alignments
		RAxml_parallel = "find %s -type f  -name 'mask_align_AA.aln' | parallel -j %s ./raxml-ng --all --msa {} --model PROTGAMMAAUTO --prefix raxml --seed 888 --threads 1 --bs-metric tbe --tree pars{20},rand{20} ".format(path_to_merged_alignment + "../macsed_alignments", args.cpu)
		RAxml_parallel = "find %s -type f  -name 'mask_align_AA.aln_Gb.fas' | parallel -j %s ./raxml-ng --all --msa {} --model PROTGAMMAAUTO --prefix raxml --seed 888 --threads 1 --bs-metric tbe --tree pars{20},rand{20} ".format(path_to_merged_alignment + "../macsed_alignments", args.cpu)

		logging.info("************************************************************************************")
		logging.info("* PERFORMING ALIGNMENTS CONCATENATION WITH Fasconcat  *")
		logging.info("************************************************************************************")
		#Fasconcat = " 
		logging.info("*******************************************************************************")
		logging.info("* RECONSTRUCTING SUPERMATRIX TREE WITH RAxML-NG  *")
		logging.info("*******************************************************************************")
		
		logging.info("******************************************************************")
		logging.info("* RECONSTRUCTING SUPERTREE WITH ASTRAL  *")
		logging.info("******************************************************************")
if __name__=='__main__':
	logger = logging.getLogger(__name__)
	logFormatter = '%(message)s'
	logging.basicConfig(filename='pipe_logger.log', format=logFormatter, level=logging.DEBUG)
	handler = logging.StreamHandler(sys.stdout)
	main()


