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

def run_gblocks(DNAextension, AAextension, path, path_to_gblocks): 
	"""USE: Launch Gblocks getting the relaxed setting from the alignments characteristics (defaults are relaxed setting from Talavera & Castresana 2007)"""  
	for gene_file in os.listdir(path):
		if gene_file.endswith(DNAextension):
			fraction1=0.56
			fraction2=0.56
			b3=10
			b4=5
			print( 'File being processed: %s' %gene_file)
			count = 0
			my_file = open(path +"/"+ gene_file, "r")
			my_file_content = my_file.readlines()
			for line in my_file_content:	
				if line.startswith(">"):
					count = count + 1
			print("The alignment has: ",count," sequences")
			b1 = str(int(count * fraction1))
			b2 = str(int(count * fraction2))	
			print("Number of char in a column of the alignment to be considered conserved and flanking regions, respectively: ", b1, b2)		
			start_Gblocks = "{} {} -t=c -b1={} -b2={} -b3=10 -b4=5 -b5=h -e=-gb".format(path_to_gblocks, path + gene_file, b1, b2)
			print(start_Gblocks)	
			os.system(start_Gblocks) 
		elif gene_file.endswith(AAextension):
			fraction1=0.56
			fraction2=0.56
			b3=10
			b4=5
			print( 'File being processed: %s' %gene_file)
			count = 0
			my_file = open(path +"/"+ gene_file, "r")
			my_file_content = my_file.readlines()
			for line in my_file_content:	
				if line.startswith(">"):
					count = count + 1
			print("The alignment has: ",count," sequences")
			b1 = str(int(count * fraction1))
			b2 = str(int(count * fraction2))	
			print("Number of char in a column of the alignment to be considered conserved and flanking regions, respectively: ", b1, b2)		
			start_Gblocks = "{} {} -t=p -b1={} -b2={} -b3=10 -b4=5 -b5=h -e=-gb".format(path_to_gblocks, path + gene_file, b1, b2) 		
			print(start_Gblocks)
			os.system(start_Gblocks) 

def check_arg(args=None):
	''' 
	Argument input for Wrapper: target enrichment fastq files, assemblies fasta, WGS fasta AND target genes fasta 
	
	'''
	parser = argparse.ArgumentParser(description='Run the whole pipeline to raw data to phylogenetic tree')
	parser.add_argument('-b', '--target_markers', default= '',
						help=' Path to fasta files containg all the sequences used to design the bait set, IT MUST BE A PROTEIN FASTA, USE  AN ABSOLUTE PATH!'
						)
	parser.add_argument('-c', '--cpu', default= '4',
						help='CPU number used by Hybpiper or parallel run of Exonerate, MACSE, RAxML etc.' 
						)				
	parser.add_argument('-t', '--target_enrichment_data', default= '',
						help='Path to target enriched data, USE AN ABSOLUTE PATH (include "/" after the folder name)! Folder must be named "target_enrichment". Files must end with "R1.fastq.gz" and "R2.fastq.gz"',
						)
#	parser.add_argument('-w', '--whole_genome_data', default= '',
#						help='Input path of de novo whole genome sequence data.'
#						)
	parser.add_argument('-a', '--assemblies', default= '',
						help='Path to assemblies, USE AN ABSOLUTE PATH (include "/" after the folder name)!  Folder must be named "assemblies", Files must end with ".fna.gz" or ".fna"',
						)	
	parser.add_argument('-f', '--first_use', action= 'store_true', 
						help='Clones Hybpiper and MACSE pipelines to your script directory from Github, use this argument only if is the first time you run the pipeline',
						)		
	parser.add_argument('-g', '--gblocks_relaxed', action= 'store_true', 
						help='Applies Gblocks with relaxed parameters (Talavera et al. 2007), both Gblock filtered and unfiltered data will be used to build phylogenies',
						)	
	parser.add_argument('-n', '--ncbi_assemblies', nargs = '+', 
						help='Extracts the pre-mined NCBI assemblies genes, then it uses the selected taxonomic rank to only include the needed samples. Can take a list of ranks (e.g. Morchella Tuber Fuffaria)'
						)	
	parser.add_argument('--nargs', nargs='+')																				
	return parser.parse_args(args)
args = check_arg(sys.argv[1:])


def main():
	#print(args)
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
		MACSE_dir = main_script_dir + "MACSE_V2_PIPELINES/OMM_MACSE/"
		MACSE_script = MACSE_dir + "S_OMM_MACSE_V10.02.sh"
		MACSE_utils_dir = main_script_dir + "MACSE_V2_PIPELINES/UTILS"
		for filename in (os.listdir(MACSE_dir)):
			if filename == "S_OMM_MACSE_V10.02.sh":
					fin = open(MACSE_dir + "/" + filename, "rt")
					data = fin.read()
					data = data.replace('LG_UTILS=${LG_UTILS_PATH}','LG_UTILS=' + MACSE_utils_dir)
					data = data.replace('mafft="${LG_MAFFT} --quiet $ALIGNER_EXTRA_OPTION"', 'mafft="mafft --quiet $ALIGNER_EXTRA_OPTION"')
					data = data.replace('muscle="${LG_MUSCLE} $ALIGNER_EXTRA_OPTION"', 'muscle="muscle $ALIGNER_EXTRA_OPTION"')
					data = data.replace('prank="${LG_PRANK} $ALIGNER_EXTRA_OPTION"', 'prank="prank $ALIGNER_EXTRA_OPTION"')
					data = data.replace('hmmcleaner="perl ${LG_HMMCLEANER}"', 'hmmcleaner="perl ' + MACSE_utils_dir + '/HMMcleanerV1_8_VR2/HMMcleanAA_VR.pl"')
					data = data.replace('macse="java -jar -Xmx${JAVA_MEM} ${LG_MACSE}"', 'macse="java -jar -Xms4g -Xmx8g ' + MACSE_utils_dir + '/macse_v2.03.jar"')
					fin.close()
					fin = open(MACSE_dir + "/" + filename, "wt")
					fin.write(data)
					fin.close()
		for filename in (os.listdir(MACSE_utils_dir + "/LGS_Fasta")):
			if filename == "S_fasta1L.sh":
				fin1 = open(MACSE_utils_dir + "/LGS_Fasta" + "/" + filename, "rt")
				data1 = fin1.read()
				data1 = data1.replace('LG_UTILS=${LG_UTILS_PATH}','LG_UTILS=' + MACSE_utils_dir)
				fin1.close()
				fin1 = open(MACSE_utils_dir + "/LGS_Fasta" + "/" + filename, "wt")
				fin1.write(data1)
				fin1.close()		
		clone_astral = "git clone https://github.com/smirarab/ASTRAL.git"
		os.system(clone_astral)
		os.chdir(main_script_dir + "ASTRAL/")
		os.system("./make.sh")
		os.chdir(main_script_dir)
		logging.info("ASTRAL cloned")
	path_to_sequences = args.target_enrichment_data
	
	if args.ncbi_assemblies:
		logging.info("Extracting pre-mined NCBI assemblies genes")
		path_to_premined = main_script_dir + "pre_mined_assemblies.tar.gz"
		unzip_premined_assemblies = "tar -zxf {}".format(path_to_premined)
		# ~ os.system(unzip_premined_assemblies)
		path_to_premined = main_script_dir + "pre_mined_assemblies/"
		path_to_taxonomy = path_to_premined + "Accession_plus_taxonomy_Pezizomycotina.txt"
		for item in args.ncbi_assemblies:
			# using the commas that are present in the csv file containing the taxonomy prevents to include ranks with similar name to the one requested in the command line (e.g. Fuffaria and Fuffarialongis)
			item1 = "," + item + ","
			with open(path_to_taxonomy.replace('Accession_plus_taxonomy_Pezizomycotina.txt','Accession_plus_taxonomy_reduced.txt'), 'a') as requested_rank:
				with open(path_to_taxonomy, 'r') as taxonomy:
					for line in taxonomy:
						if item1 in line:
							requested_rank.write(line)
						
				
	"""if argument is whole genome input data:
			# SAME AS TARGET ENRICHMENT DATA OR MAYBE ASSEMBLY WITH SPADES THEN EXONERATE?? 
			TEST THIS!!"""
			
	#user input: assemblies
	if args.assemblies:
		path_to_assemblies = args.assemblies
		logging.info("***********************************************************************")
		logging.info("PERFORMING ASSEMBLIES DATA ANALYSIS WITH Exonerate")
		logging.info("***********************************************************************")
		logging.info('Path to assemblies '+path_to_assemblies)
		
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
		#print(list_of_list)
		logging.info("Running exonerate using exonerate_hits.py script from Hybpiper..")	
		args.cpu = int(args.cpu)
		pool = multiprocessing.Pool(processes=args.cpu)
		pool.starmap(run_exonerate_hits, list_of_list)
	
	if args.target_enrichment_data:
		logging.info("*********************************************************************************")
		logging.info("TRIMMING TARGET ENRICHMENT FASTQ FILES  WITH TRIMMOMATC")
		logging.info("*********************************************************************************")
		# ~ logging.info('Path to TE data: '+path_to_sequences)
		# ~ trimming_cmd = "python3 {}/trimmer.py -f {}".format(main_script_dir, args.target_enrichment_data)
		# ~ os.system(trimming_cmd)
		# ~ #Get namelist.txt from dataset directory
		# ~ namelist_cmd = 'python3 {}/getNameList.py -f {}'.format(main_script_dir, args.target_enrichment_data)
		# ~ os.system(namelist_cmd)
		# ~ namelist = 'namelist.txt'
		# ~ path_to_namelist = os.path.join(path_to_sequences,namelist)
		# ~ logging.info("Gunzipping paired reads trimmed fastq archives")
		# ~ gunzip_fastq =' parallel gunzip ::: {}*_paired.fastq.gz'.format(path_to_sequences) 
		# ~ os.system(gunzip_fastq)
		# ~ logging.info("**************************************************************************************")
		# ~ logging.info("EXTRACTING GENES FROM TARGET ENRICHMENT DATA  WITH Hybpiper")
		# ~ logging.info("**************************************************************************************")
		# ~ os.chdir(path_to_sequences)
		# ~ with open(path_to_namelist, 'r') as f:
			# ~ for line in f:
				# ~ logging.info("Processing sample:" + line)
				# ~ sample_path = path_to_sequences + '/' + line.rstrip('\n') + '_R*.trimmed_paired.fastq'
				# ~ run_Hybpiper =  '{}HybPiper/reads_first.py -b {} -r {}  --prefix {} --cpu {} '.format(main_script_dir, args.target_markers, sample_path, line, args.cpu)
				# ~ os.system(run_Hybpiper)
		os.chdir(main_script_dir)		
		
		logging.info("****************************")
		logging.info("BUILDING FASTA FILES")
		logging.info("****************************")
		# ~ logging.info("Building alignments from assemblies data")
		# ~ get_alignment(args.assemblies)
		# ~ logging.info("Building alignments from target enrichment data")
		# ~ get_alignment(args.target_enrichment_data)
		# ~ merge_alignments(args.assemblies, args.target_enrichment_data)
		if args.ncbi_assemblies:
			logging.info("************************************************************************************************************")
			logging.info("ADDING SELECTED TAXONOMIC RANKS GENES FROM PRE-MINED ASSEMBLY DATABASE")
			logging.info("************************************************************************************************************")
			path_to_merged_alignments = args.target_enrichment_data.replace('target_enrichment/', 'alignments_merged/')
			#print(path_to_merged_alignments)
			path_to_ranks = path_to_taxonomy.replace('Accession_plus_taxonomy_Pezizomycotina.txt','Accession_plus_taxonomy_reduced.txt')
			#print(path_to_ranks)
			# make the folder for the selected sample from NCBI database
			path_to_premined = main_script_dir + "pre_mined_assemblies/"
			path_to_premined_selected = args.target_enrichment_data.replace('target_enrichment/', 'ncbi_selected_genes/')
			mkdir_ncbi_data ="mkdir {}".format(path_to_premined_selected)
			os.system(mkdir_ncbi_data)
			# make a list of the accessions selected according to the taxonomic rank chosen
			accession_from_ncbi_list = []
			with open(path_to_ranks, 'r') as ranks:
				for line in ranks:
					regex_accession = re.search('(GCA_[0-9]+.[0-9]),\w+',line)
					if regex_accession != None:
						accession_from_ncbi_list.append(regex_accession.group(1))
			# copy the gene folders from the "pre-mined assembly" database to the data folder
			for i in accession_from_ncbi_list:
				for root, dirs, files in os.walk(path_to_premined, topdown=True):
					for d in dirs:
						if i in d:
							copy_folder = "cp -r {} {}".format(path_to_premined + "/" + d, path_to_premined_selected )
							os.system(copy_folder) 
					
			# add the genes from the database	to the alignments		
			for fi in os.listdir(path_to_merged_alignments):
				if fi.endswith("_protein_merged.fasta"):
					regex_pattern = re.search("Alignment_([0-9]+at[0-9]+)_(\w+)_merged.fasta", fi)
					if regex_pattern != None:
						for root, dirs, files in os.walk(path_to_premined_selected, topdown=True):
							for f in files:
								if f.endswith(".FAA"):
									if regex_pattern.group(1) in f:
										with open(path_to_merged_alignments + fi,'a') as merged_ali:
											with open(root +"/"+ f, 'r') as gene_file:
												gene_file_content = gene_file.read()
												#print(gene_file_content)
												merged_ali.write(gene_file_content)
			for fi in os.listdir(path_to_merged_alignments):
				if fi.endswith("_nucleotide_merged.fasta"):
					regex_pattern = re.search("Alignment_([0-9]+at[0-9]+)_(\w+)_merged.fasta", fi)
					if regex_pattern != None:
						for root, dirs, files in os.walk(path_to_premined_selected, topdown=True):
							for f in files:
								if f.endswith(".FNA"):
									if regex_pattern.group(1) in f:
										with open(path_to_merged_alignments + fi,'a') as merged_ali:
											with open(root +"/"+ f, 'r') as gene_file:
												gene_file_content = gene_file.read()
												#print(gene_file_content)
												merged_ali.write(gene_file_content)


		logging.info("********************************************************")
		logging.info("PERFORMING ALIGNMENT WITH OMM_MACSE")
		logging.info("********************************************************")
		logging.info("*******************************************")
		logging.info("(mafft version) with HMMcleaner filtering")
		logging.info("*******************************************")
		# ~ path_to_merged_alignments = args.target_enrichment_data.replace('target_enrichment/', 'alignments_merged/')
		# ~ MACSE_dir = main_script_dir + "MACSE_V2_PIPELINES/OMM_MACSE/"
		# ~ MACSE_script = MACSE_dir + "S_OMM_MACSE_V10.02.sh"
		# ~ run_OMM_MACSE = 'find %s -type f -name "*_nucleotide_merged.fasta" | parallel -j %s %s --out_dir {}_out --out_file_prefix macsed --in_seq_file {} --no_prefiltering --no_postfiltering --alignAA_soft MAFFT  --min_percent_NT_at_ends 0.01 ' %(path_to_merged_alignments, args.cpu, MACSE_script)
		# ~ os.system(run_OMM_MACSE)
		# ~ print(run_OMM_MACSE)

		path_to_macsed_align = path_to_merged_alignments.replace('alignments_merged/','macsed_alignments/')
		# ~ make_align_fold = "mkdir {}".format(path_to_macsed_align)
		# ~ os.system(make_align_fold)
		# ~ for root, dirs, files in os.walk(path_to_merged_alignments, topdown=True):
			# ~ for f in files:
				# ~ if f.endswith("_final_align_NT.aln") or f.endswith("_final_align_AA.aln"):
					# ~ #print(root)
					# ~ #print(f)
					# ~ file_path = root +"/"+ f   
					# ~ regex1 =re.search("Alignment_([0-9]+at[0-9]+)_nucleotide_merged.fasta_out",root)
					# ~ #os.rename both renames and moves files
					# ~ os.rename(file_path,  path_to_macsed_align + regex1.group(1) + f + ".fas")
		
		gblocks_path= main_script_dir + "Gblocks"
		if args.gblocks_relaxed:
			logging.info("********************************************************************************")
			logging.info("PERFORMING ALIGNMENT FILTERING WITH Gblocks (relaxed param.)")
			logging.info("********************************************************************************")
			#print(path_to_macsed_align)
			# ~ run_gblocks("_final_align_NT.aln.fas","_final_align_AA.aln.fas", path_to_macsed_align, gblocks_path)	
			
		logging.info("**************************************************************************")
		logging.info("RECONSTRUCTING SINGLE MARKER TREES WITH RAxML-NG")
		logging.info("**************************************************************************")
		# ~ raxml_script = main_script_dir + "raxml-ng"
		# ~ print(path_to_macsed_align)
		# ~ raxml_parallel = "find %s -type f  -name '*_final_align_NT.aln.fas' | parallel -j %s %s --all --msa {} --model GTR+G --prefix {} --seed 888 --threads 1 --bs-metric tbe " %(path_to_macsed_align, args.cpu, raxml_script)
		# ~ os.system(raxml_parallel)
		# ~ if args.gblocks_relaxed:
			# ~ raxml_parallel1 = "find %s -type f  -name '*_final_align_NT.aln.fas-gb' | parallel -j %s %s --all --msa {} --model GTR+G --prefix {} --seed 888 --threads 1 --bs-metric tbe" %(path_to_macsed_align, args.cpu, raxml_script)
			# ~ os.system(raxml_parallel1)
		# ~ #amino acid alignments
		# ~ raxml_parallel2 = "find %s -type f  -name '*_final_align_AA.aln.fas' | parallel -j %s %s --all --msa {} --model PROTGTR+G --prefix {} --seed 888 --threads 1 --bs-metric tbe" %(path_to_macsed_align, args.cpu, raxml_script)
		# ~ os.system(raxml_parallel2)
		# ~ if args.gblocks_relaxed:
			# ~ raxml_parallel3 = "find %s -type f  -name '*_final_align_AA.aln.fas-gb' | parallel -j %s %s --all --msa {} --model PROTGTR+G --prefix {} --seed 888 --threads 1 --bs-metric tbe" %(path_to_macsed_align, args.cpu, raxml_script)
			# ~ os.system(raxml_parallel3)
		# ~ path_to_single_trees = path_to_macsed_align.replace('macsed_alignments/', 'single_locus_trees/')
		# ~ mkdir_raxml_out = "mkdir {}".format(path_to_single_trees) 
		# ~ os.system(mkdir_raxml_out)
		# ~ for root, dirs, files in os.walk(path_to_macsed_align, topdown=True):
			# ~ for f in files:
				# ~ if "raxml" in f:
					# ~ os.rename(path_to_macsed_align + f, path_to_single_trees + f)
		
		logging.info("****************************************************************************")
		logging.info("PERFORMING ALIGNMENTS CONCATENATION WITH Fasconcat")
		logging.info("****************************************************************************")
		# ~ path_to_supermatrix= path_to_macsed_align.replace('macsed_alignments/', 'supermatrix/')
		# ~ make_supermatrix_folder="mkdir {} ".format(path_to_supermatrix)
		# ~ make_supermatrix_dna= "mkdir {} ".format(path_to_supermatrix+ "supermatrix_dna/")
		# ~ make_supermatrix_aa="mkdir {} ".format(path_to_supermatrix + "supermatrix_aa/")
		# ~ os.system(make_supermatrix_folder)
		# ~ os.system(make_supermatrix_dna)
		# ~ os.system(make_supermatrix_aa)
		# ~ if args.gblocks_relaxed:
			# ~ make_supermatrix_dna= "mkdir {} ".format(path_to_supermatrix+ "supermatrix_gblocked_dna/")
			# ~ make_supermatrix_aa="mkdir {} ".format(path_to_supermatrix + "supermatrix_gblocked_aa/")
			# ~ os.system(make_supermatrix_dna)
			# ~ os.system(make_supermatrix_aa)
		# ~ path_to_supermatrix_dna = path_to_supermatrix +"supermatrix_dna/"
		# ~ path_to_supermatrix_aa = path_to_supermatrix + "supermatrix_aa/"
		# ~ copy_fasconcat = "cp {} {}".format(main_script_dir + "FASconCAT-G_v1.04.pl", path_to_supermatrix_dna)
		# ~ os.system(copy_fasconcat)
		# ~ copy_fasconcat = "cp {} {}".format(main_script_dir + "FASconCAT-G_v1.04.pl", path_to_supermatrix_aa)
		# ~ os.system(copy_fasconcat)
		# ~ copy_alignments_dna= "cp -r {}*macsed_final_align_NT.aln.fas {}".format(path_to_macsed_align, path_to_supermatrix_dna)
		# ~ os.system(copy_alignments_dna)
		# ~ copy_alignments_aa= "cp -r {}*macsed_final_align_AA.aln.fas {}".format(path_to_macsed_align, path_to_supermatrix_aa)
		# ~ os.system(copy_alignments_aa)
		# ~ if args.gblocks_relaxed:
			# ~ path_to_supermatrix_gblocked_dna = path_to_supermatrix +"/supermatrix_gblocked_dna/"
			# ~ path_to_supermatrix_gblocked_aa = path_to_supermatrix + "/supermatrix_gblocked_aa/"
			# ~ copy_fasconcat = "cp {} {}".format(main_script_dir + "FASconCAT-G_v1.04.pl", path_to_supermatrix_gblocked_dna)
			# ~ os.system(copy_fasconcat)
			# ~ copy_fasconcat = "cp {} {}".format(main_script_dir + "FASconCAT-G_v1.04.pl", path_to_supermatrix_gblocked_aa)
			# ~ os.system(copy_fasconcat)
			# ~ copy_alignments_dna= "cp -r {}*macsed_final_align_NT.aln.fas-gb {}".format(path_to_macsed_align, path_to_supermatrix_gblocked_dna)
			# ~ os.system(copy_alignments_dna)
			# ~ copy_alignments_aa= "cp -r {}*macsed_final_align_AA.aln.fas-gb {}".format(path_to_macsed_align, path_to_supermatrix_gblocked_aa)
			# ~ os.system(copy_alignments_aa)
			
		# ~ run_fasconcat = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(path_to_supermatrix_dna) 
		# ~ os.chdir(path_to_supermatrix_dna)
		# ~ os.system(run_fasconcat)
		# ~ os.rename('FcC_supermatrix.fas','FcC_supermatrix_NT.fasta')
		# ~ os.system("rm *.fas")
		# ~ run_fasconcat = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(path_to_supermatrix_aa) 
		# ~ os.chdir(path_to_supermatrix_aa)
		# ~ os.system(run_fasconcat)
		# ~ os.rename('FcC_supermatrix.fas','FcC_supermatrix_AA.fasta')
		# ~ os.system("rm *.fas")
		# ~ if args.gblocks_relaxed:
			# ~ for f in os.listdir(path_to_supermatrix_gblocked_dna):
				# ~ if f.endswith("-gb"):
					# ~ alignment_file = open(path_to_supermatrix_gblocked_dna + f,'rt')
					# ~ alignment_file_content = alignment_file.read()
					# ~ #delete spaces created by Gblocks in the alignemtns
					# ~ alignment_file_content = alignment_file_content.replace(' ','')
					# ~ alignment_file.close()
					# ~ alignment_file = open(path_to_supermatrix_gblocked_dna + f,'wt')
					# ~ alignment_file.write(alignment_file_content)
					# ~ alignment_file.close()
			# ~ for f in os.listdir(path_to_supermatrix_gblocked_dna):
				# ~ if f.endswith("-gb"):		
					# ~ os.rename(path_to_supermatrix_gblocked_dna + f, path_to_supermatrix_gblocked_dna + f +".fas" )
			# ~ run_fasconcat = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(path_to_supermatrix_gblocked_dna) 
			# ~ os.chdir(path_to_supermatrix_gblocked_dna)
			# ~ os.system(run_fasconcat)
			# ~ os.rename('FcC_supermatrix.fas','FcC_supermatrix_gblocked_NT.fasta')
			# ~ os.system("rm *.fas")
			# ~ for f in os.listdir(path_to_supermatrix_gblocked_aa):
				# ~ if f.endswith("-gb"):
					# ~ alignment_file = open(path_to_supermatrix_gblocked_aa + f,'rt')
					# ~ alignment_file_content = alignment_file.read()
					# ~ alignment_file_content = alignment_file_content.replace(' ','')
					# ~ alignment_file.close()
					# ~ alignment_file = open(path_to_supermatrix_gblocked_aa + f,'wt')
					# ~ alignment_file.write(alignment_file_content)
					# ~ alignment_file.close()
			# ~ for f in os.listdir(path_to_supermatrix_gblocked_aa):
				# ~ if f.endswith("-gb"):			
					# ~ os.rename(path_to_supermatrix_gblocked_aa + f, path_to_supermatrix_gblocked_aa + f +".fas" )
			# ~ run_fasconcat = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(path_to_supermatrix_gblocked_aa) 
			# ~ os.chdir(path_to_supermatrix_gblocked_aa)
			# ~ os.system(run_fasconcat)
			# ~ os.rename('FcC_supermatrix.fas','FcC_supermatrix_gblocked_AA.fasta')
			# ~ os.system("rm *.fas")
		# ~ os.chdir(main_script_dir)	
		logging.info("*********************************************************************")
		logging.info("RECONSTRUCTING SUPERMATRIX TREE WITH IQTREE2")
		logging.info("*********************************************************************")
		# ~ iqtree_script=main_script_dir + "iqtree2"
		# ~ iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s" %(iqtree_script, path_to_supermatrix_dna + 'FcC_supermatrix_NT.fasta' , path_to_supermatrix_dna + 'FcC_supermatrix_partition.txt', args.cpu)
		# ~ os.system(iqtree_on_supermatrix)
		# ~ iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s" %(iqtree_script, path_to_supermatrix_aa + 'FcC_supermatrix_AA.fasta' , path_to_supermatrix_aa + 'FcC_supermatrix_partition.txt', args.cpu)
		# ~ os.system(iqtree_on_supermatrix)
		# ~ iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s" %(iqtree_script, path_to_supermatrix_gblocked_dna + 'FcC_supermatrix_gblocked_NT.fasta' , path_to_supermatrix_gblocked_dna + 'FcC_supermatrix_partition.txt', args.cpu)
		# ~ os.system(iqtree_on_supermatrix)
		# ~ iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s" %(iqtree_script, path_to_supermatrix_gblocked_aa + 'FcC_supermatrix_gblocked_AA.fasta' , path_to_supermatrix_gblocked_aa + 'FcC_supermatrix_partition.txt', args.cpu)
		# ~ os.system(iqtree_on_supermatrix)
	
		# ~ ##### MAYBE RUN RAXML-NG USING THE IQTREE TOPOLOGY TO GET mbe BOOTSRAP VALUES???
		
		logging.info("***********************************************************")
		logging.info("RECONSTRUCTING SUPERTREE WITH ASTRAL")
		logging.info("***********************************************************")
		# ~ path_to_supertree = path_to_supermatrix.replace( 'supermatrix/','supertree/')
		# ~ make_supertree_folder ="mkdir {}".format(path_to_supertree)
		# ~ os.system(make_supertree_folder)
		# ~ make_supertree_dna = "mkdir {}".format(path_to_supertree + "supertree_dna")
		# ~ make_supertree_aa = "mkdir {}".format(path_to_supertree + "supertree_aa")
		# ~ make_supertree_gblocked_dna = "mkdir {}".format(path_to_supertree + "supertree_gblocked_dna")
		# ~ make_supertree_gblocked_aa = "mkdir {}".format(path_to_supertree + "supertree_gblocked_aa")
		# ~ os.system(make_supertree_dna)
		# ~ os.system(make_supertree_aa)
		# ~ os.system(make_supertree_gblocked_dna)
		# ~ os.system(make_supertree_gblocked_aa)

		# ~ copy_trees_dna= "cp -r {}*macsed_final_align_NT.aln.fas.raxml.bestTree {}".format(path_to_single_trees, path_to_supertree + "supertree_dna")
		# ~ copy_trees_aa= "cp -r {}*macsed_final_align_AA.aln.fas.raxml.bestTree {}".format(path_to_single_trees, path_to_supertree + "supertree_aa")
		# ~ copy_trees_gblocked_dna= "cp -r {}*macsed_final_align_NT.aln.fas-gb.raxml.bestTree {}".format(path_to_single_trees, path_to_supertree + "supertree_gblocked_dna")
		# ~ copy_trees_gblocked_aa= "cp -r {}*macsed_final_align_AA.aln.fas-gb.raxml.bestTree {}".format(path_to_single_trees, path_to_supertree + "supertree_gblocked_aa")
		# ~ os.system(copy_trees_dna)
		# ~ os.system(copy_trees_aa)
		# ~ os.system(copy_trees_gblocked_dna)
		# ~ os.system(copy_trees_gblocked_aa)
		
		# ~ os.chdir(path_to_supertree + "supertree_dna")
		# ~ cat_trees = "cat *.bestTree > cat_trees_dna.tre"
		# ~ os.system(cat_trees)
		# ~ os.system("rm -r *.bestTree")
		# ~ os.chdir(path_to_supertree + "supertree_aa")
		# ~ cat_trees = "cat *.bestTree > cat_trees_aa.tre"
		# ~ os.system(cat_trees)
		# ~ os.system("rm -r *.bestTree")
		# ~ os.chdir(path_to_supertree + "supertree_gblocked_dna")
		# ~ cat_trees = "cat *.bestTree > cat_trees_gblocked_dna.tre"
		# ~ os.system(cat_trees)
		# ~ os.system("rm -r *.bestTree")
		# ~ os.chdir(path_to_supertree + "supertree_gblocked_aa")
		# ~ cat_trees = "cat *.bestTree > cat_trees_gblocked_aa.tre"
		# ~ os.system(cat_trees)
		# ~ os.system("rm -r *.bestTree")
		# ~ os.chdir(main_script_dir)
		
		# ~ path_to_astral = main_script_dir + "ASTRAL/astral.5.7.7.jar"
		# ~ run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, path_to_supertree + "supertree_dna/cat_trees_dna.tre", path_to_supertree + "supertree_dna/astral_species_tree_dna.tree")
		# ~ os.system(run_astral)
		# ~ run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, path_to_supertree + "supertree_aa/cat_trees_aa.tre", path_to_supertree + "supertree_aa/astral_species_tree_aa.tree")
		# ~ os.system(run_astral)
		# ~ run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, path_to_supertree + "supertree_gblocked_dna/cat_trees_gblocked_dna.tre", path_to_supertree + "supertree_gblocked_dna/astral_species_tree_gblocked_dna.tree")
		# ~ os.system(run_astral)
		# ~ run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, path_to_supertree + "supertree_gblocked_aa/cat_trees_gblocked_aa.tre", path_to_supertree + "supertree_gblocked_aa/astral_species_tree_gblocked_aa.tree")
		# ~ os.system(run_astral)



if __name__=='__main__':
	logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	main()


