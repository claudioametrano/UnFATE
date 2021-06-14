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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import pandas

def select_best_reference_seq(prot_file_path, assemblies_path, cpu):
	"""in the assembly mode thase use exonerate_hits.py script from Hybpiper, there is no way to know what sequences from the protein file is the best for each markes,
		as we do not have the BLAST or the BWA results mappig the reads from target enrichment to those sequences, so we just BLAST all the reference for each gene
		on the assemblyes, the one with best bitscore is used to run exonerate_hits.py for that gene.
	"""
	ref_gene_list = []
	# Get a list of genes without redundancy
	for seq in SeqIO.parse(prot_file_path,"fasta"):
		regex_id = re.search("^GCA_[0-9]+\.[0-9]-([0-9]+at4890)", seq.id)
		if regex_id != None:
			ref_gene_list.append(regex_id.group(1))
	ref_gene_list = list(set(ref_gene_list))
	logging.info("Gene list: ")
	logging.info(ref_gene_list)
	# 	Generate a separate fasta for each gene in the reference sequnces file	
	for gene in ref_gene_list:
		with open(assemblies_path + gene+ "_ref.fasta", "a+") as gene_file:
			for seq in SeqIO.parse(prot_file_path,"fasta"):
				regex_id = re.search("^GCA_[0-9]+\.[0-9]-([0-9]+at4890)", seq.id)
				if regex_id != None:
					if regex_id.group(1) == gene:
						SeqIO.write(seq, gene_file, "fasta")
					else:
						pass	
	# build BLST databases for each alignment and then run every gene reference file against it usung tBLASTn									
	for f in os.listdir(assemblies_path):
		if f.endswith(".fna"):
			# Build a BLAST database for each of the assemblies
			# As assembly names ends in various characters os.pat.splitext is safer than .rstrip() that truncates some of the sample names
			make_db = "makeblastdb -in {} -dbtype nucl -out {}".format(assemblies_path + f, assemblies_path + os.path.splitext(f)[0])
			os.system(make_db)
			parallel_tblastn = 'find %s -type f -name "*_ref.fasta" | parallel -j %s tblastn -query {} -db %s -out {}%s -num_threads 1 -outfmt 6'%(assemblies_path, cpu, assemblies_path + os.path.splitext(f)[0], "___" + os.path.splitext(f)[0] + "_blastout.tsv")
			os.system(parallel_tblastn)
	# from the BLAST output file in .tsv extract the reference sequence with the best bitscore (better than e-value as it does not depend on database sequence number, easier than filter by multiple things e.g. query length percentid e-value etc.)
	for f in os.listdir(assemblies_path):
		# check if is one of the blast file output and if not empty
		if f.endswith("_blastout.tsv") and os.stat(assemblies_path + f).st_size != 0:
			#create dataframe
			df = pandas.read_table(assemblies_path + f , header=None)
			# assign columns names
			blast_columns = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
			df.columns = blast_columns
			# sort by bitscore (descending)
			df.sort_values(by='bitscore', ascending=False, inplace=True)
			# only retain the best hit
			top1 = df['qaccver'][0:1]
			regex = re.search("([0-9]+at4890)_ref.fasta___(.+?)_blastout.tsv",f)
			logging.info("For the assembly %s and gene %s reference sequence will be: "%(regex.group(2), regex.group(1)))
			logging.info(top1)
			# Write the best hit for every gene in an output file, taking the sequence from the original protein file
			output_file = open(assemblies_path + regex.group(2) + "_best_blast_scoring_reference_Hybpiper_format_aa.fas","a+")
			for rec in SeqIO.parse(prot_file_path, 'fasta'):
				if rec.id in str(top1):
					SeqIO.write(rec, output_file, 'fasta')				
			output_file.close()
			os.system("rm {}".format(assemblies_path + f))
		else:
			pass			
	# remove all the garbage needed to blast
	# find is used in order to pass to the remove command one file at a time, otherwise if there are too many files the rm command throws the error:"-bash: /bin/rm: Argument list too long" 
	extension = "*_ref.fasta"
	remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	os.system(remove_lot_of_files)
	#extension = "*_blastout.tsv"
	#remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	#os.system(remove_lot_of_files)
	extension = "*.nin"
	remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	os.system(remove_lot_of_files)
	extension = "*.nsq"
	remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	os.system(remove_lot_of_files)
	extension = "*.nhr"
	remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	os.system(remove_lot_of_files)
	extension = "*.ndb"
	remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	os.system(remove_lot_of_files)
	extension = "*.nto"
	remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	os.system(remove_lot_of_files)
	extension = "*.not"
	remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	os.system(remove_lot_of_files)
	extension = "*.ntf"
	remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
	os.system(remove_lot_of_files)
	return()

def run_exonerate_hits(file_, ref_seq_file):
	logging.info("Extracting genes from: " +file_)
	fline=open(file_).readline()
	regex_spades_header =re.search("^>NODE_[0-9]+_length_[0-9]+_cov_[0-9]+",fline)
	#if spades assembly, run exonerate_hits from HybPiper
	if regex_spades_header != None:
		os.system("python3 exonerate_hits.py {} --prefix {} {} ".format(ref_seq_file, os.path.splitext(file_)[0], file_))
	# else use the script version not using coverage information
	else:
		os.system("python3 exonerate_alt.py {} --prefix {} {} ".format(ref_seq_file, os.path.splitext(file_)[0], file_))
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
	logging.info("List of genes found: ")
	logging.info( genes_list)
	os.chdir(path_to_data)
	for g in genes_list:
		#logging.info("Building nucleotide alignment for gene {}".format(g))
		with open("Alignment_" + g + "_nucleotide.FAS",'a+') as alignment:
			for root, dirs, files in os.walk(path_to_data, topdown=True):
				for f in files:
					if f == g + ".FNA":
						os.chdir(root)
						#print("root: " + root +"/"+ f)
						with open(f, 'r') as gene:
							f_content = gene.read()
							os.chdir(path_to_data)
							#print("path to data:", path_to_data)
							alignment.write(f_content)	
	for g in genes_list:
		#logging.info("Building aminoacid alignment for gene {}".format(g))
		with open("Alignment_" + g + "_protein.FAS",'a+') as alignment:
			for root, dirs, files in os.walk(path_to_data, topdown=True):
				for f in files:
					if f == g + ".FAA":
						os.chdir(root)
						#print("root" + root +"/"+ f)
						with open(f, 'r') as gene:
							f_content = gene.read()
							os.chdir(path_to_data)
							#print("path to data: ", path_to_data)
							alignment.write(f_content)
	return()

def merge_alignments(path_to_assemblies, path_to_target_enrichment):
	if os.path.exists(path_to_target_enrichment):
		output_path = path_to_target_enrichment + '../'
	else:	
		output_path = path_to_assemblies + '../'
	os.makedirs(output_path + 'alignments')
	os.makedirs(output_path + 'alignments_merged')
	output_path = output_path + 'alignments/'
	output_path_merged = output_path.rstrip("/") + "_merged/"
	#print(output_path)
	if os.path.exists(path_to_assemblies):
		for filename in (os.listdir(path_to_assemblies)):
			if filename.endswith(".FAS"):
				os.rename(path_to_assemblies + filename , path_to_assemblies + filename.rstrip("FAS") + "assembly.FA")
	if os.path.exists(path_to_target_enrichment):
		for filename in (os.listdir(path_to_target_enrichment)):
			if filename.endswith(".FAS"):
				os.rename(path_to_target_enrichment + filename , path_to_target_enrichment + filename.rstrip("FAS") + "targenrich.FA")
	if os.path.exists(path_to_assemblies):
		for filename in (os.listdir(path_to_assemblies)):
			if filename.endswith("assembly.FA"):
				shutil.move(path_to_assemblies + filename , output_path + filename)
	if os.path.exists(path_to_target_enrichment):			 
		for filename in (os.listdir(path_to_target_enrichment)):
			if filename.endswith("targenrich.FA"):
				shutil.move(path_to_target_enrichment + filename , output_path + filename)
	file_list_total = os.listdir(output_path)
	# after the list of files is generated only retain in the name the gene name and if protein or nucleotide
	for index, item in  enumerate(file_list_total):
		file_list_total[index] = (file_list_total[index]).replace(".targenrich.FA", "")
		file_list_total[index] = (file_list_total[index]).replace(".assembly.FA", "")	
	#print(file_list_total)
	# then delete the double name doing a set
	file_list_total = list(set(file_list_total))
	#print(file_list_total)
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

def from_accession_to_species(csv_file, treefile):
	"""USE: takes as input a csv file with species name and accession: "GCA_009732865.1,Chalara longipes" and substitutes accession numbers with species names in the tree"""
	def substitute_tips(table, tree_content):
		#print(tree_content)
		count = 0
		for l in table:
			#print("element of the table", l[0])
			reg_tree = re.search(l[0],tree_content)
			if reg_tree is not None:
				#print("match in the tree: ",reg_tree.group())
				tree_content = tree_content.replace(reg_tree.group(), l[1] + "_" + l[0])
				count = count + 1
			else:
				print("Not found: ",l[0])
				tree_content = tree_content
		print("Substitutions done: ",count)		
		return(tree_content)

	"""main script"""			
	path=os.getcwd()
	# read the csv file, you get as many list as the rows in the .csv file
	table = open(csv_file) 
	csv_table = csv.reader(table, delimiter=',')
	my_tree = open(treefile, "r")
	my_tree_content = my_tree.read()
	output_tree = substitute_tips(csv_table, my_tree_content)
	output_file = open(treefile + "_SPECIES_NAME.tre", "w")
	output_file.write(output_tree)	
	output_file.close()
	return()

def check_arg(args=None):
	parser = argparse.ArgumentParser(description='UnFATE: the wrapper script that brings YOU from target enrichment sequencing data straight to phylogenetic tree inference! BE CAREFUL: At least one argument between assemblies and target enrichment is mandatory! See the readme file for data structure and additional info.')
	parser.add_argument('-bb', '--target_markers', default= '',
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
						help='Modifies some path in MACSE pipeline folder, use this argument only if is the first time you run the pipeline, then do not move the UnFATE folder, if you have to: git clone the repository again and use this arguement for the first run',
						)		
	parser.add_argument('-g', '--gblocks_relaxed', action= 'store_true', 
						help='Applies Gblocks with relaxed parameters (Talavera et al. 2007)',
						)	
	parser.add_argument('-n', '--ncbi_assemblies', nargs = '+', 
						help='Extracts the pre-mined NCBI assemblies genes, then it uses the selected taxonomic rank to only include the needed samples. Can take a list of ranks (e.g. Morchella Tuber Fuffaria)'
						)	
	parser.add_argument('--nargs', nargs='+')																				
	return parser.parse_args(args)
args = check_arg(sys.argv[1:])
print("Arguments are: ", args)


def main():
	#print(args)
	main_script_dir = os.path.realpath(__file__)
	main_script_dir = main_script_dir.rstrip("main_wrap.py")
	#print(main_script_dir)
	#print(args.target_enrichment_data)
	#print(args.assemblies)
	os.system('ulimit -n 1024000')
	
	
	if args.first_use == True:	
		# this git clones are not needed anymore as this software are now included directly in the UnFATE repository to avoid the wrapper to crash when a new version is released (could use github 
		# releases, but MACSE does not have releases). Kept here to easily clone them if needed
		#clone_hybpiper = 'git clone https://github.com/mossmatters/HybPiper.git'
		#os.system(clone_hybpiper)
		#logging.info("Hybpiper cloned")
		#clone_MACSE = 'git clone https://github.com/ranwez/MACSE_V2_PIPELINES.git'
		#os.system(clone_MACSE)
		#logging.info("MACSE alignment pipeline cloned")
		
		# commands to modity OMM_MACSE main script and utilities with the right path to them instead of installing Singularity to use it
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
					data = data.replace('macse="java -jar -Xmx${JAVA_MEM} ${LG_MACSE}"', 'macse="java -jar -Xms1g -Xmx2g ' + MACSE_utils_dir + '/macse_v2.03.jar"')
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
		# ASTRAL is now included in the repository, no need to clone (moreover the make script gives an error related to java version)				
		#clone_astral = "git clone https://github.com/smirarab/ASTRAL.git"
		#os.system(clone_astral)
		#os.chdir(main_script_dir + "ASTRAL/")
		#os.system("./make.sh")
		#os.chdir(main_script_dir)
		#logging.info("ASTRAL cloned")
	
	if args.ncbi_assemblies:
		logging.info("Obtaining target genes from pre-extracted assembly database ")
		path_to_premined = main_script_dir + "pre_mined_assemblies.tar.gz"
		unzip_premined_assemblies = "tar -zxf {}".format(path_to_premined)
		os.system(unzip_premined_assemblies)
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
		logging.info("********************************************************************************************************************")
		logging.info("          PERFORMING ASSEMBLIES DATA ANALYSIS WITH Exonerate (Slater & Birney 2005)       ")
		logging.info("********************************************************************************************************************")
		logging.info("... it can be time consuming, it depends on assembly dimension")
		logging.info('Path to assemblies '+path_to_assemblies)
		logging.info('Selecting the best reference sequence for each assembly by BLAST...')
		select_best_reference_seq(args.target_markers, args.assemblies, args.cpu)
		
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
		ref_list = []
		for k in os.listdir(path_to_assemblies):
			if k.endswith("_best_blast_scoring_reference_Hybpiper_format_aa.fas"):
				ref_list.append(path_to_assemblies + k)
		print(ref_list)
		list_of_list = []
		for z in pezizo_list:
			regex_fna = re.search("(.+?)\.fna", z)
			empty_list = []
			for v in ref_list:
				regex_ref = re.search("(.+?)_best_blast_scoring_reference_Hybpiper_format_aa\.fas", v)
				if regex_fna.group(1) == regex_ref.group(1):
					empty_list.append(z)
					empty_list.append(v)
					list_of_list.append(empty_list)
		print(list_of_list)
		logging.info("Running exonerate using exonerate_hits.py script from Hybpiper..")	
		args.cpu = int(args.cpu)
		pool = multiprocessing.Pool(processes=args.cpu)
		pool.starmap(run_exonerate_hits, list_of_list)
	
	if args.target_enrichment_data:
		logging.info("*****************************************************************************************************************************")
		logging.info("          TRIMMING TARGET ENRICHMENT FASTQ FILES  WITH TRIMMOMATC (Bolger et al. 2014)          ")
		logging.info("*****************************************************************************************************************************")
		logging.info('Path to TE data: '+args.target_enrichment_data)
		trimming_cmd = "python3 {}/trimmer.py -f {}".format(main_script_dir, args.target_enrichment_data)
		os.system(trimming_cmd)
		#Get namelist.txt from dataset directory
		namelist_cmd = 'python3 {}/getNameList.py -f {}'.format(main_script_dir, args.target_enrichment_data)
		os.system(namelist_cmd)
		namelist = 'namelist.txt'
		path_to_namelist = os.path.join(args.target_enrichment_data,namelist)
		logging.info("Gunzipping paired reads trimmed fastq archives")
		gunzip_fastq =' parallel gunzip ::: {}*_paired.fastq.gz'.format(args.target_enrichment_data) 
		os.system(gunzip_fastq)
		logging.info("******************************************************************************************************************************************")
		logging.info("           EXTRACTING GENES FROM TARGET ENRICHMENT DATA  WITH Hybpiper (Johnson et al. 2016)")
		logging.info("******************************************************************************************************************************************")
		os.chdir(args.target_enrichment_data)
		with open(path_to_namelist, 'r') as f:
			for line in f:
				logging.info("Processing sample:" + line)
				sample_path = args.target_enrichment_data + '/' + line.rstrip('\n') + '_R*.trimmed_paired.fastq'
				run_Hybpiper =  '{}HybPiper/reads_first.py -b {} -r {}  --prefix {} --cpu {} '.format(main_script_dir, args.target_markers, sample_path, line, args.cpu)
				os.system(run_Hybpiper)
		os.chdir(main_script_dir)
				
	if args.target_enrichment_data or args.assemblies:	
		logging.info("*********************************************")
		logging.info("          BUILDING FASTA FILES          ")
		logging.info("*********************************************")
		if args.assemblies:
			logging.info("Building alignments from assemblies data")
			get_alignment(args.assemblies)
		if args.target_enrichment_data:
			logging.info("Building alignments from target enrichment data")
			get_alignment(args.target_enrichment_data)
		merge_alignments(args.assemblies, args.target_enrichment_data)
		if 	args.target_enrichment_data:
			path_to_merged_alignments = args.target_enrichment_data.replace('target_enrichment/', 'alignments_merged/')
		if args.assemblies:
			path_to_merged_alignments = args.assemblies.replace('assemblies/', 'alignments_merged/')	
		
		if args.ncbi_assemblies:
			logging.info("****************************************************************************************************************************")
			logging.info("           ADDING SELECTED TAXONOMIC RANKS GENES FROM PRE-MINED ASSEMBLY DATABASE           ")
			logging.info("****************************************************************************************************************************")
			path_to_ranks = path_to_taxonomy.replace('Accession_plus_taxonomy_Pezizomycotina.txt','Accession_plus_taxonomy_reduced.txt')
			#print(path_to_ranks)
			# make the folder for the selected sample from NCBI database
			path_to_premined = main_script_dir + "pre_mined_assemblies/"
			path_to_premined_selected = path_to_merged_alignments.replace('alignments_merged/', 'ncbi_selected_genes/')
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
			# add the genes from the database to the alignments		
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
												
			# remove the temporary text file storing the selected taxa taxonomy and assembly, from the pre-calculated database folder		
			remove_taxonomy_file = "rm {}".format(path_to_ranks)
			os.system(remove_taxonomy_file)

		# Get rid of the trash strings after the accession number to be able to replace with speciens name later	
		# As OMM_MACSE will use soft masking to align and trim better get rig of all small case letters in the alignments before running MACSE pipeline								
		logging.info("Cleaning sequences names to only retain accession numbers...")
		logging.info("Converting all nucleotides to uppercase...")
		for f in os.listdir(path_to_merged_alignments):
					if f.endswith("_merged.fasta"):
						output_file = open(path_to_merged_alignments + f.rstrip("\.fasta") + "_headmod.fas","a")
						for seq in SeqIO.parse(path_to_merged_alignments + f,"fasta"):
							regex_id = re.search("(^GCA_[0-9]+.[0-9])_", seq.id)
							if regex_id is not None:
								# this strips the old header out (.id is only the accession in theory .description is the whole header instead
								seq.description = ""
								seq.id = regex_id.group(1)
								# replace sequence with the same but uppercase
								sequence = str(seq.seq).upper()
								# arrange sequence and id in a format that SeqIO can write to file
								record = SeqRecord(Seq(sequence), seq.id, "","")
								print(record)
								SeqIO.write(record, output_file,"fasta")
							else:
								sequence = str(seq.seq).upper()
								record = SeqRecord(Seq(sequence), seq.id, "","")
								print(record)
								SeqIO.write(record, output_file,"fasta")						
						output_file.close()
		
		logging.info("**********************************************************************************************************")
		logging.info("          COMPARING RETRIEVED GENES TO REFERENCE SEQUENCES LENGTH          ")
		logging.info("**********************************************************************************************************")
		markers_retrieved_percentage_script=main_script_dir + "markers_retrieved_percentage.py"		
		run_markers_retrieved_percentage = "python3 {} -b {} -f {} ".format(markers_retrieved_percentage_script, args.target_markers, path_to_merged_alignments)						
		os.system(run_markers_retrieved_percentage)
		
		logging.info("************************************************************************************************************************************************************")
		logging.info("          PERFORMING ALIGNMENT WITH OMM_MACSE with HMMcleaner filtering (Ranwez et al. 2018; Di Franco et al. 2019)         ")
		logging.info("************************************************************************************************************************************************************")
		if args.target_enrichment_data:
			path_to_merged_alignments = args.target_enrichment_data.replace('target_enrichment/', 'alignments_merged/')
		if args.assemblies:
			path_to_merged_alignments = args.assemblies.replace('assemblies/', 'alignments_merged/')
		MACSE_dir = main_script_dir + "MACSE_V2_PIPELINES/OMM_MACSE/"
		MACSE_script = MACSE_dir + "S_OMM_MACSE_V10.02.sh"
		# As OMM_MACSE uses soft masking put all the sequences in upper case before alignment and filtering 
		run_OMM_MACSE = 'find %s -type f -name "*_nucleotide_merged_headmod.fas" | parallel -j %s %s --out_dir {}_out --out_file_prefix macsed --in_seq_file {} --no_prefiltering --no_postfiltering --alignAA_soft MAFFT  --min_percent_NT_at_ends 0.01 ' %(path_to_merged_alignments, args.cpu, MACSE_script)
		os.system(run_OMM_MACSE)
		logging.info(run_OMM_MACSE)
		# move aligned files
		path_to_macsed_align = path_to_merged_alignments.replace('alignments_merged/','macsed_alignments/')
		make_align_fold = "mkdir {}".format(path_to_macsed_align)
		os.system(make_align_fold)
		for root, dirs, files in os.walk(path_to_merged_alignments, topdown=True):
			for f in files:
				if f.endswith("_final_align_NT.aln") or f.endswith("_final_align_AA.aln"):
					#print(root)
					#print(f)
					file_path = root +"/"+ f   
					regex1 =re.search("Alignment_([0-9]+at[0-9]+)_nucleotide_merged_headmod.fas_out",root)
					#os.rename both renames and moves files
					os.rename(file_path,  path_to_macsed_align + regex1.group(1) + f + ".fas")
		
		gblocks_path= main_script_dir + "Gblocks"
		if args.gblocks_relaxed:
			logging.info("***********************************************************************************************************")
			logging.info("          PERFORMING ALIGNMENT FILTERING WITH Gblocks (Castresana, 2000)            ")
			logging.info("***********************************************************************************************************")
			#print(path_to_macsed_align)
			run_gblocks("_final_align_NT.aln.fas","_final_align_AA.aln.fas", path_to_macsed_align, gblocks_path)	
			
		logging.info("******************************************************************************************************************")
		logging.info("          RECONSTRUCTING SINGLE MARKER TREES WITH IQTREE2 (Minh et al. 2020)           ")
		logging.info("******************************************************************************************************************")
		iqtree_script = main_script_dir + "iqtree2"
		#print(path_to_macsed_align)
		# DNA alignments
		if args.gblocks_relaxed:
			iqtree_parallel1 = "find %s -type f  -name '*_final_align_NT.aln.fas-gb' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
			os.system(iqtree_parallel1)
		else:
			iqtree_parallel = "find %s -type f  -name '*_final_align_NT.aln.fas' | parallel -j %s %s -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
			os.system(iqtree_parallel)
		# Amino acid alignments
		if args.gblocks_relaxed:
			iqtree_parallel3 = "find %s -type f  -name '*_final_align_AA.aln.fas-gb' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
			os.system(iqtree_parallel3)
		else:
			iqtree_parallel2 = "find %s -type f  -name '*_final_align_AA.aln.fas' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
			os.system(iqtree_parallel2)
		# move trees and other iqtree files to the dedicated folder		
		path_to_single_trees = path_to_macsed_align.replace('macsed_alignments/', 'single_locus_trees/')
		mkdir_raxml_out = "mkdir {}".format(path_to_single_trees) 
		os.system(mkdir_raxml_out)
		for root, dirs, files in os.walk(path_to_macsed_align, topdown=True):
			for f in files:
				if  f.endswith("treefile") or f.endswith("nex") or f.endswith("parttrees") or f.endswith("gz") or f.endswith("mldist") or f.endswith("log") or f.endswith("iqtree") or f.endswith("contree") or f.endswith("bionj") or f.endswith("best_scheme"):
					os.rename(path_to_macsed_align + f, path_to_single_trees + f)
	
		logging.info("************************************************************************************************************************")
		logging.info("          PERFORMING ALIGNMENTS CONCATENATION WITH Fasconcat (Kück & Longo, 2014)          ")
		logging.info("*************************************************************************************************************************")
		path_to_supermatrix= path_to_macsed_align.replace('macsed_alignments/', 'supermatrix/')
		make_supermatrix_folder="mkdir {} ".format(path_to_supermatrix)
		os.system(make_supermatrix_folder)
		if args.gblocks_relaxed:
			make_supermatrix_dna= "mkdir {} ".format(path_to_supermatrix+ "supermatrix_gblocked_dna/")
			make_supermatrix_aa="mkdir {} ".format(path_to_supermatrix + "supermatrix_gblocked_aa/")
			os.system(make_supermatrix_dna)
			os.system(make_supermatrix_aa)
		else:
			make_supermatrix_dna= "mkdir {} ".format(path_to_supermatrix+ "supermatrix_dna/")
			make_supermatrix_aa="mkdir {} ".format(path_to_supermatrix + "supermatrix_aa/")	
			os.system(make_supermatrix_dna)
			os.system(make_supermatrix_aa)
				
		if args.gblocks_relaxed:
			path_to_supermatrix_gblocked_dna = path_to_supermatrix +"supermatrix_gblocked_dna/"
			path_to_supermatrix_gblocked_aa = path_to_supermatrix + "supermatrix_gblocked_aa/"
			copy_fasconcat = "cp {} {}".format(main_script_dir + "FASconCAT-G_v1.04.pl", path_to_supermatrix_gblocked_dna)
			os.system(copy_fasconcat)
			copy_fasconcat = "cp {} {}".format(main_script_dir + "FASconCAT-G_v1.04.pl", path_to_supermatrix_gblocked_aa)
			os.system(copy_fasconcat)
			copy_alignments_dna= "cp -r {}*macsed_final_align_NT.aln.fas-gb {}".format(path_to_macsed_align, path_to_supermatrix_gblocked_dna)
			os.system(copy_alignments_dna)
			copy_alignments_aa= "cp -r {}*macsed_final_align_AA.aln.fas-gb {}".format(path_to_macsed_align, path_to_supermatrix_gblocked_aa)
			os.system(copy_alignments_aa)
		else:
			path_to_supermatrix_dna = path_to_supermatrix +"supermatrix_dna/"
			path_to_supermatrix_aa = path_to_supermatrix + "supermatrix_aa/"
			copy_fasconcat = "cp {} {}".format(main_script_dir + "FASconCAT-G_v1.04.pl", path_to_supermatrix_dna)
			os.system(copy_fasconcat)
			copy_fasconcat = "cp {} {}".format(main_script_dir + "FASconCAT-G_v1.04.pl", path_to_supermatrix_aa)
			os.system(copy_fasconcat)
			copy_alignments_dna= "cp -r {}*macsed_final_align_NT.aln.fas {}".format(path_to_macsed_align, path_to_supermatrix_dna)
			os.system(copy_alignments_dna)
			copy_alignments_aa= "cp -r {}*macsed_final_align_AA.aln.fas {}".format(path_to_macsed_align, path_to_supermatrix_aa)
			os.system(copy_alignments_aa)
				
		if args.gblocks_relaxed:
			for f in os.listdir(path_to_supermatrix_gblocked_dna):
				if f.endswith("-gb"):
					alignment_file = open(path_to_supermatrix_gblocked_dna + f,'rt')
					alignment_file_content = alignment_file.read()
					#delete spaces created by Gblocks in the alignemtns
					alignment_file_content = alignment_file_content.replace(' ','')
					alignment_file.close()
					alignment_file = open(path_to_supermatrix_gblocked_dna + f,'wt')
					alignment_file.write(alignment_file_content)
					alignment_file.close()
			for f in os.listdir(path_to_supermatrix_gblocked_dna):
				if f.endswith("-gb"):		
					os.rename(path_to_supermatrix_gblocked_dna + f, path_to_supermatrix_gblocked_dna + f +".fas" )
			run_fasconcat = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(path_to_supermatrix_gblocked_dna) 
			os.chdir(path_to_supermatrix_gblocked_dna)
			os.system(run_fasconcat)
			os.rename('FcC_supermatrix.fas','FcC_supermatrix_gblocked_NT.fasta')
			os.rename('FcC_supermatrix_partition.txt','FcC_supermatrix_partition_gblocked_NT.txt')
			os.rename('FcC_info.xls','FcC_info_gblocked_NT.xls')
			os.system("rm *.fas")
			for f in os.listdir(path_to_supermatrix_gblocked_aa):
				if f.endswith("-gb"):
					alignment_file = open(path_to_supermatrix_gblocked_aa + f,'rt')
					alignment_file_content = alignment_file.read()
					alignment_file_content = alignment_file_content.replace(' ','')
					alignment_file.close()
					alignment_file = open(path_to_supermatrix_gblocked_aa + f,'wt')
					alignment_file.write(alignment_file_content)
					alignment_file.close()
			for f in os.listdir(path_to_supermatrix_gblocked_aa):
				if f.endswith("-gb"):			
					os.rename(path_to_supermatrix_gblocked_aa + f, path_to_supermatrix_gblocked_aa + f +".fas" )
			run_fasconcat = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(path_to_supermatrix_gblocked_aa) 
			os.chdir(path_to_supermatrix_gblocked_aa)
			os.system(run_fasconcat)
			os.rename('FcC_supermatrix.fas','FcC_supermatrix_gblocked_AA.fasta')
			os.rename('FcC_supermatrix_partition.txt','FcC_supermatrix_partition_gblocked_AA.txt')
			os.rename('FcC_info.xls','FcC_info_gblocked_AA.xls')
			os.system("rm *.fas")
		else:			
			run_fasconcat = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(path_to_supermatrix_dna) 
			os.chdir(path_to_supermatrix_dna)
			os.system(run_fasconcat)
			os.rename('FcC_supermatrix.fas','FcC_supermatrix_NT.fasta')
			os.rename('FcC_supermatrix_partition.txt','FcC_supermatrix_partition_NT.txt')
			os.rename('FcC_info.xls','FcC_info_NT.xls')
			os.system("rm *.fas")
			run_fasconcat = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(path_to_supermatrix_aa) 
			os.chdir(path_to_supermatrix_aa)
			os.system(run_fasconcat)
			os.rename('FcC_supermatrix.fas','FcC_supermatrix_AA.fasta')
			os.rename('FcC_supermatrix_partition.txt','FcC_supermatrix_partition_AA.txt')
			os.rename('FcC_info.xls','FcC_info_AA.xls')
			os.system("rm *.fas")	
		os.chdir(main_script_dir)
			
		logging.info("*************************************************************************************************************")
		logging.info("          RECONSTRUCTING SUPERMATRIX TREE WITH IQTREE2  (Minh et al. 2020)          ")
		logging.info("*************************************************************************************************************")
		iqtree_script=main_script_dir + "iqtree2"
		if args.gblocks_relaxed:
			#path_to_supermatrix_gblocked_dna = path_to_supermatrix +"supermatrix_gblocked_dna/"
			#path_to_supermatrix_gblocked_aa = path_to_supermatrix + "supermatrix_gblocked_aa/"
			iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s" %(iqtree_script, path_to_supermatrix_gblocked_dna + 'FcC_supermatrix_gblocked_NT.fasta' , path_to_supermatrix_gblocked_dna + 'FcC_supermatrix_partition_gblocked_NT.txt', args.cpu)
			os.system(iqtree_on_supermatrix)
			iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s" %(iqtree_script, path_to_supermatrix_gblocked_aa + 'FcC_supermatrix_gblocked_AA.fasta' , path_to_supermatrix_gblocked_aa + 'FcC_supermatrix_partition_gblocked_AA.txt', args.cpu)
			os.system(iqtree_on_supermatrix)
		else:		
			#path_to_supermatrix_dna = path_to_supermatrix +"supermatrix_dna/"
			#path_to_supermatrix_aa = path_to_supermatrix + "supermatrix_aa/"						
			iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s" %(iqtree_script, path_to_supermatrix_dna + 'FcC_supermatrix_NT.fasta' , path_to_supermatrix_dna + 'FcC_supermatrix_partition_NT.txt', args.cpu)
			os.system(iqtree_on_supermatrix)
			iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s" %(iqtree_script, path_to_supermatrix_aa + 'FcC_supermatrix_AA.fasta' , path_to_supermatrix_aa + 'FcC_supermatrix_partition_AA.txt', args.cpu)
			os.system(iqtree_on_supermatrix)

		##### MAYBE RUN RAXML-NG USING THE IQTREE TOPOLOGY TO GET mbe BOOTSTRAP VALUES (very slow, as it is a complete boostrap?
		
		logging.info("******************************************************************************************************")
		logging.info("            RECONSTRUCTING SUPERTREE WITH ASTRAL (Zhang et al. 2018)")
		logging.info("******************************************************************************************************")
		path_to_supertree = path_to_supermatrix.replace( 'supermatrix/','supertree/')
		make_supertree_folder ="mkdir {}".format(path_to_supertree)
		os.system(make_supertree_folder)
		if args.gblocks_relaxed:
			make_supertree_gblocked_dna = "mkdir {}".format(path_to_supertree + "supertree_gblocked_dna")
			make_supertree_gblocked_aa = "mkdir {}".format(path_to_supertree + "supertree_gblocked_aa")
			os.system(make_supertree_gblocked_dna)
			os.system(make_supertree_gblocked_aa)
			copy_trees_gblocked_dna= "cp -r {}*macsed_final_align_NT.aln.fas-gb.treefile {}".format(path_to_single_trees, path_to_supertree + "supertree_gblocked_dna")
			copy_trees_gblocked_aa= "cp -r {}*macsed_final_align_AA.aln.fas-gb.treefile {}".format(path_to_single_trees, path_to_supertree + "supertree_gblocked_aa")
			os.system(copy_trees_gblocked_dna)
			os.system(copy_trees_gblocked_aa)
			os.chdir(path_to_supertree + "supertree_gblocked_dna")
			cat_trees = "cat *.treefile > cat_trees_gblocked_dna.tre"
			os.system(cat_trees)
			os.system("rm -r *.treefile")
			os.chdir(path_to_supertree + "supertree_gblocked_aa")
			cat_trees = "cat *.treefile > cat_trees_gblocked_aa.tre"
			os.system(cat_trees)
			os.system("rm -r *.treefile")
			os.chdir(main_script_dir)
			path_to_astral = main_script_dir + "ASTRAL/astral.5.7.7.jar"
			run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, path_to_supertree + "supertree_gblocked_dna/cat_trees_gblocked_dna.tre", path_to_supertree + "supertree_gblocked_dna/astral_species_tree_gblocked_dna.tree")
			os.system(run_astral)
			run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, path_to_supertree + "supertree_gblocked_aa/cat_trees_gblocked_aa.tre", path_to_supertree + "supertree_gblocked_aa/astral_species_tree_gblocked_aa.tree")
			os.system(run_astral)
		else:
			make_supertree_dna = "mkdir {}".format(path_to_supertree + "supertree_dna")
			make_supertree_aa = "mkdir {}".format(path_to_supertree + "supertree_aa")
			os.system(make_supertree_dna)
			os.system(make_supertree_aa)
			copy_trees_dna= "cp -r {}*macsed_final_align_NT.aln.fas.treefile {}".format(path_to_single_trees, path_to_supertree + "supertree_dna")
			copy_trees_aa= "cp -r {}*macsed_final_align_AA.aln.fas.treefile {}".format(path_to_single_trees, path_to_supertree + "supertree_aa")
			os.system(copy_trees_dna)
			os.system(copy_trees_aa)
			os.chdir(path_to_supertree + "supertree_dna")
			cat_trees = "cat *.treefile > cat_trees_dna.tre"
			os.system(cat_trees)
			os.system("rm -r *.treefile")
			os.chdir(path_to_supertree + "supertree_aa")
			cat_trees = "cat *.treefile > cat_trees_aa.tre"
			os.system(cat_trees)
			os.system("rm -r *.treefile")
			path_to_astral = main_script_dir + "ASTRAL/astral.5.7.7.jar"
			run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, path_to_supertree + "supertree_dna/cat_trees_dna.tre", path_to_supertree + "supertree_dna/astral_species_tree_dna.tree")
			os.system(run_astral)
			run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, path_to_supertree + "supertree_aa/cat_trees_aa.tre", path_to_supertree + "supertree_aa/astral_species_tree_aa.tree")
			os.system(run_astral)
				
		logging.info("****************************************************************")
		logging.info("          COPYING TREES TO final_trees FOLDER           ")
		logging.info("****************************************************************")
		# Copy the tree file from IQTREE and ASTRAL to a new folder
		path_to_finaltrees = path_to_supertree.replace('supertree/','final_trees/')
		make_finaltrees_folder ="mkdir {}".format(path_to_finaltrees)
		os.system(make_finaltrees_folder)
		for root, dirs, files in os.walk(path_to_supermatrix, topdown=True):
			for f in files:
				if f.endswith("treefile"):
					#print(root +f)
					copy_tree = 'cp {} {}'.format(root +"/"+f, path_to_finaltrees)
					os.system(copy_tree)
		for root, dirs, files in os.walk(path_to_supertree, topdown=True):
			for f in files:
				if f.startswith("astral_species_tree"):
					copy_tree = 'cp {} {}'.format(root +"/"+ f, path_to_finaltrees)
					os.system(copy_tree)	
				
		logging.info("*************************************************************************************************")
		logging.info("          CONVERTIG ACCESSIONS IN THE TREES, if any, TO SPECIES NAME")
		logging.info("*************************************************************************************************")
		# Open one of the suprematrices, making a list of accession sample name
		if args.gblocks_relaxed:
			supermatrix_file = path_to_supermatrix_gblocked_dna + 'FcC_supermatrix_gblocked_NT.fasta'
		else:	
			supermatrix_file = path_to_supermatrix_dna + 'FcC_supermatrix_NT.fasta'
		supermatrix_accession_file = path_to_finaltrees + 'Accessions.csv'
		with open(supermatrix_accession_file, 'w') as accessions:
			with open(supermatrix_file, 'r') as supermatrix:
				supermatrix_content = supermatrix.readlines()
				for line in supermatrix_content:
					regex = re.search("^>(GCA_[0-9]+\.[0-9])", line)
					if regex:
						accessions.write(regex.group(1)+"\n")	
					else:
						pass	
		# Add taxonomy to the accessions retrieved (get_taxonomy_with edirect script), select species name  and format the .csv file 
		accessions_plus_taxonomy_file = path_to_finaltrees + 'Accessions_plus_taxonomy.csv'
		get_taxonomy_script = main_script_dir + "get_taxonomy_with_edirect.py" 
		get_taxonomy = "python3 {} --accession_file {} --out_file {}".format(get_taxonomy_script, supermatrix_accession_file,  accessions_plus_taxonomy_file)
		os.system(get_taxonomy)
		# Clean the taxonomy file to get only "Accession,speciesname"
		accession_species_file = path_to_finaltrees + 'Accessions_plus_species.csv'
		with open(accession_species_file, 'w') as species:
			with open(accessions_plus_taxonomy_file, 'r') as acc_taxo:
				acc_taxo_cont = acc_taxo.readlines()
				for line in acc_taxo_cont:
					regex = re.search("^(GCA_[0-9]+\.[0-9],\w+ \w+)", line)
					if regex != None:
						#print(regex.group(1))
						species.write(regex.group(1) +"\n")	
					else:
						logging.warning("The following line does not have the expected format for species name, weird strain name format!")
						logging.warning(line)
		# Use the "Speciesname, Accession" csv file to substitute the Accession numbers with species names using the funcion "from_accession_to_species"
		# Final names after the substitution will be: "speciesname_accessionnumber"
		for treefile in os.listdir(path_to_finaltrees):
			if treefile.endswith("treefile") or treefile.endswith("tree"):		
				from_accession_to_species(accession_species_file, path_to_finaltrees + treefile)
	logging.info("PIPELINE COMPLETED!")
if __name__=='__main__':
	logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	main()


