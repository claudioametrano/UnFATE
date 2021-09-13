import shutil
import re
import os
import sys
import argparse
#import shlex
import multiprocessing
#import itertools
import subprocess
import logging
from os import path
from Bio import SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import pandas
from glob import glob

def blast_individual_sample(sample, genes, assemblies_path):
	fna = os.path.join(assemblies_path, sample + ".fna")
	# Build a BLAST database for each of the assemblies
	make_db = "makeblastdb -in {} -dbtype nucl -out {}".format(fna, os.path.join(assemblies_path, sample), os.path.join(assemblies_path, ))
	os.system(make_db)
	sample_best_path = os.path.join(assemblies_path, sample + "_best_blast_scoring_reference_Hybpiper_format_aa.fas")
	for gene in genes:
		ref_path = os.path.join(assemblies_path, gene + "_ref.fasta")
		blastout_path = os.path.join(assemblies_path, gene + "_ref.fasta___" + sample + "_blastout.tsv")
		tblastn_command = "tblastn -query {} -db {} -out {} -num_threads 1 -outfmt 6".format(ref_path, os.path.join(assemblies_path, sample), blastout_path)
		os.system(tblastn_command)

		#blast_columns isn't needed anymore, but it's nice for reference.
		#blast_columns = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

		with open(blastout_path) as blastFile, open(sample_best_path, "a+") as bestFile, open(ref_path) as refFile:
			bestScoring = ""
			bestScore = 0
			for line in blastFile:
				line = line.strip().split("\t")
				if float(line[-1]) > bestScore:
					bestScore = float(line[-1]) #bitscore
					bestScoring = line[0] #qaccver
			regex = re.search(r"(.+?)_ref\.fasta___(.+?)_blastout.tsv", blastout_path.split("/")[-1])
			logging.info("For assembly:{} and gene:{}, the reference sequence will be: {}".format(regex.group(2), regex.group(1), bestScoring))
			for record in SeqIO.parse(refFile, "fasta"):
				if record.id in bestScoring:
					SeqIO.write(record, bestFile, "fasta")
		os.system("rm {}".format(blastout_path))


def select_best_reference_seq(prot_file_path, assemblies_path, cpu):
	"""in the assembly mode thase use exonerate_hits.py script from Hybpiper, there is no way to know what sequences from the protein file is the best for each markes,
		as we do not have the BLAST or the BWA results mappig the reads from target enrichment to those sequences, so we just BLAST all the reference for each gene
		on the assemblies, the one with best bitscore is used to run exonerate_hits.py for that gene.
	"""
	ref_gene_list = []
	# Get a list of genes without redundancy
	for seq in SeqIO.parse(prot_file_path,"fasta"):
#		regex_id = re.search("^GCA_[0-9]+\.[0-9]-([0-9]+at4890)", seq.id)
#		if regex_id != None:
#			ref_gene_list.append(regex_id.group(1))
		geneName = seq.id.strip().split("-")[1]
		ref_gene_list.append(geneName)
	ref_gene_list = list(set(ref_gene_list))
	logging.info("Gene list: ")
	logging.info(ref_gene_list)
	# 	Generate a separate fasta for each gene in the reference sequnces file
	for gene in ref_gene_list:
		with open(assemblies_path + gene+ "_ref.fasta", "a+") as gene_file:
			for seq in SeqIO.parse(prot_file_path,"fasta"):
				geneName = seq.id.strip().split("-")[1]
				if geneName == gene:
					SeqIO.write(seq, gene_file, "fasta")
				else:
					pass
	fnas = glob(os.path.join(assemblies_path, "*.fna"))
	samples = [fna.split("/")[-1][:-4] for fna in fnas] #fna = "/gar/abc.fna", keep "abc"
	list_of_lists = [[sample, ref_gene_list, assemblies_path] for sample in samples] #[["abc", ["1", "2"], "/gar/assemblies/"],...]

	pool = multiprocessing.Pool(processes=int(args.cpu))
	pool.starmap(blast_individual_sample, list_of_lists)

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

#add function to get sample names that mirrors get_fastas_exonerate()
def get_names(path_to_data, isAssemblies):
	if isAssemblies:
		search_location = os.path.join(path_to_data, "*", "sequences")
	else:
		search_location = os.path.join(path_to_data, "*", "exonerate_genelist.txt")
	found_samples = set()
	for result in glob(search_location):
		print(result)
		found_samples.add(result.split("/")[-2])
	print(found_samples)
	return list(found_samples)


def get_fastas_exonerate(path_to_data, isAssemblies):
	for moleculeType in ["FNA", "FAA"]:
		if isAssemblies:
			search_location = os.path.join(path_to_data, "*", "sequences", moleculeType, "*")
		else:
			search_location = os.path.join(path_to_data, "*", "*", "*", "sequences", moleculeType, "*")
		for filename in glob(search_location):
			geneName = filename.split("/")[-1][:-4]
			if moleculeType == "FNA":
				with open(filename) as inFile, open(os.path.join(args.out, "fastas", "Alignment_" + geneName + "_nucleotide_merged.fasta"), 'a') as outFile:
					outFile.write(inFile.read())
			if moleculeType == "FAA":
				with open(filename) as inFile, open(os.path.join(args.out, "fastas", "Alignment_" + geneName + "_protein_merged.fasta"), 'a') as outFile:
					outFile.write(inFile.read())

def select_scores(scores, user_samples):
	num_db_samples_to_add = 10 #per "cluster"
	non_user_count = 0
	close_samples = []
	for pair in scores: #(sample, score)
		if pair[0] in user_samples:
			close_samples.append(pair[0])
		else:
			non_user_count += 1
			close_samples.append(pair[0])
		if non_user_count > num_db_samples_to_add:
			return close_samples

def get_score(seq1, seq2, index, calc, total):
	score = calc._pairwise(seq1, seq2)
	if index % 100 == 0:
		print("Processed {} of {}".format(index, total))
	return((index, score))

def find_similar_samples(query, user_samples, data_dir, cpus):
	logging.info("Finding similar samples to {}".format(query))
	logging.info('Progress shown by "Processed <progress> of <number of sequences>", may not show progress linearly, especially when being logged')
	aln = AlignIO.read(open(os.path.join(data_dir, "FcC_supermatrix.fas")), "fasta")
	for record in aln:
		record.seq = record.seq.upper()

	#distances will be calculated with a basic 2-parameter model
	calc = DistanceCalculator("trans", skip_letters=["N", "-", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V"])

	queryIndex = [record.id for record in aln].index(query)
	sequences = [record.seq for record in aln]
	ids = [record.id for record in aln]
	querySeq = sequences[queryIndex]

	scores = [] #list of tuples to alllow sorting

	with multiprocessing.Pool(cpus) as p:
		scores=p.starmap(get_score, [(querySeq, seq, i, calc, len(sequences)) for i,seq in enumerate(sequences)])
	scores.sort(key = lambda x: x[1])

	idScores = [(ids[elements[0]], elements[1]) for elements in scores]
	return select_scores(idScores, user_samples)

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
			b1 = str(round(count * fraction1))
			b2 = str(round(count * fraction2))	
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

def set_up_directories():
	#if it exists and is not an absolute path
	if args.assemblies:
		if not args.assemblies.startswith("/") or args.assemblies.startswith("~"):
			args.assemblies = os.path.abspath(args.assemblies) + "/"

	if args.target_enrichment_data:
		if not args.target_enrichment_data.startswith("/") or args.target_enrichment_data.startswith("~"):
			args.target_enrichment_data = os.path.abspath(args.target_enrichment_data) + "/"

	if args.whole_genome_data:
		if not args.whole_genome_data.startswith("/") or args.whole_genome_data.startswith("~"):
			args.whole_genome_data = os.path.abspath(args.whole_genome_data) + "/"

	if args.target_markers:
		if not args.target_markers.startswith("/") or args.target_markers.startswith("~"):
			args.target_markers = os.path.abspath(args.target_markers)

	if args.out:
		if not args.out.startswith("/") or args.out.startswith("~"):
			args.out = os.path.abspath(args.out) + "/"

	#if the user has specified an output directory, move their assemblies and target enrichment data
	#into the output directory (symlink), then run using those directories.
	if args.out:
		#make args.out directory if it doesn't already exist
		if not os.path.isdir(args.out):
			os.mkdir(args.out)
			if args.assemblies:
				os.mkdir(os.path.join(args.out, "assemblies"))
			if args.target_enrichment_data:
				os.mkdir(os.path.join(args.out, "target_enrichment_data"))
			if args.whole_genome_data:
				os.mkdir(os.path.join(args.out, "whole_genome_data"))
		else:
			if args.target_enrichment_data:
				if not os.path.isdir(os.path.join(args.out, "target_enrichment_data")):
					os.mkdir(os.path.join(args.out, "target_enrichment_data"))
			if args.assemblies:
				if not os.path.isdir(os.path.join(args.out, "assemblies")):
					os.mkdir(os.path.join(args.out, "assemblies"))
			if args.whole_genome_data:
				if not os.path.isdir(os.path.join(args.out, "whole_genome_data")):
					os.mkdir(os.path.join(args.out, "whole_genome_data"))

		if args.target_enrichment_data:
			for file in os.listdir(args.target_enrichment_data):
				if os.path.exists(os.path.join(args.out, "target_enrichment_data", file)):
					print("Path already exists")
					continue #symlink already exists
				if file.lower().endswith(".fastq") or file.lower().endswith(".fastq.gz"):
					os.symlink(os.path.join(args.target_enrichment_data, file), os.path.join(args.out, "target_enrichment_data", file))
			args.target_enrichment_data = os.path.join(args.out, "target_enrichment_data", "") #the empty one causes a trailing /

		if args.assemblies:
			for file in os.listdir(args.assemblies):
				if os.path.exists(os.path.join(args.out, "assemblies", file)):
					print("Path already exists")
					continue
				if file.lower().endswith(".fna") or file.lower().endswith(".fasta") or file.lower().endswith(".fna.gz"):
					os.symlink(os.path.join(args.assemblies, file), os.path.join(args.out, "assemblies", file))
			args.assemblies = os.path.join(args.out, "assemblies", "")

		if args.whole_genome_data:
			for file in os.listdir(args.whole_genome_data):
				if os.path.exists(os.path.join(args.out, "whole_genome_data", file)):
					print("Path already exists")
					continue #symlink already exists
				if file.lower().endswith(".fastq") or file.lower().endswith(".fastq.gz"):
					os.symlink(os.path.join(args.whole_genome_data, file), os.path.join(args.out, "whole_genome_data", file))
			args.whole_genome_data = os.path.join(args.out, "whole_genome_data", "") #the empty one causes a trailing /

def trim_and_get_namelist(main_script_dir, data_dir):
	trimming_cmd = "python3 {}/trimmer.py -f {} -c {}".format(main_script_dir, data_dir, args.cpu)
	os.system(trimming_cmd)
	#Get namelist.txt from dataset directory
	namelist_cmd = 'python3 {}/getNameList.py -f {}'.format(main_script_dir, data_dir)
	os.system(namelist_cmd)

def run_hybpiper(main_script_dir, data_dir, namelist):
	os.chdir(data_dir)
	with open(namelist, 'r') as f:
		for line in f:
			logging.info("Processing sample:" + line)
			if len(glob(data_dir + '/' + line.rstrip('\n') + '_R*.trimmed_paired.fastq')) == 2:
				sample_path = data_dir + '/' + line.rstrip('\n') + '_R*.trimmed_paired.fastq'
			else:
				sample_path = data_dir + '/' + line.rstrip('\n') + '_SE.trimmed.fastq'
			run_Hybpiper =  '{}HybPiper/reads_first.py -b {} -r {}  --prefix {} --cpu {} '.format(main_script_dir, args.target_markers, sample_path, line.strip(), args.cpu)
			logging.info("running HybPiper with: " + run_Hybpiper)
			os.system(run_Hybpiper)
			clean_command = "{}HybPiper/cleanup.py {}".format(main_script_dir, line.strip())
			os.system(clean_command)
	os.chdir(main_script_dir)

def run_exonerate_hits(file_, ref_seq_file, memory, threshold):
	logging.info("Extracting genes from: " +file_)
	fline=open(file_).readline()
	regex_spades_header =re.search("^>NODE_[0-9]+_length_[0-9]+_cov_[0-9]+",fline)

	print("EXONERATE MEMORY PER SAMPLE IS: {}GB".format(memory))
	#if spades assembly, run exonerate_hits from HybPiper
	if regex_spades_header != None:
		exonerate_command = "python3 {}exonerate_hits.py -m {} -t {} {} --prefix {} {} ".format(main_script_dir, memory, threshold, ref_seq_file, os.path.splitext(file_)[0], file_)
		print(exonerate_command)
		os.system(exonerate_command)
	# else use the script version not using coverage information
	else:
		exonerate_command = "python3 {}exonerate_alt.py -m {} -t {} {} --prefix {} {} ".format(main_script_dir, memory, threshold, ref_seq_file, os.path.splitext(file_)[0], file_)
		print(exonerate_command)
		os.system(exonerate_command)


def run_exonerate(data_dir):
	path_to_assemblies = data_dir
	logging.info("... it can be time consuming, it depends on assembly dimension")
	logging.info('Path to assemblies ' + path_to_assemblies)
	logging.info('Selecting the best reference sequence for each assembly by BLAST...')
	for root, dirs, files in os.walk(path_to_assemblies, topdown=True):
		for name in files:
			if name.endswith(".fna.gz"): # or name.endswith(".fasta.gz"):
				os.system("gunzip -f "+ path_to_assemblies + name)

	select_best_reference_seq(args.target_markers, path_to_assemblies, args.cpu)

	assemblies_count = 0
	pezizo_list = []
	for root, dirs, files in os.walk(path_to_assemblies, topdown=True):
		for name in files:
			if name.endswith(".fna"): #or name.endswith(".fasta"):
				pezizo_list.append(root + name)
				assemblies_count += 1
	#print("Samples are: ", pezizo_list)
	ref_list = []
	for k in os.listdir(path_to_assemblies):
		if k.endswith("_best_blast_scoring_reference_Hybpiper_format_aa.fas"):
			ref_list.append(path_to_assemblies + k)

	memory = 0
	if assemblies_count < int(args.cpu):
		memory = int(args.exonerate_mem / assemblies_count)
	else:
		memory = int(args.exonerate_mem / int(args.cpu))
	#print(ref_list)
	list_of_list = []
	for z in pezizo_list:
		#regex_fna = re.search("(.+?)\.fna", z)
		basename = ".".join(z.split(".")[:-1])
		empty_list = []
		for v in ref_list:
			regex_ref = re.search("(.+?)_best_blast_scoring_reference_Hybpiper_format_aa\.fas", v)
			#if regex_fna.group(1) == regex_ref.group(1):
			if basename == regex_ref.group(1):
				empty_list.append(z)
				empty_list.append(v)
				empty_list.append(memory)
				empty_list.append(int(args.threshold))
				list_of_list.append(empty_list)
	#print(list_of_list)
	logging.info("Running exonerate using exonerate_hits.py script from Hybpiper..")
	args.cpu = int(args.cpu)
	pool = multiprocessing.Pool(processes=args.cpu)
	pool.starmap(run_exonerate_hits, list_of_list)

def check_coverage(data_dir):
	blastFiles = glob("*best_blast_scoring*")
	samples = ["_".join(file.split("_")[:-7]) for file in blastFiles] #there might be _ in sample, :-7 removes consistent suffix
	allGenes = {} #{geneName:len}
	referenceRecords = SeqIO.parse(blastFiles[0], "fasta")
	samplesToRerun = {} #{sample:[genes]} to be returned
	overallThreshold = 0.70
	geneThreshold = 0.70

	for record in referenceRecords:
		geneName = record.id.split("-")[1]
		geneLength = len(record.seq)
		allGenes[geneName] = geneLength

	for sample in samples:
		capturedGenes = {} #{name:len}
		capturedGeneFiles = glob(os.path.join(data_dir, sample, "sequences", "FAA", "*.FAA"))
		for file in caturedGeneFiles:
			with open(file) as faa:
				name = file.split("/")[-1][:-4]
				faa.readline()
				length = len(faa.readline().strip())
				capturedGenes[name] = length

		coverage = sum(capturedGenes.values()) / sum(allGenes.values())
		if coverage < overallThreshold:
			badSamples = list(set(allGenes.keys()) - set(capturedGenes.keys()))
			for gene, length in capturedGenes.items():
				if length / allGenes[gene] < geneThreshold:
					badSamples.append(gene)

			samplesToRerun[sample] = badSamples
	return samplesToRerun

def run_hybpiper_selectively(data_dir, toRerun, main_script_dir):
	#toRerun is {sample, [genes]}
	for sample, genes in toRerun.items():
		with open(args.target_markers) as inFile, open(os.path.join(data_dir, "{}_filteredTargets.fas".format(sample)), "w") as outFile:
			for line in inFile:
				if line.strip().split("-")[1] in genes:
					outFile.write(line)
					outFile.write(inFile.readline())

		pairedFiles = glob(os.path.join(data_dir, "{}*trimmed_paired.fastq.gz".format(sample)))
		if len(pairedFiles) == 2:
			os.system("parallel -j 2 gunzip {} ::: " + " ".join(pairedFiles))
			sample_path = os.path.join(data_dir, sample + "_R*.trimmed_paired.fastq")
		else:
			singleFile = glob(os.path.join(data_dir, "{}*trimmed.fastq.gz".format(sample)))
			os.system("gunzip {}".format(singleFile[0])) #There will only be one
			sample_path = singleFile[0]

		run_Hybpiper =  '{}HybPiper/reads_first.py -b {} -r {}  --prefix {} --cpu {} '.format(main_script_dir, os.path.join(data_dir, sample+"_filteredTargets.fas"), sample_path, os.path.join(data_dir, sample+"_HybPiper"), args.cpu)
		logging.info("running HybPiper with: " + run_Hybpiper)
		os.system(run_Hybpiper)
		clean_command = "{}HybPiper/cleanup.py {}".format(main_script_dir, os.path.join(data_dir, sample+"_HybPiper"))
		os.system(clean_command)



def checkTestContinue(user_input):
	if user_input == "c" or user_input == "C" or user_input == "Continue" or user_input == "continue":
		return True

def check_arg():
	parser = argparse.ArgumentParser(description='UnFATE: the wrapper script that brings YOU from target enrichment sequencing data straight to phylogenetic tree inference! See the readme file for data structure and additional info.')
	mandatory_args = parser.add_argument_group("Mandatory or Suggested", "Arguments which are either required for the proper function of main_wrap.py, or should be used.")
	mandatory_args.add_argument('-c', '--cpu', default= '4', type=int,
						help='CPU number used by Hybpiper or parallel run of Exonerate, MACSE, RAxML etc.'
						)
	mandatory_args.add_argument('-o', '--out', required=True,
						help='The directory where output will be placed upon completion. MANDATORY ARGUMENT!'
						)
	mandatory_args.add_argument('-b', '--target_markers', default= 'Target_markers_rep_seq_aa.fas',
						help=' Path to fasta files containing all the sequences used to design the bait set, IT MUST BE A PROTEIN FASTA, USE  AN ABSOLUTE PATH!'
						)
	data_args = parser.add_argument_group("Data related", "Arguments relating to input data, at least one is required.")

	data_args.add_argument('-t', '--target_enrichment_data', default= '',
						help='Path to target enriched data. Files must end with "R1.fastq.gz", "R2.fastq.gz" or "SE.fastq.gz". Each sample can either consist of paired end reads or single end reads.',
						)
	data_args.add_argument('-w', '--whole_genome_data', default= '',
						help='Input path of de novo whole genome sequence data.'
						)
	data_args.add_argument('-a', '--assemblies', default= '',
						help='Path to assemblies. Files must end with ".fna.gz" or ".fna"',
						)
	optional_args = parser.add_argument_group("Optional", "Arguments which might be useful depending on your circumstances.")
	optional_args.add_argument('-f', '--first_use', action= 'store_true',
						help='Modifies some paths in MACSE pipeline folder, use this argument only if is the first time you run the pipeline, then do not move the UnFATE folder.',
						)
	optional_args.add_argument('-g', '--gblocks_relaxed', action= 'store_true',
						help='Applies Gblocks with relaxed parameters (Talavera & Castresana, 2007)',
						)
	optional_args.add_argument('-n', '--ncbi_assemblies', nargs = '+',
						help='Extracts the requested pre-mined NCBI assemblies genes from the database, takes a list of taxonomic ranks (e.g. Morchella Tuber Fuffaria). Include AUTO to have main_wrap choose which samples to add (increases run time much more than specifying a group with few samples). More info on AUTO can be found in the README.'
						)
	optional_args.add_argument('-x', '--test', action= 'store_true',
						help='Allows the user to exit early. Each step needs to be started by the user explicitly.'
						)
	optional_args.add_argument('-l', '--low_memory', action= 'store_true',
						help='Turns off SPAdes assembling of whole genome data before extracting sequences and uses HybPiper instead. Probably not required except on a desktop/laptop.'
						)
	optional_args.add_argument('-y', '--targ_hybpiper', action='store_true',
						help='Turns off metaSPAdes assembling of target enrichment data before extracting sequences and uses HybPiper instead. Not particularly recommended, as the metaSPAdes and exonerate_hits path has led to better results in our testing.'
						)
	optional_args.add_argument('-m', '--exonerate_mem',
						help='Limits the memory usage of exonerate when running on an assembly (-a or -w without -l). This does not strictly cap memory usage. More information in the README.',
						default=256,type=int
						)
	optional_args.add_argument('-e', '--threshold', default=55, type=int,
						help='Threshold for percent identity between contigs and proteins in exonerate_hits. default=55%%'
						)

	return parser.parse_args()

args = check_arg()
def main():
	logging.info("Arguments are: {}".format(args))
	testPrompt = "Finished {}, next step is {}: (C)ontinue or (Q)uit?"

	#print(args)
	global main_script_dir

	main_script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "")
	print(main_script_dir)
	#print(args.target_enrichment_data)
	#print(args.assemblies)
	os.system('ulimit -n 1024000')
	
	set_up_directories()	
	
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

	#path_to_merged_alignments = os.path.join(args.out, "alignments_merged/")
	path_to_macsed_align = os.path.join(args.out, 'macsed_alignments/')
	path_to_supermatrix= os.path.join(args.out, 'supermatrix/')
	if args.gblocks_relaxed:
		path_to_supermatrix_gblocked_dna = path_to_supermatrix +"supermatrix_gblocked_dna/"
		path_to_supermatrix_gblocked_aa = path_to_supermatrix + "supermatrix_gblocked_aa/"
	else:
		path_to_supermatrix_dna = path_to_supermatrix +"supermatrix_dna/"
		path_to_supermatrix_aa = path_to_supermatrix + "supermatrix_aa/"

	if os.path.isdir(os.path.join(args.out, "fastas")) and \
		os.path.isdir(os.path.join(args.out, "macsed_alignments")) and \
		os.path.isdir(os.path.join(args.out, "supermatrix")):

		print("fastas/, macsed_alignments/, supermatrix/ already exist, skipping steps.")
	else:
		if os.path.isdir(os.path.join(args.out, "fastas")):
			shutil.rmtree(os.path.join(args.out, "fastas"))
		if os.path.isdir(os.path.join(args.out, "macsed_alignments")):
			shutil.rmtree(os.path.join(args.out, "macsed_alignments"))
		if os.path.isdir(os.path.join(args.out, "supermatrix")):
			shutil.rmtree(os.path.join(args.out, "supermatrix"))

		all_genes = []
		#make files for all sequences
		with open(args.target_markers) as target_file:
			for line in target_file:
				if line.startswith(">"):
					all_genes.append(line.strip().split("-")[1])
		fastas_directory = os.path.join(args.out, "fastas", "")
		os.mkdir(fastas_directory)
		for gene in all_genes:
			fasta_file = "Alignment_{}_{}_merged.fasta"
			open(fastas_directory + fasta_file.format(gene, "protein"), "w").close()
			open(fastas_directory + fasta_file.format(gene, "nucleotide"), "w").close()

		ncbi_accessions = set()
		#we might want to change this to store accessions, instead of writing to a file
		if args.ncbi_assemblies:
			logging.info("Finding accessions for samples with specified taxonomy from database ")
			path_to_premined = main_script_dir + "combined_pre_mined_assemblies.tar.xz"

			if not os.path.isdir(main_script_dir + "combined_pre_mined_assemblies/"):
				unzip_premined_assemblies = "tar -C {} -Jxf {}".format(main_script_dir, path_to_premined)
				os.system(unzip_premined_assemblies)

			path_to_premined = main_script_dir + "combined_pre_mined_assemblies/"
			path_to_taxonomy = main_script_dir + "Accession_plus_taxonomy_Pezizomycotina.txt"
			for item in args.ncbi_assemblies:
				if item == "AUTO":
					continue
				# using the commas that are present in the csv file containing the taxonomy prevents to include ranks with similar name to the one requested in the command line (e.g. Fuffaria and Fuffarialongis)
				item1 = "," + item + ","
				with open(path_to_taxonomy, 'r') as taxonomy:
					for line in taxonomy:
						if item1 in line:
							ncbi_accessions.add(line.split(",")[0])
							#don't break because there are probably multiple with each
			print(ncbi_accessions)

		if args.whole_genome_data:
			logging.info("*****************************************************************************************************************************")
			logging.info("          TRIMMING WHOLE GENOME DATA FASTQ FILES  WITH TRIMMOMATC (Bolger et al. 2014)          ")
			logging.info("*****************************************************************************************************************************")
			trim_and_get_namelist(main_script_dir, args.whole_genome_data)
			namelist = 'namelist.txt'
			path_to_namelist = os.path.join(args.whole_genome_data, namelist)

			if args.low_memory:
				logging.info("Gunzipping paired reads trimmed fastq archives")
				gunzip_fastq = 'parallel -j {} gunzip ::: {}*_paired.fastq.gz'.format(args.cpu, args.whole_genome_data) 
				os.system(gunzip_fastq)
				gunzip_fastq = 'parallel -j {} gunzip ::: {}*trimmed.fastq.gz'.format(args.cpu, args.whole_genome_data)
				os.system(gunzip_fastq)

				logging.info("******************************************************************************************************************************************")
				logging.info("           EXTRACTING GENES FROM WHOLE GENOME DATA  WITH Hybpiper (Johnson et al. 2016)")
				logging.info("******************************************************************************************************************************************")
				run_hybpiper(main_script_dir, args.whole_genome_data, path_to_namelist)

				logging.info("Gzipping paired reads trimmed fastqs")
				gzip_fastq = "parallel -j {} gzip ::: {}*_paired.fastq".format(args.cpu, args.whole_genome_data)
				os.system(gzip_fastq)
				gzip_fastq = "parallel -j {} gzip ::: {}*trimmed.fastq".format(args.cpu, args.whole_genome_data)
				os.system(gzip_fastq)
			else:
				logging.info("*****************************************************************************************************************************")
				logging.info("          ASSEMBLING WHOLE GENOME DATA WITH SPADES (***CITATION NEEDED***)          ")
				logging.info("*****************************************************************************************************************************")
				with open(path_to_namelist) as namelistFile:
					for name in namelistFile:
						name = name.strip()
						if len(glob(os.path.join(args.out, "whole_genome_data", name + "*R*trimmed_paired.fastq.gz"))) == 2:
							sample_R1_path = os.path.join(args.out, "whole_genome_data", name + "_R1.trimmed_paired.fastq.gz")
							sample_R2_path = os.path.join(args.out, "whole_genome_data", name + "_R2.trimmed_paired.fastq.gz")
							spades_out_path = os.path.join(args.out, "whole_genome_data", name + "_spades/")
							spades_command = "spades.py -1 {} -2 {} -o {} -t {} --careful".format(sample_R1_path, sample_R2_path, spades_out_path, args.cpu)
							logging.info("running spades with " + spades_command)
							os.system(spades_command)
						else:
							sample_path = os.path.join(args.out, "whole_genome_data", name + "_SE.trimmed.fastq.gz")
							spades_out_path = os.path.join(args.out, "whole_genome_data", name + "_spades/")
							spades_command = "spades.py -s {} -o {} -t {} --careful".format(sample_path, spades_out_path, args.cpu)
							logging.info("running spades with " + spades_command)
							os.system(spades_command)

						os.rename(os.path.join(spades_out_path, "scaffolds.fasta"), os.path.join(args.whole_genome_data, name + ".fna"))

				#we have assemblies now
				logging.info("********************************************************************************************************************")
				logging.info("          PERFORMING ASSEMBLED WGS DATA ANALYSIS WITH Exonerate (Slater & Birney 2005)       ")
				logging.info("********************************************************************************************************************")
				run_exonerate(args.whole_genome_data)
				#toRerun = check_coverage(args.whole_genome_data) #returns {sample:[genes to rerun]}
				#run_hybpiper_selectively(args.whole_genome_data, toRerun, main_script_dir)

		#create hybpiper output in target_enrichment/
		if args.target_enrichment_data:
			logging.info("*****************************************************************************************************************************")
			logging.info("          TRIMMING TARGET ENRICHMENT FASTQ FILES  WITH TRIMMOMATC (Bolger et al. 2014)          ")
			logging.info("*****************************************************************************************************************************")
			logging.info('Path to TE data: '+args.target_enrichment_data)
			trim_and_get_namelist(main_script_dir, args.target_enrichment_data)
			namelist = 'namelist.txt'
			path_to_namelist = os.path.join(args.target_enrichment_data,namelist)

			if args.targ_hybpiper:
				#we only need to unzip if using HybPiper
				logging.info("Gunzipping paired reads trimmed fastq archives")
				gunzip_fastq =' parallel -j {} gunzip ::: {}*_paired.fastq.gz'.format(args.cpu, args.target_enrichment_data)
				os.system(gunzip_fastq)  # We dont' care about the unpaired reads if paired end are used
				gunzip_fastq = 'parallel -j {} gunzip ::: {}*trimmed.fastq.gz'.format(args.cpu, args.target_enrichment_data)
				os.system(gunzip_fastq)  # This covers single end reads

				#hybpiper
				logging.info("******************************************************************************************************************************************")
				logging.info("           EXTRACTING GENES FROM TARGET ENRICHMENT DATA  WITH Hybpiper (Johnson et al. 2016)")
				logging.info("******************************************************************************************************************************************")
				run_hybpiper(main_script_dir, args.target_enrichment_data, path_to_namelist)

				logging.info("Gzipping paired reads trimmed fastqs")
				gzip_fastq = "parallel -j {} gzip ::: {}*_paired.fastq".format(args.cpu, args.target_enrichment_data)
				os.system(gzip_fastq)
				gzip_fastq = "parallel -j {} gzip ::: {}*trimmed.fastq".format(args.cpu, args.target_enrichment_data)
				os.system(gzip_fastq)

			else:
				logging.info("*****************************************************************************************************************************")
				logging.info("          ASSEMBLING TARGET ENRICHMENT DATA WITH METASPADES (***CITATION NEEDED***)          ")
				logging.info("*****************************************************************************************************************************")
				with open(path_to_namelist) as namelistFile:
					for name in namelistFile:
						name = name.strip()
						if len(glob(os.path.join(args.out, "target_enrichment_data", name + "*R*trimmed_paired.fastq.gz"))) == 2:
							sample_R1_path = os.path.join(args.out, "target_enrichment_data", name + "_R1.trimmed_paired.fastq.gz")
							sample_R2_path = os.path.join(args.out, "target_enrichment_data", name + "_R2.trimmed_paired.fastq.gz")
							spades_out_path = os.path.join(args.out, "target_enrichment_data", name + "_spades/")
							spades_command = "spades.py -1 {} -2 {} -o {} -t {} --meta".format(sample_R1_path, sample_R2_path, spades_out_path, args.cpu)
							logging.info("running spades with " + spades_command)
							os.system(spades_command)
						else:
							sys.exit("metaspades does not work with single end reads yet. Please re-run with the -l (--low_memory) flag to use HybPiper instead of metaSPAdes for target enrichment data.")
							#metaspades does not accept single end reads yet
							#sample_path = os.path.join(args.out, "target_enrichment_data", name + "_SE.trimmed.fastq.gz")
							#spades_out_path = os.path.join(args.out, "target_enrichment_data", name + "_spades/")
							#spades_command = "spades.py -s {} -o {} -t {} --meta".format(sample_path, spades_out_path, args.cpu)
							#logging.info("running spades with " + spades_command)
							#os.system(spades_command)

						os.rename(os.path.join(spades_out_path, "scaffolds.fasta"), os.path.join(args.target_enrichment_data, name + ".fna"))

				#we have assemblies now
				logging.info("********************************************************************************************************************")
				logging.info("          PERFORMING ASSEMBLED TARGET ENRICHMENT DATA ANALYSIS WITH Exonerate (Slater & Birney 2005)       ")
				logging.info("********************************************************************************************************************")
				run_exonerate(args.target_enrichment_data)


		#user input: assemblies
		#create exonerate output in assemblies/
		if args.assemblies:
			logging.info("********************************************************************************************************************")
			logging.info("          PERFORMING ASSEMBLIES DATA ANALYSIS WITH Exonerate (Slater & Birney 2005)       ")
			logging.info("********************************************************************************************************************")
			run_exonerate(args.assemblies)

		#start filling fastas/
		logging.info("*********************************************")
		logging.info("          BUILDING FASTA FILES          ")
		logging.info("*********************************************")

		user_samples = []

		if args.assemblies:
			logging.info("Adding sequences from assemblies data")
			#get_alignment(args.assemblies)
			get_fastas_exonerate(args.assemblies, True)
			user_samples.extend(get_names(args.assemblies, True))
		if args.target_enrichment_data:
			logging.info("Adding sequences from target enrichment data")
			if args.targ_hybpiper:	#done with HybPiper
				get_fastas_exonerate(args.target_enrichment_data, False)
				user_samples.extend(get_names(args.target_enrichment_data, False))
			else:			#done with metaspades then exonerate
				get_fastas_exonerate(args.target_enrichment_data, True)
				user_samples.extend(get_names(args.target_enrichment_data, True))
		if args.whole_genome_data:
			logging.info("Adding sequences from whole genome data")
			if args.low_memory:	#done with hybpiper
				get_fastas_exonerate(args.whole_genome_data, False)
				user_samples.extend(get_names(args.whole_genome_data, False))

			else:			#done with spades then exonerate
				get_fastas_exonerate(args.whole_genome_data, True)
				user_samples.extend(get_names(args.whole_genome_data, True))

		if args.ncbi_assemblies:
			logging.info("****************************************************************************************************************************")
			logging.info("           ADDING SELECTED TAXONOMIC RANKS GENES FROM PRE-MINED ASSEMBLY DATABASE           ")
			logging.info("****************************************************************************************************************************")

			#remove requirement for AUTO to be exclusive: add other stuff anyway, but make sure there aren't duplicates
			if "AUTO" in args.ncbi_assemblies:
				logging.info("Using mafft to align your sequences with sequences from the database")
				#align new stuff in fastas/ with pre-aligned stuff
				#make supermatrix
				#find scores (discussion on this elsewhere)
				#add stuff we want to fastas (if they can have the same name as if nothing happened, we don't need to change anything below)
				auto_dir = os.path.join(args.out, "auto_selection")
				if not os.path.isdir(auto_dir):
					os.mkdir(auto_dir)
				for file in glob(os.path.join(args.out, "fastas", "*nucleotide*.fasta")):
					geneName = file.split("/")[-1].split("_")[1] #fastas/abc_123at4980_protein_merged.fasta
					mafft_command = "mafft --retree 1 --thread {} --add {} {} > {}"
					prealignedFile = os.path.join(main_script_dir, "msa_test", "truncatedIDs_{}.fasta".format(geneName))
					newFile = os.path.join(auto_dir, "added_{}.fasta".format(geneName))
					os.system(mafft_command.format(args.cpu, file, prealignedFile, newFile))

				os.chdir(auto_dir)
				logging.info("Concatenating alignments of user samples and db samples")
				os.system("perl {}/FASconCAT-G_v1.04.pl -i -s".format(main_script_dir))
				os.chdir(main_script_dir)

				print(user_samples)
				found_user_samples = set()
				#for id that isn't in found_user_samples, find closest things
				#if id in user_samples is close, add to found_user_samples
				#get ids from stuff in db to add, then add with SeqIO
				for id in user_samples:
					if id in found_user_samples:
						continue
					similar_samples = find_similar_samples(id, user_samples, auto_dir, args.cpu)
					for sample in similar_samples:
						if sample in user_samples:
							found_user_samples.add(sample)
						else:
							ncbi_accessions.add(sample)
				print(ncbi_accessions)

				#JACOB
				#sys.exit("What more would you want?")


				#Add sequences from database to sequences from supplied data
				path_to_premined_combined = os.path.join(main_script_dir, "combined_pre_mined_assemblies")
				for fasta in glob(os.path.join(args.out, "fastas", "*.fasta")):
					#logging.info("NCBI to fasta FASTA: " + fasta)
					baseName = fasta.split("/")[-1].split("_")[1]
					moleculeType = fasta.split("/")[-1].split("_")[2]
					#logging.info("NCBI to fasta MOLECULETYPE: " + moleculeType)
					if moleculeType == "protein":
						with open(os.path.join(path_to_premined_combined, "combined_" + baseName + ".FAA")) as proteinIn, open(fasta, 'a') as proteinOut:
							for line in proteinIn:
								if line.startswith(">"):
									for accession in ncbi_accessions:
										if accession in line:
											proteinOut.write(line)
											sequenceLine = proteinIn.readline()
											proteinOut.write(sequenceLine)

					if moleculeType == "nucleotide":
						with open(os.path.join(path_to_premined_combined, "combined_" + baseName + ".FNA")) as nucIn, open(fasta, 'a') as nucOut:
							for line in nucIn:
								if line.startswith(">"):
									for accession in ncbi_accessions:
										if accession in line:
											nucOut.write(line)
											sequenceLine = nucIn.readline()
											nucOut.write(sequenceLine)


		# Remove empty fastas, just in case there are no sequences for a specific gene
		os.system("find {} -type f -empty -delete".format(os.path.join(fastas_directory, "*_merged.fasta")))

		# Get rid of the trash strings after the accession number to be able to replace with speciens name later
		# As OMM_MACSE will use soft masking to align and trim better get rig of all small case letters in the alignments before running MACSE pipeline
		logging.info("Cleaning sequences names to only retain accession numbers...")
		logging.info("Converting all nucleotides to uppercase...")

		#path_to_merged_alignments = os.path.join(args.out, "fastas", "")

		for f in glob(os.path.join(fastas_directory, "*_merged.fasta")):
			output_file = open(f.rstrip("\.fasta") + "_headmod.fas","a")
			for seq in SeqIO.parse(f,"fasta"):
				regex_id = re.search("(^GCA_[0-9]+.[0-9])_", seq.id)
				if regex_id is not None:
					# this strips the old header out (.id is only the accession in theory .description is the whole header instead
					seq.description = ""
					seq.id = regex_id.group(1)
					# replace sequence with the same but uppercase
					sequence = str(seq.seq).upper()
					# arrange sequence and id in a format that SeqIO can write to file
					record = SeqRecord(Seq(sequence), seq.id, "","")
					#print(record)
					SeqIO.write(record, output_file,"fasta")
				else:
					sequence = str(seq.seq).upper()
					record = SeqRecord(Seq(sequence), seq.id, "","")
					#print(record)
					SeqIO.write(record, output_file,"fasta")
			output_file.close()

		logging.info("**********************************************************************************************************")
		logging.info("          COMPARING RETRIEVED GENES TO REFERENCE SEQUENCES LENGTH          ")
		logging.info("**********************************************************************************************************")
		markers_retrieved_percentage_script=main_script_dir + "markers_retrieved_percentage.py"
		run_markers_retrieved_percentage = "python3 {} -b {} -f {} ".format(markers_retrieved_percentage_script, args.target_markers, fastas_directory)	
		os.system(run_markers_retrieved_percentage)

		#end commented out if

		logging.info("*************************************************************************************************************************************")
		logging.info("          PERFORMING ALIGNMENT WITH OMM_MACSE (Ranwez et al. 2018; Di Franco et al. 2019)           ")
		logging.info("*************************************************************************************************************************************")

		os.chdir(fastas_directory)

		MACSE_dir = main_script_dir + "MACSE_V2_PIPELINES/OMM_MACSE/"
		MACSE_script = MACSE_dir + "S_OMM_MACSE_V10.02.sh"
		# As OMM_MACSE uses soft masking put all the sequences in upper case before alignment and filtering 
		run_OMM_MACSE = 'find %s -type f -name "*_nucleotide_merged_headmod.fas" | parallel -j %s %s --out_dir {}_out --out_file_prefix macsed --in_seq_file {} --no_prefiltering --no_postfiltering --alignAA_soft MAFFT  --min_percent_NT_at_ends 0.01 ' %(fastas_directory, args.cpu, MACSE_script)
		os.system(run_OMM_MACSE)
		logging.info(run_OMM_MACSE)
		# move aligned files
		path_to_macsed_align = os.path.join(args.out, "macsed_alignments", "")
		make_align_fold = "mkdir {}".format(path_to_macsed_align)
		os.system(make_align_fold)
		for root, dirs, files in os.walk(fastas_directory, topdown=True):
			for f in files:
				if f.endswith("_final_align_NT.aln") or f.endswith("_final_align_AA.aln"):
					#print(root)
					#print(f)
					file_path = root +"/"+ f   
					regex1 =re.search("Alignment_(\S+)_nucleotide_merged_headmod.fas_out",root)
					#os.rename both renames and moves files
					os.rename(file_path,  path_to_macsed_align + regex1.group(1) + f + ".fas")

		gblocks_path= main_script_dir + "Gblocks"
		if args.gblocks_relaxed:
			logging.info("***********************************************************************************************************")
			logging.info("          PERFORMING ALIGNMENT FILTERING WITH Gblocks (Castresana, 2000)            ")
			logging.info("***********************************************************************************************************")
			#print(path_to_macsed_align)
			run_gblocks("_final_align_NT.aln.fas","_final_align_AA.aln.fas", path_to_macsed_align, gblocks_path)

		logging.info("Deleting spaces and sequences with > 70% gaps from alignments")
		for f in os.listdir(path_to_macsed_align):
			if f.endswith("-gb") or f.endswith(".fas"):
				logging.info("Processing %s" % f)
				alignment_file = open(path_to_macsed_align + f,'rt')
				alignment_file_content = alignment_file.read()
				#delete spaces created by Gblocks in the alignments, if it is not a Gblocked alignment it does the same just to be safe (so there will be no spaces at all in the alignment file)
				alignment_file_content = alignment_file_content.replace(' ','')
				alignment_file.close()
				alignment_file = open(path_to_macsed_align + f,'wt')
				alignment_file.write(alignment_file_content)
				alignment_file.close()

		# get rid of the sequences with > 75% of gaps (sequence which remained empty o almost empty after filtering, i.e. not well alignable and/or already short before filtering)
		for f in os.listdir(path_to_macsed_align):
			if f.endswith("-gb") or f.endswith(".fas"):
				logging.info("Processing %s" % f)
				input_file = path_to_macsed_align + f
				with open(path_to_macsed_align + f +"_cleaned.fasta", "a+") as output_file:
					for record in SeqIO.parse(input_file, 'fasta'):
						drop_cutoff = 0.75
						name = record.id
						seq = record.seq
						seqLen = len(seq)
						gap_count = seq.count("-")
						#print(seqLen)
						#print(gap_count)
						#print(gap_count/(seqLen))
						if seqLen == 0 or (gap_count/seqLen) >= drop_cutoff:
							logging.info(" %s was removed." % name)
						else:
							SeqIO.write(record, output_file , 'fasta')

		logging.info("************************************************************************************************************************")
		logging.info("          PERFORMING ALIGNMENTS CONCATENATION WITH Fasconcat (Kück & Longo, 2014)          ")
		logging.info("*************************************************************************************************************************")
		#path_to_supermatrix= path_to_macsed_align.replace('macsed_alignments/', 'supermatrix/')
		make_supermatrix_folder="mkdir {} ".format(path_to_supermatrix)
		os.system(make_supermatrix_folder)

		if args.gblocks_relaxed:
			os.system("mkdir {}".format(path_to_supermatrix_gblocked_dna))
			os.system("mkdir {}".format(path_to_supermatrix_gblocked_aa))
		else:
			os.system("mkdir {}".format(path_to_supermatrix_dna))
			os.system("mkdir {}".format(path_to_supermatrix_aa))

		if args.gblocks_relaxed:
			os.system("cp -r {}*macsed_final_align_NT.aln.fas-gb_cleaned.fasta {}".format(path_to_macsed_align, path_to_supermatrix_gblocked_dna))
			os.system("cp -r {}*macsed_final_align_AA.aln.fas-gb_cleaned.fasta {}".format(path_to_macsed_align, path_to_supermatrix_gblocked_aa))
		else:
			os.system("cp -r {}*macsed_final_align_NT.aln.fas_cleaned.fasta {}".format(path_to_macsed_align, path_to_supermatrix_dna))
			os.system("cp -r {}*macsed_final_align_AA.aln.fas_cleaned.fasta {}".format(path_to_macsed_align, path_to_supermatrix_aa))

		FASconCAT_command = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(main_script_dir)
		for moleculeType in ["NT", "AA"]:
			if args.gblocks_relaxed:
				if moleculeType == "NT":
					os.chdir(path_to_supermatrix_gblocked_dna)
				else:
					os.chdir(path_to_supermatrix_gblocked_aa)

				# Remove empty fastas, just in case there are no sequences for a specific gene
				os.system("find *.fasta -type f -empty -delete")

				os.system(FASconCAT_command)
				os.rename('FcC_supermatrix.fas','FcC_supermatrix_gblocked_{}.fasta'.format(moleculeType))
				os.rename('FcC_supermatrix_partition.txt','FcC_supermatrix_partition_gblocked_{}.txt'.format(moleculeType))
				os.rename('FcC_info.xls','FcC_info_gblocked_{}.xls'.format(moleculeType))
				os.system("rm *_cleaned.fasta")
			else:
				if moleculeType == "NT":
					os.chdir(path_to_supermatrix_dna)
				else:
					os.chdir(path_to_supermatrix_aa)

				# Remove empty fastas, just in case there are no sequences for a specific gene
				os.system("find *.fasta -type f -empty -delete")

				os.system(FASconCAT_command)
				os.rename('FcC_supermatrix.fas','FcC_supermatrix_{}.fasta'.format(moleculeType))
				os.rename('FcC_supermatrix_partition.txt','FcC_supermatrix_partition_{}.txt'.format(moleculeType))
				os.rename('FcC_info.xls','FcC_info_{}.xls'.format(moleculeType))
				os.system("rm *_cleaned.fasta")

		os.chdir(main_script_dir)

	if args.test:
		test_input = input(testPrompt.format("alignment and concatenation", "single tree reconstruction"))
		if not checkTestContinue(test_input):
			sys.exit("exited after concatenating alignments")


	logging.info("******************************************************************************************************************")
	logging.info("          RECONSTRUCTING SINGLE MARKER TREES WITH IQTREE2 (Minh et al. 2020)           ")
	logging.info("******************************************************************************************************************")
	iqtree_script = main_script_dir + "iqtree2"
	#print(path_to_macsed_align)
	# DNA alignments
	if args.gblocks_relaxed:
		iqtree_parallel1 = "find %s -type f  -name '*_final_align_NT.aln.fas-gb_cleaned.fasta' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
		os.system(iqtree_parallel1)
	else:
		iqtree_parallel = "find %s -type f  -name '*_final_align_NT.aln.fas_cleaned.fasta' | parallel -j %s %s -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
		os.system(iqtree_parallel)
	# Amino acid alignments
	if args.gblocks_relaxed:
		iqtree_parallel3 = "find %s -type f  -name '*_final_align_AA.aln.fas-gb_cleaned.fasta' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
		os.system(iqtree_parallel3)
	else:
		iqtree_parallel2 = "find %s -type f  -name '*_final_align_AA.aln.fas_cleaned.fasta' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
		os.system(iqtree_parallel2)
	# move trees and other iqtree files to the dedicated folder
	path_to_single_trees = os.path.join(args.out, 'single_locus_trees/')
	mkdir_raxml_out = "mkdir {}".format(path_to_single_trees)
	os.system(mkdir_raxml_out)
	for root, dirs, files in os.walk(path_to_macsed_align, topdown=True):
		for f in files:
			if  f.endswith("treefile") or f.endswith("nex") or f.endswith("parttrees") or f.endswith("gz") or f.endswith("mldist") or f.endswith("log") or f.endswith("iqtree") or f.endswith("contree") or f.endswith("bionj") or f.endswith("best_scheme"):
				os.rename(path_to_macsed_align + f, path_to_single_trees + f)

	if args.test:
		test_input = input(testPrompt.format("single trees", "supermatrix tree"))
		if not checkTestContinue(test_input):
			sys.exit("exited after making single marker trees")



	logging.info("*************************************************************************************************************")
	logging.info("          RECONSTRUCTING SUPERMATRIX TREE WITH IQTREE2  (Minh et al. 2020)          ")
	logging.info("*************************************************************************************************************")
	iqtree_script=main_script_dir + "iqtree2"
	if args.gblocks_relaxed:
		iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s --threads-max %s" %(iqtree_script, path_to_supermatrix_gblocked_dna + 'FcC_supermatrix_gblocked_NT.fasta' , path_to_supermatrix_gblocked_dna + 'FcC_supermatrix_partition_gblocked_NT.txt', "AUTO", args.cpu)
		os.system(iqtree_on_supermatrix)
		iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s --threads-max %s" %(iqtree_script, path_to_supermatrix_gblocked_aa + 'FcC_supermatrix_gblocked_AA.fasta' , path_to_supermatrix_gblocked_aa + 'FcC_supermatrix_partition_gblocked_AA.txt', "AUTO", args.cpu)
		os.system(iqtree_on_supermatrix)
	else:
		iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s --threads-max %s" %(iqtree_script, path_to_supermatrix_dna + 'FcC_supermatrix_NT.fasta' , path_to_supermatrix_dna + 'FcC_supermatrix_partition_NT.txt', "AUTO", args.cpu)
		os.system(iqtree_on_supermatrix)
		iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s --threads-max %s" %(iqtree_script, path_to_supermatrix_aa + 'FcC_supermatrix_AA.fasta' , path_to_supermatrix_aa + 'FcC_supermatrix_partition_AA.txt', "AUTO", args.cpu)
		os.system(iqtree_on_supermatrix)

	if args.test: 
		test_input = input(testPrompt.format("supermatrix iqtree2", "astral tree"))
		if not checkTestContinue(test_input):
			sys.exit("exited after reconstructing phylogeny from supermatrix, output may be in interesting place.")

	##### MAYBE RUN RAXML-NG USING THE IQTREE TOPOLOGY TO GET mbe BOOTSTRAP VALUES (very slow, as it is a complete boostrap?
	
	logging.info("******************************************************************************************************")
	logging.info("            RECONSTRUCTING SUPERTREE WITH ASTRAL (Zhang et al. 2018)")
	logging.info("******************************************************************************************************")
	path_to_supertree = os.path.join(args.out,'supertree/')
	os.system("mkdir {}".format(path_to_supertree))

	if args.gblocks_relaxed:
		gblocked_text = "gblocked_"
	else:
		gblocked_text = "" #don't add "gblocked_" to output files or directories
	
	dna_dir = path_to_supertree + "supertree_{}dna/".format(gblocked_text)
	aa_dir = path_to_supertree + "supertree_{}aa/".format(gblocked_text)
	os.system("mkdir {}".format(dna_dir))
	os.system("mkdir {}".format(aa_dir))
	os.system("cp {}*NT*.treefile {}".format(path_to_single_trees, dna_dir))
	os.system("cp {}*AA*.treefile {}".format(path_to_single_trees, aa_dir))
	os.chdir(dna_dir)
	os.system("cat *.treefile > cat_trees_{}dna.tre".format(gblocked_text))
	os.system("rm *.treefile")
	os.chdir(aa_dir)
	os.system("cat *.treefile > cat_trees_{}aa.tre".format(gblocked_text))
	os.system("rm *.treefile")
	os.chdir(main_script_dir)
	path_to_astral = main_script_dir + "ASTRAL/astral.5.7.7.jar"
	run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, dna_dir + "cat_trees_{}dna.tre".format(gblocked_text), dna_dir + "astral_species_tree_{}dna.tree".format(gblocked_text))
	os.system(run_astral)
	run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, aa_dir + "cat_trees_{}aa.tre".format(gblocked_text), aa_dir + "astral_species_tree_{}aa.tree".format(gblocked_text))
	os.system(run_astral)
			
	logging.info("****************************************************************")
	logging.info("          COPYING TREES TO final_trees FOLDER           ")
	logging.info("****************************************************************")
	# Copy the tree file from IQTREE and ASTRAL to a new folder
	path_to_finaltrees = os.path.join(args.out, 'final_trees/')
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
	logging.info("          CONVERTING ACCESSIONS IN THE TREES, if any, TO SPECIES NAME")
	logging.info("*************************************************************************************************")
	# Open one of the suprematrices, making a list of accession sample name
	if args.gblocks_relaxed:
		supermatrix_file = path_to_supermatrix_gblocked_dna + 'FcC_supermatrix_gblocked_NT.fasta'
	else:
		supermatrix_file = path_to_supermatrix_dna + 'FcC_supermatrix_NT.fasta'
	supermatrix_accession_file = path_to_finaltrees + 'Accessions_not_found.csv'
	accessions_plus_taxonomy_file = path_to_finaltrees + 'Accessions_plus_taxonomy.csv'
	with open(supermatrix_accession_file, 'w') as accessions, open(accessions_plus_taxonomy_file, 'w') as accessions_tax,\
			open(supermatrix_file, 'r') as supermatrix, open(os.path.join(main_script_dir, "Accession_plus_taxonomy_Pezizomycotina.txt")) as tax_in:
		supermatrix_content = supermatrix.readlines()
		all_accessions = []
		accessions_added = []
		for line in supermatrix_content:
			regex = re.search("^>(GCA_[0-9]+\.[0-9])", line)
			if regex:
				all_accessions.append(regex.group(1))
#				accessions.write(regex.group(1)+"\n")	
			else:
				pass
		for line in tax_in:
			if line.split(",")[0] in all_accessions:
				accessions_tax.write(line)
				accessions_added.append(line.split(",")[0])
			else:
				pass
		#by this point, all previously known taxonomies have been added
		for id in accessions_added:
			all_accessions.remove(id)
		accessions.write("\n".join(all_accessions))

	# Add taxonomy to the accessions retrieved (get_taxonomy_with edirect script), select species name  and format the .csv file 
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
