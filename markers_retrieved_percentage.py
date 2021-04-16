import os 
import click
from Bio import SeqIO
import sys
import pandas
import re
import seaborn
import matplotlib.pyplot as plt
import logging

@click.command()
@click.option('--bait_file_aa_path','-b', default='Hybpiper_targets_first10markers_aa.fas', help='path to bait file')
@click.option('--alignments_folder_path','-f', default='./', help='path of the folder containing the fastq.gz file to trim with trimmomatic')

def seq_percentage(bait_file_aa_path, alignments_folder_path,  ):
	""" Takes the protein fasta (the reference sequences) and the alignment files.
		 For each marker an average length values for the reference sequences is produced.
		 The reference length is compared to the sequences retrieved for each sample.   
	"""
	gene_names = []
	reference_lengths = {}
	## Generate a dictionary whit gene names associated with average length of the reference sequences
	# read fasta, split the header (.id in SeqIO) using "-" as delimiter and take the last element [-1]
	for prot in SeqIO.parse(bait_file_aa_path,"fasta"):
		protname = prot.id.split("-")[-1]
		gene_names.append(protname)
		# if the gene is in the list append the length of that seq  (.seq in SeqIO) to the dictionary
		if protname in reference_lengths:
			reference_lengths[protname].append(len(prot.seq))
		# else add the gene name and the length in a new voice of the dictionary. Output is a dictionary with gene name and list of legth for that gene
		else:
			reference_lengths[protname] = [len(prot.seq)]
	unique_names = list(set(gene_names))
	#print(unique_names)
	# for every gene make the average length (sum lengths and divide for length list length), put all in a list
	avg_ref_lengths = [int(sum(reference_lengths[gene])/len(reference_lengths[gene])) for gene in unique_names]
	#print(avg_ref_lengths)
	# create a zip interator (iterate over two values in couple) and then make a dictionary out of it. Now each reference gene is associated to the average length
	zip_iterator = zip(unique_names, avg_ref_lengths)
	avg_ref_len_dict = {}
	avg_ref_len_dict = dict(zip_iterator)
	logging.info("Dictionary of genes and their reference sequences average length:")
	logging.info(avg_ref_len_dict)
	
	# Makes a list of the samples: append all the header to a list that eliminate redundancy doing a set (not efficient but works)
	sample_list = []
	for f in os.listdir(alignments_folder_path):
		for protein in SeqIO.parse(alignments_folder_path + f ,"fasta"):
			sample_list.append(protein.id)
	sample_list = set(sample_list)				
	#print(sample_list)
	
	## Makes a table that associate samples (Genome Accessions) to gene length (pandas dataframe??)
	# a dataframe is created: columns are the gene names, rows (index) are the samples
	data_frame = pandas.DataFrame(columns = unique_names, index = sample_list )
	#print(data_frame)
	
	# iterate over the alignments
	for f in os.listdir(alignments_folder_path):
		if f.endswith("_protein_merged.fasta"):
			sample_gene_length_dict = {}
			regex = re.search("Alignment_([0-9]+at4890)_protein_merged.fasta",f)
			# every alignment is converted to a dictionary containing the sample name and its sequence length value
			for protein in SeqIO.parse(alignments_folder_path + f ,"fasta"):
				sample_gene_length_dict[protein.id] = len(protein.seq)
				# dataframe is filled column by column (regex.group(1) is the gene name)
				data_frame[regex.group(1)] = pandas.Series(sample_gene_length_dict)
				# add a row (.loc is for that purpose, otherwise is a column) with reference sequences average length
				data_frame.loc['REFERENCE_AVG'] = pandas.Series(avg_ref_len_dict)
			#print(sample_gene_lenght_dict)
		else:
			pass
		
			#print(protein.id)
			#print(protein.seq)
	logging.info("Exporting the length table table to gene_length.csv. LAST ROW REPRESENTS THE AVERAGE VALUE OF THE REFERECE SEQUENCES USED TO EXTRACT THESE GENES ")
	logging.info(data_frame)
	# Exports dataframe to .csv file
	data_frame.to_csv(alignments_folder_path + 'gene_length.csv', encoding='utf-8')
			
	# Plot an heatmap from dataframe table
	fig, ax = plt.subplots(figsize=(60, 60))
	seaborn.heatmap(data_frame, cmap="RdBu")
	plt.show()
	
# starts the function					
if __name__ == '__main__':
	seq_percentage()	
