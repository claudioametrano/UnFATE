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
@click.option('--bait_file_aa_path','-b', default='./all_genes_names_FG_2_Hybpiper_format_aa.fas', help='path to baits file')
@click.option('--alignments_folder_path','-f', default='.', help='path of the folder containing the fasta files')
@click.option('--plot_heatmap','-p', is_flag=True, help=' if used, plots the heatmap using the dataframe (WARNING: do not use on a server, it need a GUI)')

def seq_percentage(bait_file_aa_path, alignments_folder_path, plot_heatmap ):
	""" Takes the protein fasta (the reference sequences) and the alignment files.
		 For each marker an average length values for the reference sequences is produced.
		 The reference length is compared to the sequences retrieved for each sample.   
	"""
	logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
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
		if f.endswith("_protein_merged.fasta"):
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
				data_frame.loc['REFERENCE_AVG_LENGTH'] = pandas.Series(avg_ref_len_dict)
			#print(sample_gene_lenght_dict)
		else:
			pass

			#print(protein.id)
			#print(protein.seq)

	#normalized data frame
	ndf = data_frame.copy()

	#get the row with mean marker length, then divide everything else by that, columnwise
	marker_means = ndf.iloc[-1:,:]
	ndf = ndf.iloc[:,:].div(marker_means.iloc[0,:], axis='columns')

	logging.info(ndf)
	logging.info(ndf.shape)

	logging.info("Exporting the length table to gene_length.csv in 'alignments_merged' folder")
	logging.info("LAST ROW OF THE TABLE REPRESENTS THE AVERAGE VALUE OF THE REFERECE SEQUENCES USED TO EXTRACT THESE GENES")
	logging.info(data_frame)
	# Exports dataframe to .csv file
	data_frame.to_csv(alignments_folder_path + 'gene_length.csv', encoding='utf-8')
	ndf.to_csv(alignments_folder_path + 'gene_length_normalized.csv', encoding='utf-8')
			
	# Plot an heatmap from dataframe table
	fig, ax = plt.subplots(figsize=(100, 60))
	seaborn.heatmap(data_frame, cmap="Greens", vmin=0, xticklabels=1, yticklabels=1)
	if plot_heatmap:
		logging.info("Plotting the heatmap...")
		plt.show()
	else:	
		logging.info("Exporting the length table as heatmap to 'gene_lengths_heatmap.pdf' in the 'fastas' folder")
		plt.savefig(alignments_folder_path + 'gene_lengths_heatmap.pdf')

	fig, ax = plt.subplots(figsize=(100,60))
	seaborn.heatmap(ndf, cmap="Greens", vmin=0, xticklabels=1, yticklabels=1)
	if plot_heatmap:
		logging.info("Plotting the heatmap...")
		plt.show()
	else:
		logging.info("Exporting the normalized length table as heatmap to 'gene_lengths_normalized_heatmap.pdf' in the 'fastas' folder")
		plt.savefig(alignments_folder_path + 'gene_lengths_normalized_heatmap.pdf')
	
		
# starts the function					
if __name__ == '__main__':
	seq_percentage()	
