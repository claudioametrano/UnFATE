#!/usr/bin/env python3
import shutil
import re
import os
import sys
import argparse
import multiprocessing
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
        #blast_columns not needed anymore, but it's nice to have them for reference.
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
    """in the assembly mode that uses exonerate_hits.py script from Hybpiper, there is no way to know what sequences from the reference protein sequecne is the best,
        as we do not have the BLAST or the BWA results mappig the reads from target enrichment to those sequences, so we just BLAST all the reference for each gene
        on the assemblies, the one with best bitscore is used to run exonerate_hits.py for that gene.
    """
    ref_gene_list = []
    # Get a list of genes without redundancy
    for seq in SeqIO.parse(prot_file_path,"fasta"):
        geneName = seq.id.strip().split("-")[1]
        ref_gene_list.append(geneName)
    ref_gene_list = list(set(ref_gene_list))
    logging.info("Gene list: ")
    logging.info(ref_gene_list)
    #     Generate a separate fasta for each gene in the reference sequnces file
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
    extensions = ["*_ref.fasta", "*.nin", "*.nsq", "*.nhr", "*.ndb", "*.nto", "*.not", "*.ntf", "*.njs"]
    for extension in extensions:
        remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
        os.system(remove_lot_of_files)
    return()

def get_names(path_to_data, isAssemblies):
    if isAssemblies:
        search_location = os.path.join(path_to_data, "*", "sequences")
    else:
        search_location = os.path.join(path_to_data, "*", "exonerate_genelist.txt")
    found_samples = set()
    for result in glob(search_location):
        #print("get_names result: ",result)
        found_samples.add(result.split("/")[-2])
    return list(found_samples)

def get_fastas_exonerate(path_to_data, isAssemblies):
    for moleculeType in ["FNA", "FAA"]:
        if isAssemblies:
            search_location = os.path.join(path_to_data, "*", "*", "*", "sequences", moleculeType, "*")
        else:
            search_location = os.path.join(path_to_data, "*", "*", "*", "sequences", moleculeType, "*")
        for filename in glob(search_location):
            geneName = filename.split("/")[-1][:-4]
            #print("FILENAME: ",filename, "GENENAME: ", geneName)
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
        logging.info("Processed {} of {}".format(index, total))
    return((index, score))

def find_similar_samples(query, user_samples, data_dir, cpus):
    logging.info("Finding similar samples to {} in Pezizomycotina database".format(query))
    aln = AlignIO.read(open(os.path.join(data_dir, "FcC_supermatrix.fas")), "fasta")
    for record in aln:
        record.seq = record.seq.upper()

    #distances will be calculated with a basic 2-parameter model
    calc = DistanceCalculator("trans", skip_letters=["N", "-", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V"])

    queryIndex = [record.id for record in aln].index(query)
    sequences = [record.seq for record in aln]
    ids = [record.id for record in aln]
    querySeq = sequences[queryIndex]

    scores = [] # will be in form of: [(index, score), ...]

    with multiprocessing.Pool(cpus) as p:
        scores=p.starmap(get_score, [(querySeq, seq, i, calc, len(sequences)) for i,seq in enumerate(sequences)])
    scores.sort(key = lambda x: x[1])

    idScores = [(ids[elements[0]], elements[1]) for elements in scores]
    return select_scores(idScores, user_samples)

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
                logging.info("Not found: ",l[0])
                tree_content = tree_content
        logging.info(f"Substitutions done:  {count}")        
        return(tree_content)

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
                    logging.info("Path already exists")
                    continue #symlink already exists
                if file.lower().endswith(".fastq") or file.lower().endswith(".fastq.gz"):
                    os.symlink(os.path.join(args.target_enrichment_data, file), os.path.join(args.out, "target_enrichment_data", file))
            args.target_enrichment_data = os.path.join(args.out, "target_enrichment_data", "") #the empty one causes a trailing /

        if args.assemblies:
            for file in os.listdir(args.assemblies):
                if os.path.exists(os.path.join(args.out, "assemblies", file)):
                    logging.info("Path already exists")
                    continue
                if file.lower().endswith(".fna") or file.lower().endswith(".fasta") or file.lower().endswith(".fna.gz"):
                    os.symlink(os.path.join(args.assemblies, file), os.path.join(args.out, "assemblies", file))
            args.assemblies = os.path.join(args.out, "assemblies", "")

        if args.whole_genome_data:
            for file in os.listdir(args.whole_genome_data):
                if os.path.exists(os.path.join(args.out, "whole_genome_data", file)):
                    logging.info("Path already exists")
                    continue #symlink already exists
                if file.lower().endswith(".fastq") or file.lower().endswith(".fastq.gz"):
                    os.symlink(os.path.join(args.whole_genome_data, file), os.path.join(args.out, "whole_genome_data", file))
            args.whole_genome_data = os.path.join(args.out, "whole_genome_data", "") #the empty one causes a trailing /

def trim_and_get_namelist(exec_dir, data_dir):
    trimming_cmd = "python3 {}/trimmer.py -f {} -c {}".format(exec_dir, data_dir, args.cpu)
    os.system(trimming_cmd)
    #Get namelist.txt from dataset directory
    namelist_cmd = 'python3 {}/getNameList.py -f {}'.format(exec_dir, data_dir)
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
            run_Hybpiper =  'hybpiper assemble -t_aa {} -r {}  --prefix {} --cpu {}'.format(args.target_markers, sample_path, line.strip(), args.cpu)    
            # old command for hybpiper1
            #run_Hybpiper =  '{}HybPiper/reads_first.py -b {} -r {}  --prefix {} --cpu {} '.format(main_script_dir, args.target_markers, sample_path, line.strip(), args.cpu)
            logging.info("running HybPiper with: " + run_Hybpiper)
            os.system(run_Hybpiper)
            #new hybpiper cleans up by defaulf
            #clean_command = "{}HybPiper/cleanup.py {}".format(main_script_dir, line.strip())
            #os.system(clean_command)
    os.chdir(main_script_dir)

#def run_exonerate_hits(file_, ref_seq_file, memory, threshold, depth_threshold): 
# memory removed, new exonerate_hits.py does not have that
def run_exonerate_hits(folder,file_, ref_seq_file, threshold, prefix, coverage_depth_multiplier):    
    logging.info("Extracting genes from: " +file_)
    fline=open(file_).readline()
    regex_spades_header =re.search("^>NODE_[0-9]+_length_[0-9]+_cov_[0-9]+",fline)

    #print("EXONERATE MEMORY PER SAMPLE IS: {}GB".format(memory))
    #if spades assembly, run exonerate_hits from HybPiper

    if regex_spades_header != None:
        os.chdir(folder)
        #exonerate_command = "python3 {}exonerate_hits.py -m {} -t {} {} --prefix {} {} ".format(dependencies_dir, memory, threshold, ref_seq_file, os.path.splitext(file_)[0], file_) ###~~~
        exonerate_command = "python3 {}exonerate_hits_orig.py  {} {} --thresh {} --prefix {} --keep_intermediate_files --verbose_logging --no_stitched_contig --depth_multiplier {} ".format(dependencies_dir,  ref_seq_file,  file_, threshold, prefix, coverage_depth_multiplier)
        logging.info(exonerate_command)
        os.system(exonerate_command)
    # else use the script version not using coverage information
    else:
        os.chdir(folder)
        #exonerate_command = "python3 {}exonerate_alt.py -m {} -t {} {} --prefix {} {} ".format(dependencies_dir, memory, threshold, ref_seq_file, os.path.splitext(file_)[0], file_) ###~~~
        exonerate_command = "python3 {}exonerate_hits_orig_alt.py  {} {} --thresh {} --prefix {} --keep_intermediate_files --verbose_logging --no_stitched_contig --depth_multiplier 0".format(dependencies_dir,  ref_seq_file,  file_, threshold, prefix)
        logging.info(exonerate_command)
        os.system(exonerate_command)
    
def run_exonerate(data_dir):
    path_to_assemblies = data_dir
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
        #print(path_to_assemblies)
        for name in files:
            if name.endswith(".fna"): #or name.endswith(".fasta"):
                #print(name)
                pezizo_list.append(root + name)
                assemblies_count += 1
    #print(assemblies_count)                        
    #print("Samples are: ", pezizo_list)
    ref_list = []
    for k in os.listdir(path_to_assemblies):
        if k.endswith("_best_blast_scoring_reference_Hybpiper_format_aa.fas"):
            ref_list.append(path_to_assemblies + k)
    #print(ref_list)
    
    #memory = 0
    #if assemblies_count < int(args.cpu):
    #    if assemblies_count == 0:
    #        assemblies_count = 1
    #    memory = int(args.exonerate_mem / assemblies_count)
    #else:
    #    memory = int(args.exonerate_mem / int(args.cpu))
    #if memory == 0:
    #    memory = 1

    #print("REF LISTS ARE: ",ref_list)
    for z in pezizo_list:
        for i in ref_list:
            regex_fna = re.search(f"{data_dir}(.+?)\.fna", z)
            regex_ref = re.search(f"{data_dir}(.+?)_best_blast_scoring_reference_Hybpiper_format_aa.fas", i)
            if regex_fna.group(1) == regex_ref.group(1):                
                os.mkdir(data_dir + regex_fna.group(1))
                shutil.move(z,data_dir + regex_fna.group(1) + "/" + regex_fna.group(1) + ".fna")
                with open(i, 'r') as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        with open(data_dir + "/" + regex_fna.group(1) + "/" + record.id + "_best_blast_scoring_reference_Hybpiper_format_aa.fas", "w") as output_hand:
                            SeqIO.write(record, output_hand,"fasta")
                for d in os.listdir(data_dir + regex_fna.group(1)):
                    regexx = re.search("(.+?)-(.+?)_best_blast_scoring_reference_Hybpiper_format_aa.fas", d)
                    if regexx:
                        os.mkdir(data_dir + "/" + regex_fna.group(1) +"/"+ regexx.group(2))
                        shutil.move(data_dir + regex_fna.group(1) +"/"+ d, data_dir + "/" + regex_fna.group(1) +"/"+ regexx.group(2))
                        for file in os.listdir(data_dir + regex_fna.group(1)):
                            name, ext = os.path.splitext(file)
                            if ext == ".fna":
                                # instead of coying the assembly in each folder and rename, better to symlink to save space
                                #os.copy(os.path.join(data_dir + regex_fna.group(1),file), os.path.join(data_dir, regex_fna.group(1), regexx.group(2)))
                                #os.rename(os.path.join(data_dir, regex_fna.group(1), regexx.group(2), file), os.path.join(data_dir, regex_fna.group(1), regexx.group(2), regexx.group(2) + ".fna"  ))
                                os.symlink(os.path.join(data_dir, regex_fna.group(1), file), os.path.join(data_dir, regex_fna.group(1), regexx.group(2), regexx.group(2) + ".fna"  ))
    for subdir in os.listdir(data_dir):
        if os.path.isdir(os.path.join(data_dir,subdir)) and not os.path.join(data_dir,subdir).endswith("_spades"):
            current_dir = os.path.join(data_dir,subdir)
            list_of_lists = []
            for dirs in os.listdir(current_dir):
                if os.path.isdir(os.path.join(current_dir,dirs)):
                    files = os.listdir(os.path.join(current_dir,dirs))
                    for i in files:
                        if i.endswith(".fna"):
                            for j in files:
                                if j.endswith(".fas"):
                                    empty_list = []
                                    empty_list.append(os.path.join(current_dir,dirs))
                                    empty_list.append(os.path.join(current_dir,dirs,i))
                                    empty_list.append(os.path.join(current_dir,dirs,j))
                                    empty_list.append(int(args.threshold))
                                    #take the name of the sample as prefix, which is the name of the folder
                                    empty_list.append(current_dir.split("/") [-1])
                                    empty_list.append(int(args.depth_multiplier)) 
                                    list_of_lists.append(empty_list)                        
            logging.info(list_of_lists)
            
            logging.info("Running exonerate using exonerate_hits.py script from Hybpiper..")
            args.cpu = int(args.cpu)
            pool = multiprocessing.Pool(processes=args.cpu)
            pool.starmap(run_exonerate_hits, list_of_lists)
            os.chdir(data_dir)

    # as exonerate hits is launched as standalone and was also modified (prefix was not working as expected)
    # my FNA and FAA final files have the name of the file empty as i only have the sample name as prefix and exonerate_hits.py uses [-2] element which does not exist
    # so I need to rename them after they are created, with the name of the gene
    for root, dirs, files in os.walk('.'):
        for name in files:
            if name.endswith(".FNA") and not name.startswith("exonerate_hits_trimmed"):
                gene_name = root.split("/") [-4]
                os.rename(os.path.join(root,name), os.path.join(root, gene_name +".FNA") )
            elif name.endswith(".FAA") and not name.startswith("exonerate_hits_trimmed"):
                gene_name = root.split("/") [-4]
                os.rename(os.path.join(root,name), os.path.join(root, gene_name + ".FAA") )
    # rename also the file in "intron" folder
    for root, dirs, files in os.walk('.'):
        for name in files:
            if name.endswith("_supercontig.fasta"):
                gene_name = root.split("/") [-4]
                os.rename(os.path.join(root,name), os.path.join(root, gene_name +"_supercontig.fasta") )    
    # rename also paralog output
    for root, dirs, files in os.walk('.'):
        for name in files:
            if name.endswith("_paralogs.fasta"): 
                gene_name = root.split("/") [-3]
                os.rename(os.path.join(root,name), os.path.join(root, gene_name + "_paralogs.fasta") )
    # rename also intronerate output
    for root, dirs, files in os.walk('.'):
        for name in files:
            if name.endswith("_intronerate_fasta_and_gff.txt"): 
                gene_name = root.split("/") [-3]
                os.rename(os.path.join(root,name), os.path.join(root, gene_name + "_intronerate_fasta_and_gff.txt") )
    for root, dirs, files in os.walk('.'):
        for name in files:
            if name.endswith("_intronerate_supercontig_individual_contig_hits.fasta"):
                gene_name = root.split("/") [-3]
                os.rename(os.path.join(root,name), os.path.join(root, gene_name + "_intronerate_supercontig_individual_contig_hits.fasta") )
    for root, dirs, files in os.walk('.'):
        for name in files:
            if name.endswith("_supercontig_without_Ns.fasta"):
                gene_name = root.split("/") [-3]
                os.rename(os.path.join(root,name), os.path.join(root, gene_name + "_supercontig_without_Ns.fasta") )

    # generate the exonerate_genelist.txt file and the genes_with_seqs.txt file from hybpiper
    for subdir in os.listdir(data_dir):
        if os.path.isdir(os.path.join(data_dir,subdir)) and not os.path.join(data_dir,subdir).endswith("_spades"):
            current_dir = os.path.join(data_dir,subdir)
            list_of_exonerate_genes = []
            for dir in os.listdir(current_dir):
                # is a dir and FNA dir has .FNA file 
                if os.path.isdir(os.path.join(current_dir,dir)):
                    for root, dirs, files in os.walk(os.path.join(current_dir,dir)):
                        for name in files:
                            if name.endswith(".FNA"):         
                                list_of_exonerate_genes.append(dir)        
                list_of_exonerate_genes = list(set(list_of_exonerate_genes))

            with open(os.path.join(current_dir, "exonerate_genelist.txt"), 'w') as exonerategenes_handle:
                for i in list_of_exonerate_genes:
                    exonerategenes_handle.write(i + "\n")
            # now generate genes_with_seqs.txt
            with open(os.path.join(current_dir, "genes_with_seqs.txt"), 'w') as geneswithseqs_handle:
                for dir in os.listdir(current_dir):
                    for i in list_of_exonerate_genes:
                        for root, dirs, files in os.walk(os.path.join(current_dir,dir)):
                            for name in files:
                                if name.endswith(".FNA") and not name.startswith("exonerate_hits_trimmed"):
                                    if i in name:
                                        with open(root +"/"+ name, 'r') as fastafile:
                                            lines = fastafile.readlines()
                                            gene_length = 0
                                            for j in lines:
                                                if not j.startswith(">"):
                                                    line_length = (len(j))
                                                    gene_length = gene_length + line_length
                                                    
                                            #print(f"For gene {i} the gene length is {str(gene_length)}")
                                            geneswithseqs_handle.write(i + "\t" + str(gene_length) + "\n")

    ## code which does what Hybpiper2 assemble.py does to create reports about stitched contigs and paralog warnings
    # for now I'm skipping this internal stop codon check performed by Hypiper assemble (no warning will be in output), it should not matter if then alignment and trimming for the phylogeny is performed
    # Also Unfate genes were checked for internal stop codons
    # Stitched contigs reports concatenation
    for subdir in os.listdir(data_dir):
        if os.path.isdir(os.path.join(data_dir,subdir)) and not os.path.join(data_dir,subdir).endswith("_spades"):
            current_dir = os.path.join(data_dir,subdir)
            with open(current_dir +"/"+ current_dir.split("/") [-1] + "_genes_with_stitched_contig.csv", 'w') as stitched_cat:
                for root, dirs, files in os.walk(os.path.join(current_dir)):
                    for name in files:
                        if name.startswith("genes_with_stitched_contig"):
                            with open(root+"/"+name, 'r') as filetocat:
                                lines = filetocat.readline()
                                genename_to_add = root.split("/") [-2]
                                lines = lines.replace(",,","," + genename_to_add + ",")
                                stitched_cat.write(lines)
            stitched_cat.close()
    # Putative chimeras reports (IT WILL BE EMPTY if using "assembly first" method, chimeras are checked using the discordant mapping reads (no mapping in assembly first strategy)
    for subdir in os.listdir(data_dir):
        if os.path.isdir(os.path.join(data_dir,subdir)) and not os.path.join(data_dir,subdir).endswith("_spades"):
            current_dir = os.path.join(data_dir,subdir)
            with open(current_dir +"/"+ current_dir.split("/") [-1] + "_genes_derived_from_putative_chimeric_stitched_contig.csv", 'w') as chimera_cat:
                for root, dirs, files in os.walk(os.path.join(current_dir)):
                    for name in files:
                        if name.startswith("putative_chimeric_stitched_contig"):
                            with open(root+"/"+name, 'r') as filetocat1:
                                lines = filetocat1.readline()
                                ###to be checked
                                genename_to_add = root.split("/") [-2]
                                lines = lines.replace(",,","," + genename_to_add + ",")
                                chimera_cat.write(lines)
            chimera_cat.close()
    # Collate report paralog by contig depth (check what happens is depth is not used... assembly not by spades in assembly as initial data)
    for subdir in os.listdir(data_dir):
        if os.path.isdir(os.path.join(data_dir,subdir)) and not os.path.join(data_dir,subdir).endswith("_spades"):
            current_dir = os.path.join(data_dir,subdir)
            with open(current_dir +"/"+ current_dir.split("/") [-1] + "_genes_with_paralog_warnings_by_contig_depth.csv", 'w') as coveragedepth_cat:
                for root, dirs, files in os.walk(os.path.join(current_dir)):
                    for name in files:
                        if name.startswith("paralog_warning_by_contig_depth"):
                            with open(root+"/"+name, 'r') as filetocat2:
                                lines = filetocat2.readline()
                                ###to be checked
                                genename_to_add = root.split("/") [-2]
                                lines = lines.replace(", gene ,",", gene " + genename_to_add + ",")
                                coveragedepth_cat.write(lines)
            coveragedepth_cat.close()                    
    # Collate report for long paralogs
    for subdir in os.listdir(data_dir):
        if os.path.isdir(os.path.join(data_dir,subdir)) and not os.path.join(data_dir,subdir).endswith("_spades"):
            current_dir = os.path.join(data_dir,subdir)
            with open(current_dir +"/"+ current_dir.split("/") [-1] + "_genes_with_long_paralog_warnings.txt", 'w') as longparalog_cat:
                gene_list = []
                for root, dirs, files in os.walk(os.path.join(current_dir)):
                    for name in files:
                        if name.startswith("paralog_warning_long"):
                            with open(root+"/"+name, 'r') as filetocat:
                                lines = filetocat.readlines()
                                for l in lines:
                                    # there can be more than one
                                    gene_list.append(l.split("\t") [0])
                # only one time per gen must be in the output file
                gene_list = list(set(gene_list))
                for j in gene_list:    
                    longparalog_cat.write(j + "\n")                

################################################################################################################################

    
def checkTestContinue(user_input):
    if user_input == "c" or user_input == "C" or user_input == "Continue" or user_input == "continue":
        return True

def check_arg():
    parser = argparse.ArgumentParser(description='UnFATE version 1.0: the wrapper script that brings YOU from target enrichment sequencing data straight to phylogenetic tree inference! See the readme file for data structure and additional info.')
    mandatory_args = parser.add_argument_group("Mandatory or Suggested", "Arguments which are either required for the proper function of main_wrap.py, or should be used.")
    mandatory_args.add_argument('-o', '--out', required=True,
                    help='The directory where output will be placed upon completion. Required.'
                    )
    mandatory_args.add_argument('-c', '--cpu', default= '4', type=int,
                    help='CPU number used by Hybpiper or parallel run of Exonerate, MACSE, etc. Defaults to 4. Recommended.'
                    )

    data_args = parser.add_argument_group("Input-data related", "Arguments relating to input data, at least one is required.")
    data_args.add_argument('-t', '--target_enrichment_data', default= '',
                   help='Path to target enriched data directory. Also metagenomic reads can be run (be careful) with this flag. Files must end with "_R1.fastq[.gz]", "_R2.fastq[.gz]" or "_SE.fastq[.gz]" '
                   )
    data_args.add_argument('-w', '--whole_genome_data', default= '',
                   help='Path to whole genome sequence data directory. Reads files must be end with "_R1.fastq[.gz]", "_R2.fastq[.gz]"',
                   )
    data_args.add_argument('-a', '--assemblies', default= '',
                   help='Path to assemblies directory. Preovide in this folder already assembled geneome (With SPAdes or any other assembler). SPAdes assemblies will exploit gonting coverage for gene copy selection. Files must end with ".fna[.gz]"',
                   )
    data_args.add_argument('-n', '--ncbi_assemblies', nargs = '+',
                   help='Adds samples from the NCBI assembly database, takes a space-delimited list of taxonomic ranks (e.g. Morchella Tuber Fuffaria). Include AUTO to have UnFATE choose samples from the database which are the most similar to yours, to be included in the analysis.'
                   )

    optional_args = parser.add_argument_group("Optional", "Arguments which might be useful depending on your data.")
    optional_args.add_argument('-f', '--first_use', action= 'store_true',
                   help='USE THIS ARGUMENT THE FIRST TIME YOU RUN UnFATE, then do not move its folder. Modifies some paths in MACSE folder.',
                   )
    optional_args.add_argument('-b', '--target_markers', default= 'UnFATE_markers_195.fas',
                   help='Path to a protein reference file, must be in the format required by HybPiper. Defaults to our reference file. Reference file can be improved if sequence from species similar to the samples are added. This will help the most, if the pipeline will be used in the Hybpiper2 (-l and/or -y) mode'
                   )
    optional_args.add_argument('-l', '--low_memory', action= 'store_true',
                   help='Turns off SPAdes assembling of whole genome data before extracting sequences and uses HybPiper2 instead (less memory intensive). In particular if reference sequecne are distantly related to your sample will perform worse (less reads will map) than the default method (assemblying reads first)'
                   )
    optional_args.add_argument('-y', '--targ_hybpiper', action='store_true',
                   help='Turns off metaSPAdes assembling of target enrichment data (default) before extracting sequences, and uses HybPiper2 instead (less memory intensive). Be careful if input fastq have their /1 /2 inn the header of paired reads, if not Hybpiper seems to be unable to distribute reads to folder before the assembly step.'
                   )
    #optional_args.add_argument('-m', '--exonerate_mem',
    #               help='Limits the memory usage of exonerate when running on an assembly (-a or -w without -l or -t without -y). This does not strictly cap memory usage. More information in the README.',
    #               default=256,type=int
    #               )
    optional_args.add_argument('-e', '--threshold', default=55, type=int,
                   help='Threshold for percent identity between contigs and proteins in exonerate_hits. default=55%%'
                   )
    optional_args.add_argument('-x', '--test', action= 'store_true',
                   help='Allows the user to exit early. Each further step of the pipeline needs then to be started by the user explicitly (run will not complete by itself).'
                   )
    optional_args.add_argument('-r', '--trimal', action= 'store_true',
                   help='Uses TrimAl block filtering method.'
                   )               
    optional_args.add_argument('-s', '--strict_filtering', action= 'store_true',
                   help='STRONGLY RECCOMENDED, especially if samples are not from pure cultures (e.g. herbarium samples), and/or samples have a lot of multiple copy markers. If active, it deletes every gene in the sample, that is in multiple copies, except those one that can be defined as single copy with confidence, having one hit with X times the coverage than the other hits (Contaminant fungi should be excluded, most of the time, this way). Without applying this method, after the coverage multiplier check, if that does not help, the copy most similar to the reference will be selected (default Hybpiper2 xonerate_hits.py behaviour. The wrong fungus can be selected if contamination are abundant and/or captured efficiently. In a future version of UnFATE, if exactly two copies of similar coverage, or two high coverage copies that both pass the coverage check are present, a paralogy resolution algorithm will be implemented to retrieve more markers, which are now excluded with this precautionary filtering).'
                   )
    optional_args.add_argument('-d', '--depth_multiplier', default=10, type=float,
                   help='Coverage depth multiplier: it means the sequence with the highest coverage is selected by exonerate_hits.py if the coverage is N times higher than the coverage of the other sequences to select among. It is reccomeneded to get rid low coverage contaminants. Samples/genes in multiple copies that do not pass this coverage check are deleted if the strict filtering (-s) is applied. Otherwise the copy with the highest similarity to the reference is retained (check carefully intermediate results).  default=10x'
                   )
    optional_args.add_argument('-g', '--gappy_out', default=90, type=int,
                   help='It deletes the samples that have more than X percent of markers missing. It operates on the fasta folder, using the data from the markers_retrieved_percentage.py report (gene_length.csv).  default=90%%'
                   )                                                                 
    return parser.parse_args()

args = check_arg()
def main():
    logging.info("Arguments are: {}".format(args))
    testPrompt = "Finished {}, next step is {}: (C)ontinue or (Q)uit? "
    #print(args)
    global main_script_dir
    main_script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "")
    global dependencies_dir
    dependencies_dir = os.path.join(main_script_dir, "dependencies/")
    # raise file limit
    os.system('ulimit -n 1024000')
    # silence Parallel annoying citation notice
    os.system('parallel --citation')
    set_up_directories()    
    
    if args.first_use == True:    
        # commands to modify OMM_MACSE main script and utilities with the right path to them instead of installing Singularity to use it
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

    path_to_macsed_align = os.path.join(args.out, 'macsed_alignments/')
    path_to_supermatrix= os.path.join(args.out, 'supermatrix/')
    if args.trimal:
        path_to_supermatrix_blocked_dna = path_to_supermatrix +"supermatrix_blocked_dna/"
        path_to_supermatrix_blocked_aa = path_to_supermatrix + "supermatrix_blocked_aa/"
    else:
        path_to_supermatrix_dna = path_to_supermatrix +"supermatrix_dna/"
        path_to_supermatrix_aa = path_to_supermatrix + "supermatrix_aa/"
    if os.path.isdir(os.path.join(args.out, "fastas")) and \
        os.path.isdir(os.path.join(args.out, "macsed_alignments")) and \
        os.path.isdir(os.path.join(args.out, "supermatrix")):
        logging.info("fastas/, macsed_alignments/, supermatrix/ already exist, skipping steps.")
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
    
    if args.whole_genome_data:
        logging.info("\033[1;32;50m ***TRIMMING WHOLE GENOME DATA FASTQ FILES  WITH TRIMMOMATC (Bolger et al. 2014)*** \033[1;37;50m")
        trim_and_get_namelist(dependencies_dir, args.whole_genome_data)
        namelist = 'namelist.txt'
        path_to_namelist = os.path.join(args.whole_genome_data, namelist)

        if args.low_memory:
            logging.info("Gunzipping paired reads trimmed fastq archives")
            gunzip_fastq = 'parallel -j {} gunzip ::: {}*_paired.fastq.gz'.format(args.cpu, args.whole_genome_data) 
            os.system(gunzip_fastq)
            gunzip_fastq = 'parallel -j {} gunzip ::: {}*trimmed.fastq.gz'.format(args.cpu, args.whole_genome_data)
            os.system(gunzip_fastq)
            logging.info("\033[1;32;50m ***EXTRACTING GENES FROM WHOLE GENOME DATA  WITH Hybpiper (Johnson et al. 2016)*** \033[1;37;50m")
            run_hybpiper(main_script_dir, args.whole_genome_data, path_to_namelist)

            logging.info("Gzipping paired reads trimmed fastqs")
            gzip_fastq = "parallel -j {} gzip ::: {}*_paired.fastq".format(args.cpu, args.whole_genome_data)
            os.system(gzip_fastq)
            gzip_fastq = "parallel -j {} gzip ::: {}*trimmed.fastq".format(args.cpu, args.whole_genome_data)
            os.system(gzip_fastq)
        else:
            logging.info("\033[1;32;50m ***ASSEMBLING WHOLE GENOME DATA WITH SPAdes (Prjibelski et al., 2020)*** \033[1;37;50m")
            with open(path_to_namelist) as namelistFile:
                for name in namelistFile:
                    name = name.strip()
                    if len(glob(os.path.join(args.out, "whole_genome_data", name + "*R*trimmed_paired.fastq.gz"))) == 2:
                        sample_R1_path = os.path.join(args.out, "whole_genome_data", name + "_R1.trimmed_paired.fastq.gz")
                        sample_R2_path = os.path.join(args.out, "whole_genome_data", name + "_R2.trimmed_paired.fastq.gz")
                        spades_out_path = os.path.join(args.out, "whole_genome_data", name + "_spades/")
                        spades_command = "spades.py -1 {} -2 {} -o {} -t {} --careful --phred-offset 33".format(sample_R1_path, sample_R2_path, spades_out_path, args.cpu)
                        logging.info("running spades with " + spades_command)
                        os.system(spades_command)
                    else:
                        sample_path = os.path.join(args.out, "whole_genome_data", name + "_SE.trimmed.fastq.gz")
                        spades_out_path = os.path.join(args.out, "whole_genome_data", name + "_spades/")
                        spades_command = "spades.py -s {} -o {} -t {} --careful --phred-offset 33".format(sample_path, spades_out_path, args.cpu)
                        logging.info("running spades with " + spades_command)
                        os.system(spades_command)

                    if os.path.isfile(os.path.join(spades_out_path, "scaffolds.fasta")):
                        #print(os.path.join(spades_out_path))
                        os.rename(os.path.join(spades_out_path, "scaffolds.fasta"), os.path.join(args.whole_genome_data, name + ".fna"))
                    else:
                        if os.path.isfile(os.path.join(spades_out_path, "contigs.fasta")):
                            os.rename(os.path.join(spades_out_path, "contigs.fasta"), os.path.join(args.whole_genome_data, name + ".fna"))

            #we have assemblies now
            logging.info("\033[1;32;50m ***PERFORMING ASSEMBLED WGS DATA ANALYSIS WITH Exonerate (Slater & Birney 2005)*** \033[1;37;50m")
            run_exonerate(args.whole_genome_data)

    #create hybpiper output in target_enrichment/
    if args.target_enrichment_data:
        logging.info("\033[1;32;50m ***TRIMMING TARGET ENRICHMENT FASTQ FILES  WITH TRIMMOMATC (Bolger et al. 2014)*** \033[1;37;50m")
        logging.info('Path to TE data: '+args.target_enrichment_data)
        trim_and_get_namelist(dependencies_dir, args.target_enrichment_data)
        namelist = 'namelist.txt'
        path_to_namelist = os.path.join(args.target_enrichment_data,namelist)

        if args.targ_hybpiper:
            #we only need to unzip if using HybPiper
            logging.info("Gunzipping paired reads trimmed fastq archives")
            gunzip_fastq =' parallel -j {} gunzip ::: {}*_paired.fastq.gz'.format(args.cpu, args.target_enrichment_data)
            os.system(gunzip_fastq)  # We don't care about the unpaired reads if paired end are used
            gunzip_fastq = 'parallel -j {} gunzip ::: {}*trimmed.fastq.gz'.format(args.cpu, args.target_enrichment_data)
            os.system(gunzip_fastq)  # This covers single end reads
            logging.info("\033[1;32;50m           EXTRACTING GENES FROM TARGET ENRICHMENT DATA  WITH Hybpiper (Johnson et al. 2016) \033[1;37;50m")
            run_hybpiper(main_script_dir, args.target_enrichment_data, path_to_namelist)
            logging.info("Gzipping paired reads trimmed fastqs")
            gzip_fastq = "parallel -j {} gzip ::: {}*_paired.fastq".format(args.cpu, args.target_enrichment_data)
            os.system(gzip_fastq)
            gzip_fastq = "parallel -j {} gzip ::: {}*trimmed.fastq".format(args.cpu, args.target_enrichment_data)
            os.system(gzip_fastq)

        else:
            logging.info("\033[1;32;50m          ASSEMBLING TARGET ENRICHMENT DATA WITH METASPADES (Nurk et al., 2017)          \033[1;37;50m")
            with open(path_to_namelist) as namelistFile:
                for name in namelistFile:
                    name = name.strip()
                    if len(glob(os.path.join(args.out, "target_enrichment_data", name + "_R*trimmed_paired.fastq.gz"))) == 2:
                        sample_R1_path = os.path.join(args.out, "target_enrichment_data", name + "_R1.trimmed_paired.fastq.gz")
                        sample_R2_path = os.path.join(args.out, "target_enrichment_data", name + "_R2.trimmed_paired.fastq.gz")
                        spades_out_path = os.path.join(args.out, "target_enrichment_data", name + "_spades/")
                        spades_command = "spades.py -1 {} -2 {} -o {} -t {} --meta --phred-offset 33".format(sample_R1_path, sample_R2_path, spades_out_path, args.cpu)
                        logging.info("running spades with " + spades_command)
                        os.system(spades_command)
                        if os.path.isfile(os.path.join(spades_out_path, "scaffolds.fasta")):
                            #print(os.path.join(spades_out_path))
                            os.rename(os.path.join(spades_out_path, "scaffolds.fasta"), os.path.join(args.target_enrichment_data, name + ".fna"))
                        else:
                            if os.path.isfile(os.path.join(spades_out_path, "contigs.fasta")):
                                os.rename(os.path.join(spades_out_path, "contigs.fasta"), os.path.join(args.target_enrichment_data, name + ".fna"))
                    
                    else:
                        sys.exit("metaspades does not work with single end reads yet. Please re-run with the -l (--low_memory) flag to use HybPiper instead of metaSPAdes for target enrichment data.")
                    


            #we have assemblies now
            logging.info("\033[1;32;50m ***PERFORMING ASSEMBLED TARGET ENRICHMENT DATA ANALYSIS WITH Exonerate (Slater & Birney 2005)*** \033[1;37;50m")
            run_exonerate(args.target_enrichment_data)
            
    if args.assemblies:
        logging.info("\033[1;32;50m ***PERFORMING ASSEMBLIES DATA ANALYSIS WITH Exonerate (Slater & Birney 2005)       \033[1;37;50m")
        run_exonerate(args.assemblies)

    ###~~~ Here the code for the strict filtering of putative paralogs/alleles/contaminations
    ###~~~ For now every multiple copies genes with no copy having at least 10x (or similar) the coverage of the rest will be removed
    ###~~~ Later on I will implement a more advanced version tha if there are two copies with similar coverage will be
    # treated as real paralogs (ParaGONE pipeline or similar will be implemented), deleting eventual low coverage contaminants before      
    if args.strict_filtering:
        logging.info("\033[1;32;50m ***PERFORMING STRICT FILTERING OF PUTATIVE PARALOGS/CONTAMINATIONS/ALLELES (no genes with multiple copies retained, except if they have at least X times the coverage of the rest)*** \033[1;37;50m")
        for root, dirs, files in os.walk(args.out, topdown=True):
            # if the paralogs folder is created by exonerate_hits.py, otherwise there must be a single copy
            if "paralogs" in dirs:
                paralogs_folder = os.path.join(root, "paralogs")
                sequences_folder = os.path.join(root, "sequences")
                sample_name = os.path.split(root)[-1]
                gene_name = os.path.split(os.path.split(root)[0])[-1]
                paralogs_file = os.path.join(paralogs_folder, gene_name + "_paralogs.fasta")
        
                # Check if the _paralogs.fasta file exists in the paralogs folder
                if os.path.isfile(paralogs_file):
                    # if first header is SPAdes style (contains '_cov_') do the coverage depth check, else remove the gene (or , later, pass it to parGONE to solve paralogy)
                    with open(paralogs_file, "r") as fasta_file:
                        # read first line
                        line = fasta_file.readline()
                        # if line contains _cov_
                        if "_cov_" in line:
                            logging.info(f"Analyzing headers in {paralogs_file}")
                            header_info = {}
                            for record in SeqIO.parse(fasta_file, "fasta"):
                                # Add the header to a dictionary
                                coverage = float(record.description.split("_cov_")[1].split(",")[0])
                                header_info[record.description] = coverage # Add the coverage to the dictionary                  
                            # here code that select the highest value entry of the dictionary
                            max_key = None
                            max_value = float('-inf')  # Initialize with negative infinity to ensure any value will be greater
                            for key, value in header_info.items():
                                if value > max_value:
                                    max_value = value
                                    max_key = key
                            # now i see if the highest coverage is 10x (depth threshold) the other coverages
                            countwins = 0 
                            for key, value in header_info.items():
                                if max_value / value >= args.depth_multiplier:
                                    countwins = countwins + 1
                                else:
                                    countwins = countwins
                           #print(countwins)
                            # if yes then the genes was selected by exonerate_hist.py by coverage depth and the gene can be retained in the sequecnes folder
                            # works only if there is just one copy, but will not be used anyway because paralogs folder should not exist in this case
                            if countwins == len(header_info) -1:
                                logging.info("Gene", gene_name, "for sample", sample_name, "was retained in the sequecnes folder")
                            else:
                                # if the coverage condition is not met the gene is removed
                                # later on i will add code to retain genes if there are two (similar) coverage copies
                                # and maybe other lower coverage (below the threshold), these will be passeded to ParaGONE pipeline
                                # and then when paralogs are resolved the right copy will be put in the sequecnes folder in place
                                # of the one chosen by exonerate_hist.py using similarity to the reference  
                                #    
                                # delete the sequences_folder
                                if os.path.isdir(sequences_folder):
                                    os.system(f"rm -r {sequences_folder}")
                                    logging.info("Gene", gene_name, "for sample", sample_name, "was deleted due to multiple copies with similar coverage (e.g. contaminations, paralogs, alleles)")
                        else:
                            # removed due to more than one copy and no coverage info (can happen if assembly mode is used and the assembly have no SPAdes style header)
                            if os.path.isdir(sequences_folder):
                                    os.system(f"rm -r {sequences_folder}")
                                    logging.info("Gene", gene_name, "for sample", sample_name, "was deleted due to multiple copies and no coverage info available")
                else:
                    logging.info(f"No _paralogs.fasta file found in {sample_name} for {gene_name}")


    
    logging.info("\033[1;32;50m ***BUILDING FASTA FILES*** \033[1;37;50m")
    user_samples = []

    if args.assemblies:
        logging.info("Adding sequences from assemblies data")
        #get_alignment(args.assemblies)
        get_fastas_exonerate(args.assemblies, True)
        user_samples.extend(get_names(args.assemblies, True))
    if args.target_enrichment_data:
        logging.info("Adding sequences from target enrichment data")
        if args.targ_hybpiper:    #done with HybPiper
            get_fastas_exonerate(args.target_enrichment_data, False)
            user_samples.extend(get_names(args.target_enrichment_data, False))
        else:            #done with metaspades then exonerate, but now with same folder structure of hybpiper so i put False also here below
            get_fastas_exonerate(args.target_enrichment_data, False)
            user_samples.extend(get_names(args.target_enrichment_data, False))
    if args.whole_genome_data:
        logging.info("Adding sequences from whole genome data")
        if args.low_memory:    #done with hybpiper
            get_fastas_exonerate(args.whole_genome_data, False)
            user_samples.extend(get_names(args.whole_genome_data, False))

        else:            #done with spades then exonerate , but now with same folder structure of hybpiper so i put False also here below
            get_fastas_exonerate(args.whole_genome_data, False)
            user_samples.extend(get_names(args.whole_genome_data, False))

    ncbi_accessions = set()
    if args.ncbi_assemblies:
        logging.info("Finding accessions for samples with specified taxonomy from database ")

        #unzip pre-mined database if still zipped
        path_to_premined_dna = main_script_dir + "pre_mined_dna.tar.xz"
        if not os.path.isdir(main_script_dir + "pre_mined_dna/"):
            unzip_premined_assemblies = "tar -C {} -Jxf {}".format(main_script_dir, path_to_premined_dna)
            os.system(unzip_premined_assemblies)
            
        path_to_premined_aa = main_script_dir + "pre_mined_aa.tar.xz"
        if not os.path.isdir(main_script_dir + "pre_mined_aa/"):
            unzip_premined_assemblies = "tar -C {} -Jxf {}".format(main_script_dir, path_to_premined_aa)
            os.system(unzip_premined_assemblies)

        #path_to_premined = main_script_dir + "combined_pre_mined_assemblies/"
        path_to_taxonomy = main_script_dir + "Accession_plus_taxonomy_Pezizomycotina.txt"
        for item in args.ncbi_assemblies:
            #print("ITEM IS:", item)
            if item == "AUTO":
                continue    
            # using the commas that are present in the csv file containing the taxonomy prevents to include ranks with similar name to the one requested in the command line (e.g. Fuffaria and Fuffarialongis)
            item1 = "," + item + ","
            with open(path_to_taxonomy, 'r') as taxonomy:
                for line in taxonomy:
                    if item1 in line:
                        ncbi_accessions.add(line.split(",")[0])
                        #don't break because there are probably multiple with each
        logging.info("\033[1;32;50m ***ADDING SELECTED TAXONOMIC RANKS GENES FROM PRE-MINED ASSEMBLY DATABASE*** \033[1;37;50m")
        if "AUTO" in args.ncbi_assemblies:
            logging.info("Using mafft to align your sequences with sequences from the database")
            #align new stuff in fastas/ with pre-aligned stuff
            #make supermatrix
            #find scores (discussion on this elsewhere)
            #add samples we want to ncbi_accessions
            auto_dir = os.path.join(args.out, "auto_selection")
            if not os.path.isdir(auto_dir):
                os.mkdir(auto_dir)
            for file in glob(os.path.join(args.out, "fastas", "*nucleotide*.fasta")):
                geneName = file.split("/")[-1].split("_")[1] #fastas/abc_123at4980_protein_merged.fasta
                mafft_command = "mafft --retree 1 --thread {} --add {} {} > {}"
                prealignedFile = os.path.join(main_script_dir, "pre_mined_dna", "combined_{}.FNA".format(geneName))
                newFile = os.path.join(auto_dir, "added_{}.fasta".format(geneName))
                os.system(mafft_command.format(args.cpu, file, prealignedFile, newFile))

            os.chdir(auto_dir)
            logging.info("Concatenating alignments of user samples and db samples")
            os.system("perl {}FASconCAT-G_v1.04.pl -i -s".format(dependencies_dir))
            os.chdir(main_script_dir)
            #print(user_samples)
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
                
        logging.info(f"ACCESSIONS RETRIEVED FROM DATABASE ARE:  {ncbi_accessions}")
        #Add sequences from database to sequences from supplied data
        path_to_premined_aa = os.path.join(main_script_dir, "pre_mined_aa")
        path_to_premined_dna = os.path.join(main_script_dir, "pre_mined_dna")
        for fasta in glob(os.path.join(args.out, "fastas", "*.fasta")):
            #logging.info("NCBI to fasta FASTA: " + fasta)
            baseName = fasta.split("/")[-1].split("_")[1]
            moleculeType = fasta.split("/")[-1].split("_")[2]
            #logging.info("NCBI to fasta MOLECULETYPE: " + moleculeType)
            if moleculeType == "protein":
                with open(os.path.join(path_to_premined_aa, "combined_" + baseName + ".FAA"),'r') as proteinIn:
                    with open(fasta, 'a+') as proteinOut:
                        for line in proteinIn:
                            if line.startswith(">"):
                                for accession in ncbi_accessions:
                                    if accession in line:
                                        proteinOut.write(line)
                                        sequenceLine = proteinIn.readline()
                                        proteinOut.write(sequenceLine)                                        
            if moleculeType == "nucleotide":
                with open(os.path.join(path_to_premined_dna, "combined_" + baseName + ".FNA"),'r') as nucIn:
                    with open(fasta, 'a+') as nucOut:
                        for line in nucIn:
                            if line.startswith(">"):
                                for accession in ncbi_accessions:
                                    if accession in line:
                                        nucOut.write(line)
                                        sequenceLine = nucIn.readline()
                                        nucOut.write(sequenceLine)
                
    # Remove empty fastas, just in case there are no sequences for a specific gene
    os.system("find {} -type f -empty -delete".format(os.path.join(fastas_directory, "*_merged.fasta")))

    # Get rid of the trash strings after the accession number to be able to replace with species name later
    # As OMM_MACSE will use soft masking to align and trim better get rid of all small case letters in the alignments before running MACSE pipeline
    logging.info("Cleaning sequences names to only retain accession numbers...")
    logging.info("Converting all nucleotides to uppercase...")
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
                SeqIO.write(record, output_file,"fasta")
            else:
                sequence = str(seq.seq).upper()
                record = SeqRecord(Seq(sequence), seq.id, "","")
                SeqIO.write(record, output_file,"fasta")
        output_file.close()

    logging.info("\033[1;32;50m ***COMPARING RETRIEVED GENES TO REFERENCE SEQUENCES LENGTH*** \033[1;37;50m")
    markers_retrieved_percentage_script = os.path.join(dependencies_dir, "markers_retrieved_percentage.py")
    run_markers_retrieved_percentage = "python3 {} -b {} -f {} ".format(markers_retrieved_percentage_script, args.target_markers, fastas_directory)
    os.system(run_markers_retrieved_percentage)

    logging.info("\033[1;32;50m ***GETTING RID OF GAPPY (low amount of markers) SAMPLES*** \033[1;37;50m")
    if args.gappy_out:
        data1 = pandas.read_csv(args.out + "/fastas/gene_length.csv")
        gappy_samples = []
        for index, row in data1.iterrows():
            # counts zeros then divide for lenght of the object minus 1, or the title will be counted
            zeros_percentage = ((row == 0).sum() / (len(row)-1))*100
            if zeros_percentage > args.gappy_out:
                gappy_samples.append(row["Unnamed: 0"])
                #print(gappy_samples)
        logging.info(f"GAPPY SAMPLES:  {gappy_samples}")
        for fasta_f in glob(os.path.join(args.out, "fastas","*_merged_headmod.fas")):
            for item_to_del in gappy_samples:
                # Read the sequences from the FASTA file
                records = [record for record in SeqIO.parse(fasta_f, "fasta") if item_to_del not in record.id]
                #print(records)
                #input("press enter to continue")
                # Write the filtered sequences back to the file (overwrites the existing file)
                with open(fasta_f, "w") as output_handle:
                    SeqIO.write(records, output_handle, "fasta")

    logging.info("\033[1;32;50m ***PERFORMING ALIGNMENT WITH OMM_MACSE (Ranwez et al. 2018; Di Franco et al. 2019)*** \033[1;37;50m")
    os.chdir(fastas_directory)
    MACSE_dir = main_script_dir + "MACSE_V2_PIPELINES/OMM_MACSE/"
    MACSE_script = MACSE_dir + "S_OMM_MACSE_V10.02.sh"
    # As OMM_MACSE uses soft masking put all the sequences in upper case before alignment and filtering (just to be safe)
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
                file_path = root +"/"+ f
                regex1 =re.search("Alignment_(\S+)_nucleotide_merged_headmod.fas_out",root)
                #os.rename both renames and moves files
                os.rename(file_path,  path_to_macsed_align + regex1.group(1) + f + ".fas")

    #gblocks_path = os.path.join(dependencies_dir, "Gblocks")
    if args.trimal:
        logging.info("\033[1;32;50m ***PERFORMING ALIGNMENT FILTERING WITH TrimAl (Capella-Gutiérrez et al. 2009)*** \033[1;37;50m")
        # not needed if gblocks is not used
    if args.trimal:
        for f in os.listdir(path_to_macsed_align):
            if f.endswith(".aln.fas"):
                trimal_command = "trimal -in {} -out {}trimal.fas -automated1".format(os.path.join(path_to_macsed_align,f),(os.path.join(path_to_macsed_align,f)).rstrip("fas"))
                os.system(trimal_command)
    # get rid of the sequences with > 75% of gaps (sequence which remained empty o almost empty after filtering, i.e. not well alignable and/or already short before filtering)
    logging.info("Deleting sequences with > 75% gaps from alignments")
    for f in os.listdir(path_to_macsed_align):
        if f.endswith("-gb") or f.endswith(".fas"):
            #logging.info("Processing %s" % f)
            input_file = path_to_macsed_align + f
            with open(path_to_macsed_align + f +"_cleaned.fasta", "a+") as output_file:
                for record in SeqIO.parse(input_file, 'fasta'):
                    drop_cutoff = 0.75
                    name = record.id
                    seq = record.seq
                    seqLen = len(seq)
                    gap_count = seq.count("-")
                    if seqLen == 0 or (gap_count/seqLen) >= drop_cutoff:
                        logging.info(" %s was removed." % name)
                    else:
                        SeqIO.write(record, output_file , 'fasta')

    logging.info("\033[1;32;50m ***PERFORMING ALIGNMENTS CONCATENATION WITH Fasconcat (Kück & Longo, 2014)*** \033[1;37;50m")
    #path_to_supermatrix= path_to_macsed_align.replace('macsed_alignments/', 'supermatrix/')
    make_supermatrix_folder="mkdir {} ".format(path_to_supermatrix)
    os.system(make_supermatrix_folder)
    if args.trimal:
        os.system("mkdir {}".format(path_to_supermatrix_blocked_dna))
        os.system("mkdir {}".format(path_to_supermatrix_blocked_aa))
    else:
        os.system("mkdir {}".format(path_to_supermatrix_dna))
        os.system("mkdir {}".format(path_to_supermatrix_aa))
    if args.trimal:
        os.system("cp -r {}*macsed_final_align_NT.aln.trimal.fas_cleaned.fasta {}".format(path_to_macsed_align, path_to_supermatrix_blocked_dna))
        os.system("cp -r {}*macsed_final_align_AA.aln.trimal.fas_cleaned.fasta {}".format(path_to_macsed_align, path_to_supermatrix_blocked_aa))
    else:
        os.system("cp -r {}*macsed_final_align_NT.aln.fas_cleaned.fasta {}".format(path_to_macsed_align, path_to_supermatrix_dna))
        os.system("cp -r {}*macsed_final_align_AA.aln.fas_cleaned.fasta {}".format(path_to_macsed_align, path_to_supermatrix_aa))

    FASconCAT_command = 'perl {}FASconCAT-G_v1.04.pl -l -s'.format(dependencies_dir)
    for moleculeType in ["NT", "AA"]:
        if args.trimal:
            if moleculeType == "NT":
                os.chdir(path_to_supermatrix_blocked_dna)
            else:
                os.chdir(path_to_supermatrix_blocked_aa)
            # Remove empty fastas, just in case there are no sequences for a specific gene
            os.system("find *.fasta -type f -empty -delete")
            os.system(FASconCAT_command)
            os.rename('FcC_supermatrix.fas','FcC_supermatrix_blocked_{}.fasta'.format(moleculeType))
            os.rename('FcC_supermatrix_partition.txt','FcC_supermatrix_partition_blocked_{}.txt'.format(moleculeType))
            os.rename('FcC_info.xls','FcC_info_blocked_{}.xls'.format(moleculeType))
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

    logging.info("\033[1;32;50m ***RECONSTRUCTING SINGLE MARKER TREES WITH IQTREE2 (Minh et al. 2020)*** \033[1;37;50m")
    iqtree_script = os.path.join(dependencies_dir, "iqtree2")
    # DNA alignments
    if args.trimal:
        iqtree_parallel1 = "find %s -type f  -name '*_final_align_NT.aln.trimal.fas_cleaned.fasta' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
        os.system(iqtree_parallel1)
    else:
        iqtree_parallel = "find %s -type f  -name '*_final_align_NT.aln.fas_cleaned.fasta' | parallel -j %s %s -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
        os.system(iqtree_parallel)
    # Amino acid alignments
    if args.trimal:
        iqtree_parallel3 = "find %s -type f  -name '*_final_align_AA.aln.trimal.fas_cleaned.fasta' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
        os.system(iqtree_parallel3)
    else:
        iqtree_parallel2 = "find %s -type f  -name '*_final_align_AA.aln.fas_cleaned.fasta' | parallel -j %s %s  -s {} -m MFP -B 1000 -T 1" %(path_to_macsed_align, args.cpu, iqtree_script)
        os.system(iqtree_parallel2)
    # move trees and other iqtree files to the dedicated folder
    path_to_single_trees = os.path.join(args.out, 'single_locus_trees/')
    mkdir_single_trees = "mkdir {}".format(path_to_single_trees)
    os.system(mkdir_single_trees)
    for root, dirs, files in os.walk(path_to_macsed_align, topdown=True):
        for f in files:
            if  f.endswith("treefile") or f.endswith("nex") or f.endswith("parttrees") or f.endswith("gz") or f.endswith("mldist") or f.endswith("log") or f.endswith("iqtree") or f.endswith("contree") or f.endswith("bionj") or f.endswith("best_scheme"):
                os.rename(path_to_macsed_align + f, path_to_single_trees + f)
    if args.test:
        test_input = input(testPrompt.format("single trees", "supermatrix tree"))
        if not checkTestContinue(test_input):
            sys.exit("exited after making single marker trees")
    logging.info("\033[1;32;50m ***RECONSTRUCTING SUPERMATRIX TREE WITH IQTREE2  (Minh et al. 2020)*** \033[1;37;50m")
    iqtree_script = os.path.join(dependencies_dir, "iqtree2")
    if args.trimal:
        iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s --threads-max %s" %(iqtree_script, path_to_supermatrix_blocked_dna + 'FcC_supermatrix_blocked_NT.fasta' , path_to_supermatrix_blocked_dna + 'FcC_supermatrix_partition_blocked_NT.txt', "AUTO", args.cpu)
        os.system(iqtree_on_supermatrix)
        iqtree_on_supermatrix =  "%s -s %s -Q %s -m MFP -B 1000 -T %s --threads-max %s" %(iqtree_script, path_to_supermatrix_blocked_aa + 'FcC_supermatrix_blocked_AA.fasta' , path_to_supermatrix_blocked_aa + 'FcC_supermatrix_partition_blocked_AA.txt', "AUTO", args.cpu)
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
    logging.info("\033[1;32;50m ***RECONSTRUCTING SUPERTREE WITH ASTRAL (Zhang et al. 2018)*** \033[1;37;50m")
    path_to_supertree = os.path.join(args.out,'supertree/')
    os.system("mkdir {}".format(path_to_supertree))
    if args.trimal:
        blocked_text = "blocked_"
    else:
        blocked_text = "" #don't add "blocked_" to output files or directories
    dna_dir = path_to_supertree + "supertree_{}dna/".format(blocked_text)
    aa_dir = path_to_supertree + "supertree_{}aa/".format(blocked_text)
    os.system("mkdir {}".format(dna_dir))
    os.system("mkdir {}".format(aa_dir))
    os.system("cp {}*NT*.treefile {}".format(path_to_single_trees, dna_dir))
    os.system("cp {}*AA*.treefile {}".format(path_to_single_trees, aa_dir))
    os.chdir(dna_dir)
    os.system("cat *.treefile > cat_trees_{}dna.tre".format(blocked_text))
    os.system("rm *.treefile")
    os.chdir(aa_dir)
    os.system("cat *.treefile > cat_trees_{}aa.tre".format(blocked_text))
    os.system("rm *.treefile")
    os.chdir(main_script_dir)
    path_to_astral = main_script_dir + "ASTRAL/astral.5.7.7.jar"
    run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, dna_dir + "cat_trees_{}dna.tre".format(blocked_text), dna_dir + "astral_species_tree_{}dna.tree".format(blocked_text))
    os.system(run_astral)
    run_astral = "java -jar %s -i %s -o %s" %(path_to_astral, aa_dir + "cat_trees_{}aa.tre".format(blocked_text), aa_dir + "astral_species_tree_{}aa.tree".format(blocked_text))
    os.system(run_astral)        
    logging.info("\033[1;32;50m ***COPYING TREES TO final_trees FOLDER*** \033[1;37;50m")
    # Copy the tree file from IQTREE and ASTRAL to a new folder
    path_to_finaltrees = os.path.join(args.out, 'final_trees/')
    make_finaltrees_folder ="mkdir {}".format(path_to_finaltrees)
    os.system(make_finaltrees_folder)
    for root, dirs, files in os.walk(path_to_supermatrix, topdown=True):
        for f in files:
            if f.endswith("treefile"):
                copy_tree = 'cp {} {}'.format(root +"/"+f, path_to_finaltrees)
                os.system(copy_tree)
    for root, dirs, files in os.walk(path_to_supertree, topdown=True):
        for f in files:
            if f.startswith("astral_species_tree"):
                copy_tree = 'cp {} {}'.format(root +"/"+ f, path_to_finaltrees)
                os.system(copy_tree)    
            
    logging.info("\033[1;32;50m ***CONVERTING ACCESSIONS IN THE TREES, if any, TO SPECIES NAME*** \033[1;37;50m")
    # Open one of the suprematrices, making a list of accession sample name
    if args.trimal:
        supermatrix_file = path_to_supermatrix_blocked_dna + 'FcC_supermatrix_blocked_NT.fasta'
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
    get_taxonomy_script = os.path.join(dependencies_dir, "get_taxonomy_with_edirect.py")
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
                    logging.warning("\033[1;31;50mThe following line does not have the expected format for species name, weird strain name format!")
                    logging.warning(line)
    # Use the "Speciesname, Accession" csv file to substitute the Accession numbers with species names using the funcion "from_accession_to_species"
    # Final names after the substitution will be: "speciesname_accessionnumber"
    for treefile in os.listdir(path_to_finaltrees):
        if treefile.endswith("treefile") or treefile.endswith("tree"):        
            from_accession_to_species(accession_species_file, path_to_finaltrees + treefile)

    logging.info("\033[1;32;50m ***PIPELINE COMPLETED!*** \033[1;37;50m")

if __name__=='__main__':
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    main()
