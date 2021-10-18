#!/usr/bin/env python3

#workflow:
#trim fastqs (trimmomatic)
#HybPiper
#mafft adds captured sequences to pre-aligned database sequences for captured loci
#gblocks filters previous alignment (Do we actually want to gblock here?)
#FASconCAT-G concatenates filtered alignments
#2-parameter substitution pairwise scores are calculated between query and all other samples
#20 closest samples (with a maximum of 4 per species) are selected from the database
#mafft aligns the captured loci with the loci from only those 20 samples
#gblocks runs on the small alignment
#iqtree2 runs on the small gblocked alignment

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import SeqIO
import argparse
import os
from glob import glob
import multiprocessing
from re import findall, sub, match, search
from sys import exit
from math import ceil

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="List the .fastq(.gz) file(s) for a single sample separated by a space, or a single .fna(.gz). Wildcards are allowed.", nargs="+") #expect readfile for now, allow assembly later
parser.add_argument("-o", "--out", help="output directory, this will be made automatically")
parser.add_argument("-c", "--cpu", default=1, type=int, help="number of threads to use")
parser.add_argument("-b", "--target_markers", help="baitfile used to extract genes. Only our file or a subset of it will work with the pre-mined database. Defaults to our database, needs to be set manually if run from outside the UnFATE directory.", default="UnFATE_markers_195.fas")
parser.add_argument("-y", "--use_hybpiper", help="use hybpiper instead of metaspades and exonerate_hits.py, reduces memory usage but may increase runtime.", action="store_true")
parser.add_argument("-e", "--exonerate_mem", help="memory for exonerate to try to limit itself to. Set this to a number of GB if you see 'Killed' in the output and there are many sequences unexpectedly missing.", type=int, default=256)
parser.add_argument("-t", "--exonerate_threshold", default=55, type=int, help="Threshold for percent identity between contigs and proteins in exonerate_hits. default=55%%")
parser.add_argument("-n", "--num_tips", help="number of samples to compare to the query sample, defaults to 20", type=int, default=20)
parser.add_argument("-m", "--max_per_sp", help="maximum number of samples from one species, defaults to 4", type=int, default=4)
args = parser.parse_args()

def trim_and_get_namelist(exec_dir, data_dir):
  trimming_cmd = "python3 {}/trimmer.py -f {} -c {}".format(exec_dir, data_dir, args.cpu)
  os.system(trimming_cmd)
  #Get namelist.txt from dataset directory
  namelist_cmd = 'python3 {}/getNameList.py -f {}'.format(exec_dir, data_dir + "/")
  os.system(namelist_cmd)

def run_hybpiper(main_script_dir, data_dir):
  with open(os.path.join(data_dir, "namelist.txt"), 'r') as f:
    for line in f:
      print("Processing sample:" + line)
      sample_path = ""
      if len(args.input) == 2:
        sample_path = data_dir + '/' + line.rstrip('\n') + '_R*.trimmed_paired.fastq'
      else:
        sample_path = data_dir + '/' + line.rstrip('\n') + '*trimmed.fastq'
      run_Hybpiper = '{}/HybPiper/reads_first.py -b {} -r {}  --prefix {}/{} --cpu {} '.format(main_script_dir, args.target_markers, sample_path, data_dir, line.strip(), args.cpu)
      print("running HybPiper with: " + run_Hybpiper)
      os.system(run_Hybpiper)
      clean_command = "{}/HybPiper/cleanup.py {}".format(main_script_dir, os.path.join(data_dir, line.strip(), ""))
      os.system(clean_command)
      os.chdir(main_script_dir)

      return line.strip()
      #this only allows for one sample for now

def blast_individual_sample(sample, genes, assemblies_path, cpu):
  #This only runs once in barcode_wrap
  fna = os.path.join(assemblies_path, sample + ".fna")
  # Build a BLAST database for each of the assemblies
  make_db = "makeblastdb -in {} -dbtype nucl -out {}".format(fna, os.path.join(assemblies_path, sample), os.path.join(assemblies_path, ))
  os.system(make_db)
  sample_best_path = os.path.join(assemblies_path, sample + "_best_blast_scoring_reference_Hybpiper_format_aa.fas")
  for gene in genes:
    ref_path = os.path.join(assemblies_path, gene + "_ref.fasta")
    blastout_path = os.path.join(assemblies_path, gene + "_ref.fasta___" + sample + "_blastout.tsv")
    tblastn_command = "tblastn -query {} -db {} -out {} -num_threads {} -outfmt 6".format(ref_path, os.path.join(assemblies_path, sample), blastout_path, cpu)
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
      regex = search(r"(.+?)_ref\.fasta___(.+?)_blastout.tsv", blastout_path.split("/")[-1])
      print("For assembly:{} and gene:{}, the reference sequence will be: {}".format(regex.group(2), regex.group(1), bestScoring))
      for record in SeqIO.parse(refFile, "fasta"):
        if record.id in bestScoring:
          SeqIO.write(record, bestFile, "fasta")
    os.system("rm {}".format(blastout_path))


def select_best_reference_seq(prot_file_path, assemblies_path, cpu):
  """in the assembly mode this uses exonerate_hits.py script from Hybpiper, there is no way to know what sequences from the protein file is the best for each marker,
     as we do not have the BLAST or the BWA results mapping the reads from target enrichment to those sequences, so we just BLAST all the reference for each gene
     on the assemblies, the one with best bitscore is used to run exonerate_hits.py for that gene.
  """
  assemblies_path = os.path.join(assemblies_path, "")
  ref_gene_list = []
  # Get a list of genes without redundancy
  for seq in SeqIO.parse(prot_file_path,"fasta"):
    geneName = seq.id.strip().split("-")[1]
    ref_gene_list.append(geneName)

  ref_gene_list = list(set(ref_gene_list))
  print("Gene list: ")
  print(ref_gene_list)
  #  Generate a separate fasta for each gene in the reference sequnces file
  for gene in ref_gene_list:
    with open(assemblies_path + gene + "_ref.fasta", "a+") as gene_file:
      for seq in SeqIO.parse(prot_file_path,"fasta"):
        geneName = seq.id.strip().split("-")[1]
        if geneName == gene:
          SeqIO.write(seq, gene_file, "fasta")
        else:
          pass

  fnas = glob(os.path.join(assemblies_path, "*.fna"))
  samples = [fna.split("/")[-1][:-4] for fna in fnas] #fna = "/gar/abc.fna", keep "abc"

  blast_individual_sample(samples[0], ref_gene_list, assemblies_path, cpu)

  # remove all the garbage needed to blast
  # find is used in order to pass to the remove command one file at a time, otherwise if there are too many files the rm command throws the error:"-bash: /bin/rm: Argument list too long"
  extensions = ["*_ref.fasta", "*.nin", "*.nsq", "*.nhr", "*.ndb", "*.nto", "*.not", "*.ntf"]
  for extension in extensions:
    remove_lot_of_files = "find %s -type f -name '%s' -exec rm {} \;" %(assemblies_path, extension)
    os.system(remove_lot_of_files)

#modified from https://biopython.org/wiki/Split_large_file
def batch_iterator(iterator, batch_size):
  """Returns lists of length batch_size.

  This can be used on any iterator, for example to batch up
  SeqRecord objects from Bio.SeqIO.parse(...), or to batch
  Alignment objects from Bio.AlignIO.parse(...), or simply
  lines from a file handle.

  This is a generator function, and it returns lists of the
  entries from the supplied iterator.  Each list will have
  batch_size entries, although the final list may be shorter.
  """
  entry = True  # Make sure we loop once
  while entry:
    batch = []
    while len(batch) < batch_size:
      try:
        entry = next(iterator)
      except StopIteration:
        entry = None
      if entry is None:
        # End of file
        break
      batch.append(entry)
    if batch:
      yield batch

def split_markers_file(data_dir, markersFile, cpu):
  # we will split the _best_blast* file into multiple segments to parallelize exonerate,
  # but each part needs to be in its own file so they can have the same sample name
  #make dirs for each part
  for i in range(cpu):
    os.mkdir(os.path.join(data_dir, "part{}".format(i)))

  #find number of seqs in markersFile in case it isn't 195
  with open(markersFile) as inFile:
    contents = inFile.read()
    numGenes = contents.count(">")

  #split the file
  record_iter = SeqIO.parse(open(markersFile), "fasta")
  for i, batch in enumerate(batch_iterator(record_iter, ceil(numGenes / cpu))):
    filename = os.path.join(data_dir, "part{}".format(i), "part{}_ref.fas".format(i))
    with open(filename, "w") as outFile:
      count = SeqIO.write(batch, outFile, "fasta")
    print("Wrote {} records to {}".format(count, filename))

def run_exonerate_part(data_dir, fna_file, dependencies_dir, partNum, sampleName, mem, thresh, useHits):
  print("MAX EXONERATE MEMORY PER PART IS: {}GB. It will probably use much much less.".format(mem))
  #if spades assembly, run exonerate_hits from HybPiper
  if useHits:
    exonerate_command = "python3 {}exonerate_hits.py -m {} -t {} {} --prefix {} {} ".format(dependencies_dir, mem, thresh, os.path.join(data_dir, "part{}".format(partNum), "part{}_ref.fas".format(partNum)), os.path.join(data_dir, "part{}".format(partNum), sampleName), fna_file)
    print(exonerate_command) #os.path.join(data_dir, "part{}".format(partNum), "part{}_ref.fas".format(partNum)), os.path.join(data_dir, "part{}".format(partNum), sampleName), fna_file)
    os.system(exonerate_command)
  # else use the script version not using coverage information
  else:
    exonerate_command = "python3 {}exonerate_alt.py -m {} -t {} {} --prefix {} {} ".format(dependencies_dir, mem, thresh, os.path.join(data_dir, "part{}".format(partNum), "part{}_ref.fas".format(partNum)), os.path.join(data_dir, "part{}".format(partNum), sampleName), fna_file)
    print(exonerate_command)
    os.system(exonerate_command)

def run_exonerate(data_dir, dependencies_dir, markers, cpu, mem, thresh):
  #to parallelize this step, the best_blast_hits are split into groups of ceil(195/cpu)
  #we make part{0..cpu-1} directories, then use part#/regularPrefix as prefix so that sampleID is kept
  #add_fasta then just needs to look one wildcard deeper than before
  for file in glob(os.path.join(data_dir, "*.fna.gz")):
    os.system("gunzip -f {}".format(file))

  select_best_reference_seq(markers, data_dir, cpu)

  pezizo_list = glob(os.path.join(data_dir, "*.fna"))
  ref_list = glob(os.path.join(data_dir, "*_best_blast_*.fas"))
  fna_file = pezizo_list[0]
  fline = open(fna_file).readline()

  split_markers_file(data_dir, ref_list[0], cpu)

  regex_spades_header = search("^>NODE_[0-9]+_length_[0-9]+_cov_[0-9]+",fline)
  if regex_spades_header != None:
    useHits = True
  else:
    useHits = False

  sampleName = os.path.splitext(fna_file.split("/")[-1])[0]
  partMem = ceil(mem / cpu)

  #this is the list of lists of arguments for running exonerate on each part of the input.
  #everything stays the same except for the part number, and the number of parts is equal to the number of cpus
  exonerate_part_args = [[data_dir, fna_file, dependencies_dir, i, sampleName, partMem, thresh, useHits] for i in range(cpu)]

  pool = multiprocessing.Pool(processes=cpu)
  pool.starmap(run_exonerate_part, exonerate_part_args)

def get_taxes(main_script_dir):
  taxes = {} #{accession: listed sp.}
  fullTaxes = {} #{accession: sp. plus strain number}
  with open(os.path.join(main_script_dir, "Accession_plus_taxonomy_Pezizomycotina.txt")) as taxFile:
    for line in taxFile:
      line = line.strip().split(",")
      if line[1].split(" ")[1] == "cf.":
        #This is because of issues with fusarium, where a lot of samples don't have firm ids
        #There might be a better way to show species that works in all cases, but I haven't found it yet.
        taxes[line[0]] = " ".join(line[1].split(" ")[:3])
      else:
        taxes[line[0]] = " ".join(line[1].split(" ")[:2]) #this works for the majority of entries, where standard bionomials are used
      fullTaxes[line[0]] = line[1]
  return taxes, fullTaxes

def add_fastas(path_to_data, isAssemblies, main_script_dir):
  if isAssemblies:
    search_location = os.path.join(path_to_data, "*", "*", "sequences", "FNA", "*")
  else:
    search_location = os.path.join(path_to_data, "*", "*", "*", "sequences", "FNA", "*")
  mafft_command = "mafft --retree 1 --thread {} --add {} {} > {}"
  for filename in glob(search_location):  #this only processes genes that were recovered, not all genes in the marker file or database
    geneName = filename.split("/")[-1][:-4]
    inFasta = os.path.join(main_script_dir, "pre_mined_dna", "combined_" + geneName + ".FNA")
    outFasta = os.path.join(args.out, "fastas", "added_" + geneName + ".fasta")
    os.system(mafft_command.format(args.cpu, filename, inFasta, outFasta))

def set_up_dirs():
  if not os.path.isdir(args.out):
    os.mkdir(args.out)
  directories = ["input", "fastas", "final_fastas", "trees"]
  for dir in directories:
    if not os.path.isdir(os.path.join(args.out, dir)):
      os.mkdir(os.path.join(args.out, dir))

  #get data
  for file in args.input:  #args.input is now a list
    if not os.path.exists(os.path.join(args.out, "input", file.split("/")[-1])):
      os.symlink(file, os.path.join(args.out, "input", file.split("/")[-1]))

#we don't want all of the recovered database samples to be from the same species
#so we go through scores and take the top totalNum where there are <= spNum of each species (not including query, of course)
#Seperating by species works in most groups, but "<Genus> sp." appears as one species in this stage.
def get_scores_low_replication(scores, ids, taxes, totalNum, spNum):
  scoresToKeep = []
  speciesKept = {}
  accRe = r"GCA_[\d]+\.\d"
  for scoreTuple in scores: #scores must already be sorted
    score = scoreTuple[1]
    id = scoreTuple[0] #eg. 123
    id = ids[id]       #eg. GCA_4321.1_more_stuff
    try:
      id = match(accRe, id)[0]	#eg. GCA_4321.1
    except:
      print("id does not start with GCA, so it must be the query.")

    #start logic of selecting samples
    if id in taxes.keys(): #if not query, essentially
      if taxes[id] in speciesKept.keys(): #if at least one of that species has been kept
        if speciesKept[taxes[id]] < spNum: #if fewer than spNum have been kept
          scoresToKeep.append(scoreTuple)
          speciesKept[taxes[id]] += 1
        else: 	#more than spNum have been kept, so we don't care
          pass
      else:	#first of this species
        scoresToKeep.append(scoreTuple)
        speciesKept[taxes[id]] = 1
    else: 	#almost certainly the query, unless the taxonomy stuff got messed up
      scoresToKeep.append(scoreTuple)

    if len(scoresToKeep) >= totalNum:
      return scoresToKeep

def get_score(seq1, seq2, index, calc, total):
  score = calc._pairwise(seq1, seq2)
  if index % 100 == 0:
    print("Processed {} of {}".format(index, total))
  return((index, score))

def find_similar_samples(query, taxes, totalNum, spNum):
  aln = AlignIO.read(open(os.path.join(args.out, "fastas", "FcC_supermatrix.fas")), "fasta")
  for record in aln:
    record.seq = record.seq.upper()

  #distances will be calculated with a basic 2-parameter model
  calc = DistanceCalculator("trans", skip_letters=["N", "-", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V"])

  queryIndex = [record.id for record in aln].index(query)
  sequences = [record.seq for record in aln]
  ids = [record.id for record in aln] #accessions
  querySeq = sequences[queryIndex]

  scores = [] #list of tuples with index and score [(index, score), ...]

  with multiprocessing.Pool(int(args.cpu)) as p:
    scores=p.starmap(get_score, [(querySeq, seq, i, calc, len(sequences)) for i,seq in enumerate(sequences)])

  scores.sort(key = lambda x: x[1])
  #print(scores)
  #topScores = scores[0:21]

  topScores = get_scores_low_replication(scores, ids, taxes, totalNum, spNum)

  #print(topScores) #not terribly useful or human readable
  print([(taxes["_".join(ids[i[0]].split("_")[:2])], i[1]) for i in topScores if "_".join(ids[i[0]].split("_")[:2]) in taxes]) #prints the listed species name and score
  topSamples = [ids[i[0]] for i in topScores]
  print(topSamples)
  return topSamples

def extract_similar_samples(samples):
  print(samples)
  for file in glob(os.path.join(args.out, "fastas", "added_*.fasta")):
    baseName = file.split("/")[-1]
    dirs = "/".join(file.split("/")[:-1])
    with open(os.path.join(dirs, "few_" + baseName[6:]), 'w') as outFile, open(file) as inFile:
      aln = SeqIO.parse(inFile, "fasta")
      recordsToKeep = []
      for record in aln:
        if record.id in samples:
          recordsToKeep.append(record)
        #this is because the ids in msa/ were truncated to allow for gblocks
	#we should change pre-mined stuff to just have accessions
        #print(record.id)
        elif len(record.id) > 30:
          if record.id[:30] in samples:
            recordsToKeep.append(record)

      SeqIO.write(recordsToKeep, outFile, "fasta")

def align_similar_samples():
  mafft_command = "mafft --thread {} --auto {} > {}"
  for file in glob(os.path.join(args.out, "fastas", "few*.fasta*")):
    baseName = file.split("/")[-1]
    os.system(mafft_command.format(args.cpu, file, os.path.join(args.out, "fastas", "aligned_" + baseName)))

def run_gblocks_individual(exec_dir, file):
  fraction1=0.56
  fraction2=0.56
  b3=10
  b4=5
  print('File being processed: {}'.format(file))

  count = 0
  with open(file) as inFile:
    count = inFile.read().count(">")
  print("The alignment has: ", count ," sequences")
  b1 = str(round(count * fraction1))
  b2 = str(round(count * fraction2))
  print("Number of char in a column of the alignment to be considered conserved and flanking regions, respectively: ", b1, b2)
  start_Gblocks = "{} {} -t=d -b1={} -b2={} -b3=10 -b4=5 -b5=h -e=-gb".format(os.path.join(exec_dir, "Gblocks"), file, b1, b2)

  print(start_Gblocks)
  os.system(start_Gblocks)

def run_gblocks(isFew, exec_dir, cpu):
  files = []
  if isFew:
    files = glob(os.path.join(args.out, "fastas", "aligned_few*.fasta"))
  else:
    files = glob(os.path.join(args.out, "fastas", "added_*.fasta"))

  exec_and_files = [[exec_dir, file] for file in files]

  pool = multiprocessing.Pool(processes=cpu)
  pool.starmap(run_gblocks_individual, exec_and_files)

  for file in files:
    #Gblocks output is *fasta-gb, so we can remove file (*fasta) and replace it with the filtered file
    #os.system(start_Gblocks)
    os.remove(file)

    #*.fasta-gb -> *.fasta (without spaces)
    removeSpaces = "cat {} | sed 's/ //g' > {}"
    os.system(removeSpaces.format(file + "-gb", file))

  if isFew:
    for file in glob(os.path.join(args.out, "fastas", "aligned_few*.fasta")):
      baseName = file.split("/")[-1]
      os.rename(file, os.path.join(args.out, "final_fastas", "gb_" + baseName))
  for file in glob(os.path.join(args.out, "fastas", "*.fasta-gb.htm")):
    #we want to remove *.fasta-gb.htm
    os.remove(file)

def make_trees(exec_dir):
  iqtree_command = "{}/iqtree2 -m TEST -s {} -T AUTO --threads-max {} --prefix {}/final_fastas"
  os.system(iqtree_command.format(exec_dir, os.path.join(args.out, "final_fastas"), args.cpu, os.path.join(args.out, "trees")))

def rename_terminals(fullTaxes):
  #this finds accessions in the tree and replaces them with the full species information listed in the taxonomy file
  with open(os.path.join(args.out, "trees", "final_fastas.treefile")) as tree_in, \
      open(os.path.join(args.out, "trees", "final_fastas_named.treefile"), "w") as tree_out:
    tree = tree_in.read()
    accessions = findall(r"GCA_[\d]+\.\d", tree)

    #{accession: "Genus species (strain info if present)", ...}
    accessionToSp = {accession: fullTaxes[accession] for accession in accessions if accession in fullTaxes.keys()}

    accPlusRestRegex = r"{}[^:]*" #{} and .format works with regex
    for accession, sp in accessionToSp.items():
      tree = sub(accPlusRestRegex.format(accession), sp, tree)
    tree_out.write(tree)

def main():
  main_script_dir = os.path.realpath(__file__)
  main_script_dir = "/".join(main_script_dir.split("/")[:-1]) + "/"
  dependencies_dir = os.path.join(main_script_dir, "dependencies", "")

  args.out = os.path.realpath(args.out)
  args.target_markers = os.path.realpath(args.target_markers)
  args.input = [os.path.realpath(read) for read in args.input]

  print(args.input)
  if len(args.input) > 2:
    exit("\nERROR: Too many input files, must be 1 (unpaired or assembly) or 2 (paired).")
  if args.input[0].endswith(".fna") or args.input[0].endswith(".fna.gz"):
    if len(args.input) > 1:
      exit("\nERROR: Too many input files, must be 1 (unpaired or assembly) or 2 (paired).")

  #args.read_dir = os.path.realpath(args.read_dir)  #args.reads is now a list

  #this isn't strictly necessary, but FASconCAT-G needs to be run from the
  #out/fastas directory, so we will change to main_script_dir, then to out/fastas,
  #then back to main_script_dir
  os.chdir(main_script_dir)

  #unpackage ncbi assemblies if not already done
  if not os.path.isdir(os.path.join(main_script_dir, "pre_mined_dna/")):
    unzip_premined_assemblies = "tar -C {} -Jxf {}".format(main_script_dir, os.path.join(main_script_dir, "pre_mined_dna.tar.xz"))
    os.system(unzip_premined_assemblies)
  taxes, fullTaxes = get_taxes(main_script_dir)
  set_up_dirs()
  input_dir = os.path.join(args.out, "input")

  if args.input[0].endswith("fastq") or args.input[0].endswith("fastq.gz"):
    print("\n***  Trimming  ***\n")
    trim_and_get_namelist(dependencies_dir, input_dir)

    if args.use_hybpiper:
      if len(args.input) == 2:
        gunzip_fastq = 'parallel -j {} gunzip ::: {}*_paired.fastq.gz'.format(args.cpu, input_dir + "/")
        os.system(gunzip_fastq)
      else:
        gunzip_fastq = 'gunzip {}*trimmed.fastq.gz'.format(input_dir + "/")
        os.system(gunzip_fastq)

      print("\n***  Running HybPiper ***\n")
      #run_hybpiper reads namelist.txt and returns the top entry
      query_sample = run_hybpiper(main_script_dir, input_dir)
      add_fastas(input_dir, False, main_script_dir)

    else:
      if len(args.input) == 2:
        sample_R1_path = glob(os.path.join(input_dir, "*_R1.trimmed_paired.fastq.gz"))[0]
        sample_R2_path = glob(os.path.join(input_dir, "*_R2.trimmed_paired.fastq.gz"))[0]
        name = sample_R1_path.strip("_R1.trimmed_paired.fastq.gz")
        spades_out_path = os.path.join(input_dir, name + "_spades/")
        spades_command = "spades.py -1 {} -2 {} -o {} -t {} --meta".format(sample_R1_path, sample_R2_path, spades_out_path, args.cpu)
        print("running spades with " + spades_command)
        os.system(spades_command)
        query_sample = name.split("/")[-1]
      else:
        exit("ERROR: metaspades can't work with SE reads yet.")

      os.rename(os.path.join(spades_out_path, "scaffolds.fasta"), os.path.join(input_dir, name + ".fna"))

      print("Extracting sequences from assembly using exonerate")
      run_exonerate(input_dir, dependencies_dir, args.target_markers, int(args.cpu), args.exonerate_mem, args.exonerate_threshold)
      add_fastas(input_dir, True, main_script_dir)

  elif args.input[0].endswith("fna") or args.input[0].endswith("fna.gz"):
    print("Extracting sequences from assembly using exonerate")
    run_exonerate(input_dir, dependencies_dir, args.target_markers, int(args.cpu), args.exonerate_mem, args.exonerate_threshold)
    add_fastas(input_dir, True, main_script_dir)
    if args.input[0].endswith("fna"):
      query_sample = args.input[0].split("/")[-1][:-4] #name of fna file without extension
    if args.input[0].endswith("fna.gz"):
      query_sample = args.input[0].split("/")[-1][:-7]

  else:
    exit("Input does not end in an allowed extension (fastq(.gz) or fna(.gz)).")

  #extracted sequences have been added to pre-aligned database
  run_gblocks(False, dependencies_dir, int(args.cpu)) #filter large alignments (remove in future?)
  os.chdir(os.path.join(args.out, "fastas"))
  fasconcat_command = "perl {}/FASconCAT-G_v1.04.pl -i -s"
  os.system(fasconcat_command.format(dependencies_dir))
  os.chdir(main_script_dir)

  #use supermatrix of filtered alignments to find similar samples
  print("\n*** Finding most similar sequences in database ***\n")
  similar_samples = find_similar_samples(query_sample, taxes, args.num_tips, args.max_per_sp)
  print("\n*** Adding sequences from query and selected database samples ***\n")
  extract_similar_samples(similar_samples)

  #realign, filter, and make trees from alignments that only contain selected samples
  print("\n*** Aligning sequences ***\n")
  align_similar_samples()
  run_gblocks(True, dependencies_dir, int(args.cpu))
  print("\n*** Making a tree ***\n")
  make_trees(dependencies_dir)
  rename_terminals(fullTaxes)

if __name__=="__main__":
  main()

