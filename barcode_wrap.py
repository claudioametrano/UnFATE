#!/usr/bin/env python

#reasons this exists:
#Using multiple loci most likely results in higher resolution phylogenies than just ITS (for example)
#

#pipeline:
#trim fastqs (trimmomatic)
#HybPiper
#mafft adds captured sequences to pre-aligned database sequences for captured loci
#gblocks filters previous alignment
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
from re import findall, sub, match
from sys import exit

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--reads", help="List the .fastq(.gz) file(s) for a single sample separated by a space", nargs="+") #expect readfile for now, allow assembly later
parser.add_argument("-o", "--out", help="output directory, this will be made automatically")
parser.add_argument("-c", "--cpu", default=1, help="number of threads to use")
parser.add_argument("-b", "--target_markers", help="file with baits for HybPiper")
args = parser.parse_args()

def trim_and_get_namelist(main_script_dir, data_dir):
  trimming_cmd = "python3 {}/trimmer.py -f {} -c {}".format(main_script_dir, data_dir, args.cpu)
  os.system(trimming_cmd)
  #Get namelist.txt from dataset directory
  namelist_cmd = 'python3 {}/getNameList.py -f {}'.format(main_script_dir, data_dir + "/")
  os.system(namelist_cmd)

def run_hybpiper(main_script_dir, data_dir):
  with open(os.path.join(data_dir, "namelist.txt"), 'r') as f:
    for line in f:
      print("Processing sample:" + line)
      sample_path = ""
      if len(args.reads) == 2:
        sample_path = data_dir + '/' + line.rstrip('\n') + '_R*.trimmed_paired.fastq'
      else:
        sample_path = data_dir + '/' + line.rstrip('\n') + '*trimmed.fastq'
      run_Hybpiper = '{}/HybPiper/reads_first.py -b {} -r {}  --prefix {}/{} --cpu {} '.format(main_script_dir, args.target_markers, sample_path, data_dir, line.strip(), args.cpu)
      print("running HybPiper with: " + run_Hybpiper)
      os.system(run_Hybpiper)
      clean_command = "{}/HybPiper/cleanup.py {}".format(main_script_dir, line.strip())
      os.system(clean_command)
      os.chdir(main_script_dir)

      return line.strip()
      #this only allows for one sample for now

def get_taxes(main_script_dir):
  taxes = {}
  with open(os.path.join(main_script_dir, "Accession_plus_taxonomy_Pezizomycotina.txt")) as taxFile:
    for line in taxFile:
      line = line.strip().split(",")
      taxes[line[0]] = " ".join(line[1].split(" ")[:2])
  return taxes

def add_fastas(path_to_data, isAssemblies, main_script_dir):
#  for moleculeType in ["FNA"]:
  if isAssemblies:
    search_location = os.path.join(path_to_data, "*", "sequences", "FNA", "*")
  else:
    search_location = os.path.join(path_to_data, "*", "*", "*", "sequences", "FNA", "*")
  mafft_command = "mafft --retree 1 --thread {} --add {} {} > {}"
  for filename in glob(search_location):
    geneName = filename.split("/")[-1][:-4]
    inFasta = os.path.join(main_script_dir, "msa_test", "truncatedIDs_" + geneName + ".fasta")
    outFasta = os.path.join(args.out, "fastas", "added_" + geneName + ".fasta")
    os.system(mafft_command.format(args.cpu, filename, inFasta, outFasta))

#    if moleculeType == "FNA":
#    with open(filename) as inFile, open(os.path.join(args.out, "fastas", "Alignment_" + geneName + "_nucleotide_merged.fasta"), 'a') as outFile:
#      outFile.write(inFile.read())
#  if moleculeType == "FAA":
#    with open(filename) as inFile, open(os.path.join(args.out, "fastas", "Alignment_" + geneName + "_protein_merged.fasta"), 'a') as outFile:
#      outFile.write(inFile.read())

def set_up_dirs():
  if not os.path.isdir(args.out):
    os.mkdir(args.out)
  directories = ["reads", "fastas", "final_fastas", "trees"]
  for dir in directories:
    if not os.path.isdir(os.path.join(args.out, dir)):
      os.mkdir(os.path.join(args.out, dir))

  #get data
  for file in args.reads:  #args.reads is now a list
    if not os.path.exists(os.path.join(args.out, "reads", file.split("/")[-1])):
      os.symlink(file, os.path.join(args.out, "reads", file.split("/")[-1]))

#we don't want all of the recovered database samples to be from the same species
#so we go through scores and take the top 20 where there are <= 4 of each species (not including query, of course)
def get_scores_low_replication(scores, ids, taxes):
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
        if speciesKept[taxes[id]] < 4: #if fewer than 4 have been kept
          scoresToKeep.append(scoreTuple)
          speciesKept[taxes[id]] += 1
        else: 	#more than 4 have been kept, so we don't care
          pass
      else:	#first of this species
        scoresToKeep.append(scoreTuple)
        speciesKept[taxes[id]] = 1
    else: 	#almost certainly the query, unless the taxonomy stuff got messed up
      scoresToKeep.append(scoreTuple)

    if len(scoresToKeep) >= 20:
      return scoresToKeep

def get_score(seq1, seq2, index, calc, total):
  score = calc._pairwise(seq1, seq2)
  if index % 100 == 0:
    print("Processed {} of {}".format(index, total))
  return((index, score))

def find_similar_samples(query, taxes):
  aln = AlignIO.read(open(os.path.join(args.out, "fastas", "FcC_supermatrix.fas")), "fasta")
  for record in aln:
    record.seq = record.seq.upper()
    #record.seq = record.seq.replace("W", "N")
    #record.seq = record.seq.replace("S", "N")
    #record.seq = record.seq.replace("M", "N")
    #record.seq = record.seq.replace("K", "N")
    #record.seq = record.seq.replace("R", "N")
    #record.seq = record.seq.replace("Y", "N")
    #record.seq = record.seq.replace("B", "N")
    #record.seq = record.seq.replace("D", "N")
    #record.seq = record.seq.replace("H", "N")
    #record.seq = record.seq.replace("V", "N")

  #distances will be calculated with a basic 2-parameter model
  calc = DistanceCalculator("trans", skip_letters=["N", "-", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V"])

  queryIndex = [record.id for record in aln].index(query)
  sequences = [record.seq for record in aln]
  ids = [record.id for record in aln]
  querySeq = sequences[queryIndex]

  scores = [] #list of tuples to allow sorting

  with multiprocessing.Pool(int(args.cpu)) as p:
    scores=p.starmap(get_score, [(querySeq, seq, i, calc, len(sequences)) for i,seq in enumerate(sequences)])

  scores.sort(key = lambda x: x[1])
  #print(scores)
  #topScores = scores[0:21]

  topScores = get_scores_low_replication(scores, ids, taxes)

  print(topScores)
  print([(taxes["_".join(ids[i[0]].split("_")[:2])], i[1]) for i in topScores if "_".join(ids[i[0]].split("_")[:2]) in taxes])
  topSamples = [ids[i[0]] for i in topScores]
  print(topSamples)
  return topSamples

def extract_similar_samples(samples):
  print(samples)
  for file in glob(os.path.join(args.out, "fastas", "added_*.fasta")):
    baseName = file.split("/")[-1]
    dirs = "/".join(file.split("/")[:-1])
    with open(os.path.join(dirs, "few_" + baseName[6:]), 'w') as outFile:
      aln = SeqIO.parse(open(file), "fasta")
      recordsToKeep = []
      for record in aln:
        #this is because the ids in msa/ were truncated to allow for gblocks
        #print(record.id)
        if len(record.id) > 30:
          if record.id[:30] in samples:
            recordsToKeep.append(record)
        else:
          if record.id in samples:
            recordsToKeep.append(record)

      SeqIO.write(recordsToKeep, outFile, "fasta")

def align_similar_samples():
  mafft_command = "mafft --thread {} --auto {} > {}"
  for file in glob(os.path.join(args.out, "fastas", "few*.fasta*")):
    baseName = file.split("/")[-1]
    os.system(mafft_command.format(args.cpu, file, os.path.join(args.out, "fastas", "aligned_" + baseName)))

def run_gblocks(isFew, main_script_dir):
  files = []
  if isFew:
    files = glob(os.path.join(args.out, "fastas", "aligned_few*.fasta"))
  else:
    files = glob(os.path.join(args.out, "fastas", "added_*.fasta"))
  for file in files:
    fraction1=0.56
    fraction2=0.56
    b3=10
    b4=5
    print( 'File being processed: {}'.format(file))

    count = 0
    with open(file) as inFile:
      count = inFile.read().count(">")
    print("The alignment has: ",count," sequences")
    b1 = str(round(count * fraction1))
    b2 = str(round(count * fraction2))
    print("Number of char in a column of the alignment to be considered conserved and flanking regions, respectively: ", b1, b2)
    start_Gblocks = "{} {} -t=d -b1={} -b2={} -b3=10 -b4=5 -b5=h -e=-gb".format(os.path.join(main_script_dir, "Gblocks"), file, b1, b2)
    print(start_Gblocks)

    #*.fasta -> *.fasta-gb
    os.system(start_Gblocks)

    #*.fasta-gb -> *.fasta_temp
    removeSpaces = "cat {} | sed 's/ //g' > {}"
    os.system(removeSpaces.format(file + "-gb", file + "_temp"))

    #*.fasta_temp -> *.fasta
    os.remove(file)
    os.rename(file + "_temp", file)

  if isFew:
    for file in glob(os.path.join(args.out, "fastas", "aligned_few*.fasta")):
      baseName = file.split("/")[-1]
      os.rename(file, os.path.join(args.out, "final_fastas", "gb_" + baseName))
  else:
    for file in glob(os.path.join(args.out, "fastas", "*-gb.htm")):
      #file is *.fasta-gb.htm
      #we want to remove *.fasta-gb.htm
      os.remove(file)

def make_trees(main_script_dir):
  iqtree_command = "{}/iqtree2 -m TEST -s {} -T AUTO --threads-max {} --prefix {}/final_fastas"
  os.system(iqtree_command.format(main_script_dir, os.path.join(args.out, "final_fastas"), args.cpu, os.path.join(args.out, "trees")))

def rename_terminals(main_script_dir):
  with open(os.path.join(main_script_dir, "combined_pre_mined_assemblies", "Accession_plus_taxonomy_Pezizomycotina.txt")) as tax_file, \
      open(os.path.join(args.out, "trees", "final_fastas.treefile")) as tree_in, \
      open(os.path.join(args.out, "trees", "final_fastas_named.treefile"), "w") as tree_out:
    tree = tree_in.read()
    accessions = findall(r"GCA_[\d]+\.\d", tree)
    accessionToSp = {}
    for line in tax_file:
      line = line.strip().split(",")
      if line[0] in accessions:
        accessionToSp[line[0]] = line[1]

    accPlusRestRegex = r"{}[^:]+" #{} and .format works with regex
    for accession, sp in accessionToSp.items():
      tree = sub(accPlusRestRegex.format(accession), sp, tree)
    tree_out.write(tree)

def main():
  main_script_dir = os.path.realpath(__file__)
  main_script_dir = "/".join(main_script_dir.split("/")[:-1])
  args.out = os.path.realpath(args.out)
  args.target_markers = os.path.realpath(args.target_markers)
  args.reads = [os.path.realpath(read) for read in args.reads]

  print(args.reads)
  if len(args.reads) > 2:
    exit("\nERROR: Too many read files given as input, must be 1 (unpaired) or 2 (paired).")


  #args.read_dir = os.path.realpath(args.read_dir)  #args.reads is now a list

  #this isn't strictly necessary, but FASconCAT-G needs to be run from the
  #out/fastas directory, so we will change to main_script_dir, then to out/fastas,
  #then back to main_script_dir
  os.chdir(main_script_dir)

  set_up_dirs()
  reads_dir = os.path.join(args.out, "reads")

  print("\n***  Trimming  ***\n")
  trim_and_get_namelist(main_script_dir, reads_dir)
  if len(args.reads) == 2:
    gunzip_fastq = 'parallel -j {} gunzip ::: {}*_paired.fastq.gz'.format(args.cpu, reads_dir + "/")
    os.system(gunzip_fastq)
  else:
    gunzip_fastq = 'gunzip {}*trimmed.fastq.gz'.format(reads_dir + "/")
    os.system(gunzip_fastq)

  print("\n***  Running HybPiper ***\n")
  #run_hybpiper reads namelist.txt and returns the top entry
  query_sample = run_hybpiper(main_script_dir, reads_dir)

  #unpackage ncbi assemblies if not already done
  #extracted sequences from query reads will be added to the pre mined assemblies
  if not os.path.isdir(os.path.join(main_script_dir, "combined_pre_mined_assemblies/")):
    unzip_premined_assemblies = "tar -C {} -Jxf {}".format(main_script_dir, os.path.join(main_script_dir, "combined_pre_mined_assemblies.tar.xz"))
    os.system(unzip_premined_assemblies)
  taxes = get_taxes(main_script_dir)


  add_fastas(reads_dir, False, main_script_dir)
  run_gblocks(False, main_script_dir)
  os.chdir(os.path.join(args.out, "fastas"))
  fasconcat_command = "perl {}/FASconCAT-G_v1.04.pl -i -s"
  os.system(fasconcat_command.format(main_script_dir))
  os.chdir(main_script_dir)

  #TESTING test_samples = ['Trypethelium_eluteriae_NAN5_GGAGCTAC-CTCCTTAC_L008', 'GCA_005059845.1_ASM505984v1_genomic', 'GCA_001692895.1_Cenococcum_geophilum_1.58_v2.0_genomic', 'GCA_001455585.1_Dip_scr_CMW30223_v1.0_genomic', 'GCA_001307955.1_ASM130795v1_genomic', 'GCA_001307945.1_ASM130794v1_genomic', 'GCA_015295685.1_ASM1529568v1_genomic', 'GCA_000281105.1_Coni_apol_CBS100218_V1_genomic', 'GCA_000504465.1_CryAan1.0_genomic', 'GCA_001692915.1_Glonium_stellatum_CBS_207.34_v1.0_genomic', 'GCA_002111425.1_ASM211142v1_genomic', 'GCA_010015785.1_Sacpr1_genomic', 'GCA_010093885.1_Aplpr1_genomic', 'GCA_009829795.1_ASM982979v1_genomic', 'GCA_008931885.1_ASM893188v1_genomic', 'GCA_009829455.1_ASM982945v1_genomic', 'GCA_009829845.1_ASM982984v1_genomic', 'GCA_001307885.1_ASM130788v1_genomic', 'GCA_009829855.1_ASM982985v1_genomic', 'GCA_001307935.1_ASM130793v1_genomic', 'GCA_011057605.1_IIL_Mf_1.0_genomic']
  #TESTING query_sample = "DRR234452"


  print("\n*** Finding most similar sequences in database ***\n")
  similar_samples = find_similar_samples(query_sample, taxes)
  print("\n*** Adding genes from HybPiper to genes in database  ***\n")
  extract_similar_samples(similar_samples)
  #TESTING extract_similar_samples(test_samples)

  print("\n*** Aligning sequences ***\n")
  align_similar_samples()
  run_gblocks(True, main_script_dir)
  print("\n*** Making a tree ***\n")
  make_trees(main_script_dir)
  rename_terminals(main_script_dir)

if __name__=="__main__":
  main()

