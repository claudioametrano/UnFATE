# UnFATE

### Languages and Packages
UnFATE wrapper script is written in Python3 and uses bash commands and third part softwares written in Java, Perl, C++.
For targer enrichment sequencing data it exploits Hybpiper (Jhonson et al. 2016)...
 
### Installation and data preparation
All softwares needed to run the wrapper script are already included in the repository except Hybpiper, MACSE alignment pipelines and ASTRAL.
1. Clone Hybpiper, MACSE and ASTRAl from Github into the wrapper script folder:
	$ cd path/to/UnFATE
	$ python3 main_script.py --first_use

2. Set file permission recursively for the wrapper script folder (read, write and execute rights, only for the owner):
	$ chmod -R 700 path/to/UnFATE

3. Read the help section of the main scrip to set up the command line for your analysis:
	$ python3 main_script.py --help

IMPORTANT NOTE: folder containing data from target enrichment sequencing should be placed in a folder called "target_enrichment", the assemblies in a folder called "assemblies", in order to be recognized by the script.	

4. Command line example:
	$ python3 main_wrap.py -b ~/path/to/protein/fasta/Hybpiper_markers_aa.fas -t ~/path/to/target_enrichment/ -a ~/path/to/assemblies/ -g -c 8 -n Tuber Morchella

### Pipeline Layout

##### main_wrap.py

###### Input for script:
--target_enrichment_data if target enrichment sequencing method 

--whole_genome_data if data is whole genome input data 

--assembly_data if data is either a SPAdes assembly or non-SPAdes assembly


##### getNameList.py
getNameList.py writes a list of sample names to a text file namelist.txt. This text file is then used in the amino acid and nucleotide HybPiper scripts. 


### References 
Johnson, M. G., Gardner, E. M., Liu, Y., Medina, R., Goffinet, B., Shaw, A. J., Zerega, N. J. C, and Wickett, N. J. (2016). HybPiper: Extracting Coding Sequence and Introns for Phylogenetics from High-Throughput Sequencing Reads Using Target Enrichment. Applications in Plant Sciences, 4(7), 1600016. doi:10.3732/apps.1600016

