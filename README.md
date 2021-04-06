# UnFATE

#### Universal Filamentous Ascomycetes Target Enrichement bait set and wrapper script for phylogenetics and genome-based barcoding 

## Please cite: 
The wrapper script relies on many great software developed by other people. If you use this wrapper and bait set please cite:
1. Ametrano et al. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
2. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.
3. Castresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Molecular biology and evolution, 17(4), 540-552.
4. Johnson, M. G., Gardner, E. M., Liu, Y., Medina, R., Goffinet, B., Shaw, A. J., ... & Wickett, N. J. (2016). HybPiper: Extracting coding sequence and introns for phylogenetics from high‐throughput sequencing reads using target enrichment. Applications in plant sciences, 4(7), 1600016. 
5. Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21), 4453-4455. 
6. Kück, P., & Longo, G. C. (2014). FASconCAT-G: extensive functions for multiple sequence alignment preparations concerning phylogenetic studies. Frontiers in zoology, 11(1), 1-8.
7. Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular biology and evolution, 37(5), 1530-1534.
8. Ranwez, V., Douzery, E. J., Cambon, C., Chantret, N., & Delsuc, F. (2018). MACSE v2: toolkit for the alignment of coding sequences accounting for frameshifts and stop codons. Molecular biology and evolution, 35(10), 2582-2584.
9. Slater, G. S. C., & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC bioinformatics, 6(1), 1-11. 
10. Zhang, C., Sayyari, E., & Mirarab, S. (2017, October). ASTRAL-III: increased scalability and impacts of contracting low support branches. In RECOMB international workshop on comparative genomics (pp. 53-75). Springer, cham. 
 
## Installation and use
Many softwares needed to run the wrapper script are included in the repository, except Hybpiper, OMM_MACSE and ASTRAL. As new versions of these software are released, this wrapper could stop working properly.
Im this case, please open a Github issue or contact claudiog.ametrano@gmail.com  

1. Clone the repository with:  
	$ git clone https://github.com/claudioametrano/UnFATE.git  

2. Clone Hybpiper, MACSE and ASTRAL from Github into the wrapper script folder:  
$ cd path/to/UnFATE  
$ python3 main_script.py --first_use  

3. Set file permission recursively for the wrapper script folder (read, write and execute rights, only for the owner):  
$ chmod -R 700 path/to/UnFATE  

4. Read the help section of the script to set up the command line for your analysis:  
$ python3 main_script.py --help  

IMPORTANT NOTE: folder containing data from target enrichment sequencing must be placed in a folder called "target_enrichment", the assemblies must be in a folder called "assemblies", in order to be used by the script.	

5. Run the pipeline using the comman line. For example:  
$ python3 main_wrap.py -b ~/path/to/protein/fasta/protein_markers_aa.fas -t ~/path/to/target_enrichment/ -a ~/path/to/assemblies/ -g -c 8 -n Tuber Morchella  

