# UnFATE

#### Universal Filamentous Ascomycetes Target Enrichement bait set and wrapper script for phylogenetics and genome-based barcoding 

## Please cite: 
####The wrapper script relies on many software developed by others. If you use this wrapper and bait set please cite:
1. Ametrano et al. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
2. Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21), 4453-4455. 
3. Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular biology and evolution, 37(5), 1530-1534.
if Gblock option is active:
4. Castresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Molecular biology and evolution, 17(4), 540-552.
5. Kück, P., & Longo, G. C. (2014). FASconCAT-G: extensive functions for multiple sequence alignment preparations concerning phylogenetic studies. Frontiers in zoology, 11(1), 1-8.
6. Zhang, C., Sayyari, E., & Mirarab, S. (2017, October). ASTRAL-III: increased scalability and impacts of contracting low support branches. In RECOMB international workshop on comparative genomics (pp. 53-75). Springer, cham. 
7. Ranwez, V., Douzery, E. J., Cambon, C., Chantret, N., & Delsuc, F. (2018). MACSE v2: toolkit for the alignment of coding sequences accounting for frameshifts and stop codons. Molecular biology and evolution, 35(10), 2582-2584.
8. Johnson, M. G., Gardner, E. M., Liu, Y., Medina, R., Goffinet, B., Shaw, A. J., ... & Wickett, N. J. (2016). HybPiper: Extracting coding sequence and introns for phylogenetics from high‐throughput sequencing reads using target enrichment. Applications in plant sciences, 4(7), 1600016. 
9. Slater, G. S. C., & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC bioinformatics, 6(1), 1-11. 
 
## Installation and data preparation
All softwares needed to run the wrapper script are already included in the repository except Hybpiper, MACSE alignment pipelines and ASTRAL.
As new versions of these software are developed this wrapper could stop working properly. Please contact claudiog.ametrano@gmail.com 
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

