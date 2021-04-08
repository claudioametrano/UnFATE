# UnFATE

#### Universal Filamentous Ascomycetes Target Enrichement bait set and wrapper script for phylogenetics and genome-based barcoding 

## Please cite: 
The wrapper script relies on many great software developed by other people. If you use this wrapper and bait set please cite:
1. Ametrano et al. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
2. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.
3. Castresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Molecular biology and evolution, 17(4), 540-552.
4. Di Franco, A., Poujol, R., Baurain, D., & Philippe, H. (2019). Evaluating the usefulness of alignment filtering methods to reduce the impact of errors on evolutionary inferences. BMC evolutionary biology, 19(1), 1-17.
5. Johnson, M. G., Gardner, E. M., Liu, Y., Medina, R., Goffinet, B., Shaw, A. J., ... & Wickett, N. J. (2016). HybPiper: Extracting coding sequence and introns for phylogenetics from high‐throughput sequencing reads using target enrichment. Applications in plant sciences, 4(7), 1600016. 
6. Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21), 4453-4455. 
7. Kück, P., & Longo, G. C. (2014). FASconCAT-G: extensive functions for multiple sequence alignment preparations concerning phylogenetic studies. Frontiers in zoology, 11(1), 1-8.
8. Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular biology and evolution, 37(5), 1530-1534.
9. Ranwez, V., Douzery, E. J., Cambon, C., Chantret, N., & Delsuc, F. (2018). MACSE v2: toolkit for the alignment of coding sequences accounting for frameshifts and stop codons. Molecular biology and evolution, 35(10), 2582-2584.
10. Slater, G. S. C., & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC bioinformatics, 6(1), 1-11. 
11. Zhang, Chao, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153.
 

## Workflow
The wrapper is designed to be easy to use and to provide a fast way from target enrichment data (assemblies and whole genome sequence) to phylogenetic tree.
In order to avoid the installation of dependancies and external software, that often lead to problems for inexperienced users (...well, not only for them), most of the software needed by UnFATE is already included in this repository, therefore, aknowledge their work citing them!!.
                                                                                                                                                                           
1. Data from:  **Target enrichment** sequencing (mandatory. Also Whole Genome Sequencing data are accepted), **Assemblies** (optional, also from NCBI assembly database).  
**Representative sequences** used to build the bait set (included in the repository: Target_markers_rep_seq_aa.fas)  

2. **Exonerate 2.2.0** and Exonerate_hits.py script from Hybpiper to mine gene using the amino acid fasta file
                              
3. **Trimmomatic 0.39** to trim Illumina paired reads .fastq.gz from TE or WGS 

4. **HybPiper 1.3.1** to extract the markers from target enrichment sequening reads .fastq files  

5. Fasta files are built from retrieved markers, eventually adding the markers from the pre-mined batabase of NCBI assemblies included in UnFATE repository (if you plan to use it, make sure you download the repository from the browser interface, or install the Github large file storage system before you git clone Unfate)

6. **MACSE2.03** pipline **OMM_MACSE10.02**  to perform codon-aware alignment and segment-based filtering **(HMMcleaner 1.8)**

7. **Gblocks 0.91b** to add an optional second block-based filtering step

8. **RAxML-NG 1.0.2** is used for single locus phylogenetic inference

9. **FASconCAT-G_v1.04** is used to concatenate single marker alignments into a supermatrix

10. **IQTREE2** is used for the supermatrix phylogenetic inference

 11. **ASTRAL 5.7.7** is used to buil the species tree from single locus trees
  
## Installation and use
Prerequisites/Dependancies:  
* A working Ubuntu Linux operating system or Windows10 Linux subsystem:
  *   Many guides are available online, such as: https://www.windowscentral.com/install-windows-subsystem-linux-windows-10
* Python 3 or later (already present in Ubuntu)
* GNU Parallel (already present in Ubuntu)
* Anaconda 
  *  Download and install the Anaconda installer for Linux: https://docs.anaconda.com/anaconda/install/linux/
  * $ bash ~/path/to/Anaconda3-2020.02-Linux-x86_64.sh
  * Answer "yes" to conda init
  * Restart Ubuntu, it should show the (base) conda environment at the beginning of your command line
* Install HybPiper dependancies:  
   * BIOPYTHON 1.59 or later : $ conda install biopython  
   * BLAST command line tools: $ conda install -c bioconda blast 
   * SPAdes: $ conda install -c bioconda spades 
   * EXONERATE: $ conda install -c bioconda exonerate

1. Clone the UnFATE repository with (or download the .zip file from Github browser interface). Chose a position for the UnFATE folder you like, do not move the repository after the first run (use the argument --first_use) 
$ git clone https://github.com/claudioametrano/UnFATE.git  

2. Read the help section of the script to set up the command line for your analysis:  
$ python3 main_script.py --help  

3. Set up the folder structure: sequencing data from target enrichment must be placed in a folder called "target_enrichment", the assemblies must be placed in a folder called "assemblies", in order to be used by the script.	

4. Run UnFATE using the comman line. For example:  
$ python3 main_wrap.py -b ~/path/to/protein/fasta/protein_markers_aa.fas -t ~/path/to/target_enrichment/ -a ~/path/to/assemblies/ --gblocks --cpu 8 -n Tuber Morchella --first_use
  * Remember to use absolute paths for your data folders and the reference protein fasta (e.g. ~/this/is/an/absolute/path/to/the/folder/path/)
  * Consider to run the script from a "tmux" detachable session, as the run can be very long, according to how many samples you have (this tools is usually preinstalled in Linux) 
  * Analyses with hundreds of sample should definitely run on a server-grade hardware!

5.  Cross your fingers and wait, good luck!

## Output description
The UnFATE output will be placed in many folders in the same location of your "target_enrichment" folder, will be created several ouput folders corresponding to the pipeline steps:  
* Within the "target_enrichment" folder there will be the HybPiper runs folder, one per sample  
* Within the "assemblies" folder there will be the "Exonerate_hits.py" runs folder, one per sample  
* "alignments" folder will contain DNA and aa fasta file of the sequences from assemblies (optional) and target enrichment separately, one each marker of interest.  
* "alignments_merged" folder will contain DNA and aa fasta file, one per marker of interest, and the MACSE runs folders.  
* "alignments_merged" folder will contain DNA and aa alignments, aligned and filtered with OMM_MACSE pipeline and their optional version also filtered with Gblocks  
* "single_locus_trees" folder will contain the RAxML phylogenetic analyse on single markers
* "supermatrix" will contain both the concatenation of the single marker alignments and the IQTREE2 phylogenetic inference  
* "supertree" will contain both the file with the best tree for each marker and the ASTRAL species tree
