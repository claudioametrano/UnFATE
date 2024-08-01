# UnFATE 
<img src="./UnFATE_logo.png" alt="Drawing" width="250"/>

## Universal Filamentous Ascomycetes Target Enrichment bait set and wrapper script for phylogenetics and multi-locus barcoding

The wrapper script is designed to be easy to use and to provide a fast way from target enrichment data, assemblies, and/or whole genome sequencing to phylogenetic trees and/or multi-locus barcoding.
If you use the pipline please cite this work and the tools used to build this pipeline (see the section "please cite").

## Workflow
<img src="./pipeline.png" alt="Drawing" height="480">                                                                                                                                                                           
  
## Installation and use

### With Docker container
1. Download the container ready to be run from Docker Hub
  * `docker pull claudioametrano/unfate1.0:latest`
### OR 
### Build yourself the Docker container 
1. clone UnFATE repository  
  * `git clone https://github.com/claudioametrano/UnFATE.git` 
2. Build the docker container (running the Dockerfile which is inside the UnFATE folder, you will not use the rest of the files in the folder, as you will run UnFATE within the container):
  *  `sudo docker build -t unfate:latest  Path/to/UnFATE/folder/`
3. Start the container (interactive)
  *  `sudo docker run -it unfate`
4. Start the UnFATE conda environment within the container
  *  `conda activate unfate`  
5. Quick run in your current directory with the tutorial dataset:
  *  `wget  https://raw.githubusercontent.com/claudioametrano/UnFATE/master/TUTORIAL_DATASET.tar.gz`
  *  `tar -xf TUTORIAL_DATASET.tar.gz`
  *  `python3 main_wrap.py -b ./TUTORIAL_DATASET/10_Unfate_markers_aa.fasta -a ./TUTORIAL_DATASET/assemb_tutorial/ -w ./TUTORIAL_DATASET/WGS_tutorial/ -t ./TUTORIAL_DATASET/TE_tutorial/ -n Letharia -o ./output_wgs_te_ass_letharia --cpu 4 --first_use --trimal --strict_filtering --depth_multiplier 10 --gappy_out 90`
`

### With Conda environment
Prerequisites/Dependencies: 
1. A working Linux operating system (Ubuntu 22.04 LTS and 24.04 LTS were tested; other Linux distributions could work), as the main OS or as a virual machine (e.g. https://www.linuxvmimages.com/images/ubuntu-2404/).
2.  Download and install the Miniconda installer for Linux: 
  * `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
  * `bash ./Miniconda3-latest-Linux-x86_64.sh`
  * Answer "yes" to conda init
3. Restart the Terminal, it should show the (base) conda environment at the beginning of your command line
4. Install mamba (all the conda installations will be faster)
  * `conda install conda-forge::mamba`
5. Create a conda environment which uses Python3.7
  * `mamba create -n unfate python=3.11`
6. Start the environment
  * `mamba activate unfate`  
    * `conda config --add channels conda-forge`
    * `conda config --add channels bioconda`  
7. Install dependencies using conda (mamba)
  * `mamba install -c bioconda blast=2.14.1 spades=3.15.5 exonerate=2.4.0 hmmer=3.3.2 trimal=1.4 mafft=7.520` 
  * `mamba install -c conda-forge biopython=1.80`
  * `mamba install pandas=2.1.1 seaborn=0.13.2 click=8.1.7`
  * `mamba install -c conda-forge parallel=20240722`   
8. Install java if not already installed  
  * `mamba install -n unfate --yes -c conda-forge openjdk=22.0.1`
9. Install Hybpiper2
  * `mamba install hybpiper=2.1.8`
10. clone UnFATE repository  
`git clone https://github.com/claudioametrano/UnFATE.git`
11. Quick run in your current directory with the tutorial dataset:
  * `wget  https://raw.githubusercontent.com/claudioametrano/UnFATE/master/TUTORIAL_DATASET.tar.gz`
  * `tar -xf TUTORIAL_DATASET.tar.gz`
  * `python3 ./UnFATE/main_wrap.py -b ./TUTORIAL_DATASET/10_Unfate_markers_aa.fasta -a ./TUTORIAL_DATASET/assemb_tutorial/ -w ./TUTORIAL_DATASET/WGS_tutorial/ -t ./TUTORIAL_DATASET/TE_tutorial/ -n Letharia -o ./output_wgs_te_ass_letharia --cpu 4 --first_use --trimal --strict_filtering --depth_multiplier 10 --gappy_out 90`
`

### To run Phyparts and phypartspiecharts.py with the pie_wrap.py helper scripts (Optional)
* Create a dedicated enviroment  
  * `mamba create -n phyparts python=3.9.12`
  * `mamba activate phyparts`
* Install phypartspiecharts dependencies
  * `mamba install -c cyclus java-jre=8.45.14` 
  * `mamba install -c bioconda blast=2.12.0 diamond=0.9.21`
  * `mamba install -c etetoolkit ete3=3.1.3`
  * `mamba install -c conda-forge matplotlib=3.8.4`
 
## NOTES
* Rememeber to set up file extensions for your data (if needed): Sequencing data must be in files ending with _R<1|2>.fastq(.gz) or _SE.fastq(.gz). 
Assemblies must be in fasta files ending with .fna(.gz).
* **TUTORIAL_DATASET.tar.gz**:  This reduced dataset and reference sequences file only uses 10 UnFATE genes, assemblies file which only contain the target genes, TE and WGS fastq artificially generated from the same reduced assemblies. Its olny purpose is to test the UnFATE pipeline before you start to work on your own data. It should complete the analyses in a short time, even on a laptop. Please check intermediate results, such as the pdf heatmap in the "fastas" folder. Check also the "final_trees" folder, which should contain a very simple phylogeny containing 17 tips (one sample each of the main Pezizomycotina class from the assemblies, and the same samples from simulated WGS or TE data) plus two _Letharia_ tips from the UnFATE database. For a more computationally intensive test run, wich uses real (downlasampled) data, run the content of **TEST_Data_final.tar.gz** (Running time on 10 Xeon E5-2697v3 cores is about 11').
* Consider running the script from a "tmux" detachable session, as the run can be very long, according to how many samples you have (this tools is usually preinstalled in Linux). Analyses with hundreds of samples should be run on high core number machines!  
* Consider logging the script output with `python3 main_wrap.py {params} |& tee <logfile>`. This saves the stdout and stderr from the running main_wrap.py into \<logfile\> as well as printing it to the console.
* The pre-mined UnFATE database is accessed with the `-n` argument. You can select any taxonomic rank included in Accession_plus_taxonomy_Pezizomycotina.txt. If you want to select a binomial species name, remember to put a backslash before the blank (e.g. Fuffaria\ fuffolosa). If you add AUTO to the list of taxa you want, main_wrap.py will use a similar method to barcode_wrap.py to find the closest samples in the database. Up to 10 additional samples (per user sample) can be added to the dataset this way.
* Although the script was written with our bait set in mind, it should work with any amino acid target file in the format required by HybPiper. If using an external protein file, the UnFATE database will not be helpful, as it only contains the 195 UnFATE genes.
* There are two ways to reduce the memory requirements and CPU burden of UnFATE depending on your input data. If you are supplying whole genome data, you can use the `-l` flag to run HybPiper instead of Spades assembly first and Exonerate. This greatly reduces memory requirements.
* If you wish to use HybPiper to capture sequences from your target enrichment reads, use the `-y` flag. This is kept seperate from the low memory flag which causes HybPiper to be used for WGS data, as assembly of reads from target enrichment is usually not as memory intensive as assembling a whole fungal genome.
* The script parallelizes using the `--cpu n` you specify as an argument, HybPiper and Exonerate will process n samples at a time. Up to the same number of cpu is then used to parallelize IQ-TREE runs for single locus trees and for concatenated supermatrices.  
* A run can be resumed if the script is terminated before generating trees, but after generating supermatrices. This will happen automatically if the output directory contains the assemblies/ and/or target_enrichment/, fastas/, macsed_alignments/, and supermatrix/ directories. Please remove any tree directories from the output directory (if present) before resuming to avoid errors.
* [PhypartsPieCharts](https://github.com/mossmatters/phyloscripts/tree/master/phypartspiecharts) is a nice tool for visualize the nodal conflict level  for the species tree which uses Phyparts. Consider running PhypartsPieCharts through our helper script with `python pie_wrap.py -t /path/to/single_locus_trees/ -p /path/to/species/tree`. Plese use a separate conda environement from the UnFATE one (this is due to some dependendncies version incompatibility)
*  If you need to build a quick phylogeny UnFATE also works using samples from the database only (-n), take advantage of this feature to get phylogenomic trees of any rank in Pezizomycotina with zero effort.     

## Output description
The UnFATE output will be placed in many folders within the location specified by -o, several output folders will be created corresponding to the pipeline steps:  
* The "target_enrichment" folder will contain symlinks to your supplied target enrichment data, trimmed read files, and the HybPiper or metaSPAdes and exonerate_hits.py analysis folders, one per sample.
* The "whole_genome_data" folder will contain symlinks to your supplied WGS data, trimmed read files, and the HybPiper or SPAdes and exonerate_hits.py analysis folders, one per sample.
* The "assemblies" folder will contain symlinks to your assemblies and the "Exonerate_hits.py" runs folder, one per sample.
* The "fastas" folder will contain DNA and AA fasta files, one per marker of interest, the MACSE runs folders, and summaries of the amount of data that could be captured from your data.
* The "macsed_alignments" folder will contain DNA and AA alignments, aligned and filtered with OMM_MACSE pipeline and (optionally) filtered with TrimAl.
* The "auto_selection" folder will exist if you ran main_wrap with `-n AUTO` and will contain DNA alignments of sequences from your data and the pre-mined database.
* The "single_locus_trees" folder will contain the IQ-TREE phylogenetic analyses on single markers (from both DNA and AA alignments).
* The "supermatrix" folder will contain both the concatenation of the single marker alignments and the IQ-TREE phylogenetic inference (from both DNA and AA alignments).
* The "supertree" folder will contain both the file with the best tree for each marker and the ASTRAL species tree (from both DNA and AA alignments).
* The "final_trees" folder will contain the trees generate from concatenation (IQ-TREE) and coalescence-based approach (ASTRAL) and their version renamed to species name (where NCBI accession numbers were used; e.g. when samples from the precalculated database or assemblies downloaded from NCBI are used).
* The "PhyParts" folder will be made if `pie_wrap.py` is run. The key output is pies.svg, but the full phyparts output will be present.
 
## `barcode_wrap.py`
* This script handles multilocus barcoding of target enrichment or WGS data as well as assemblies.
* It finds the most similary samples to the input data in our database (skipping samples if there are already a enough representatives from that species), then builds a tree of those samples using the 195 UnFATe genes. The precision of the taxonomy inferred from the output trees will depend on the completeness of the database. Species level identification could be possible in highly sequenced groups (e.g. Aspergillaceae) but is less likely in groups with few sequenced genomes, as the genomic resources will increase the precision will improve. 
* A possible usage of `barcode_wrap.py` is to find a closely related group to one of your samples, then run `main_wrap.py` with all of your samples and all members of that group using the `-n <taxon>` argument.
* `barcode_wrap.py` does not allow running multiple samples in one run at the moment. If you have multiple samples, consider running `main_wrap.py -n AUTO to get the closest species from the database in your phylogeny`.
* The output directories of `barcode_wrap.py` mirror the output directories of `main_wrap.py`. The "input" directory contains the raw and trimmed reads supplied, as well as the HybPiper or Spades and Exonerate output, if fastqs are supplied. If an assembly is supplied, the contents will be the Exonerate results split into multiple parts.
* The "fastas" directory contains the genes extracted from the input added to the pre-mined data. The "final_fastas" directory contains the genes extracted from the input aligned to the genes from the samples selected from the database, ran through Gblocks with relaxed parameters.
* The "trees" directory contains the IQ-TREE2 output from running on the "final_fastas" directory, in addition to a treefile where the accession numbers from NCBI have been replaced with binomials (final_fastas_named.treefile).

## Please cite: 
The wrapper script relies on many great software developed by other people. If you use this wrapper and bait set please cite the applicable papers:

#### UnFATE
Ametrano et al. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#### Trimmomatic
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.
#### Gblocks
Castresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Molecular biology and evolution, 17(4), 540-552.
### TrimAl
Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), 1972-1973.
### Hmmercleaner
Di Franco, A., Poujol, R., Baurain, D., & Philippe, H. (2019). Evaluating the usefulness of alignment filtering methods to reduce the impact of errors on evolutionary inferences. BMC evolutionary biology, 19(1), 1-17.
#### FASconCAT-G
Kück, P., & Longo, G. C. (2014). FASconCAT-G: extensive functions for multiple sequence alignment preparations concerning phylogenetic studies. Frontiers in zoology, 11(1), 1-8.
#### IQ-TREE 2
Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular biology and evolution, 37(5), 1530-1534.
#### MACSE v2
Ranwez, V., Douzery, E. J., Cambon, C., Chantret, N., & Delsuc, F. (2018). MACSE v2: toolkit for the alignment of coding sequences accounting for frameshifts and stop codons. Molecular biology and evolution, 35(10), 2582-2584.
#### Exonerate
Slater, G. S. C., & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC bioinformatics, 6(1), 1-11.
#### ASTRAL
Zhang, Chao, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153.
#### HybPiper
Johnson, M. G., Gardner, E. M., Liu, Y., Medina, R., Goffinet, B., Shaw, A. J., ... & Wickett, N. J. (2016). HybPiper: Extracting coding sequence and introns for phylogenetics from high‐throughput sequencing reads using target enrichment. Applications in plant sciences, 4(7), 1600016.
#### SPAdes
Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., Korobeynikov, A. (2020). Using SPAdes De Novo Assembler. Current Protocols in Bioinformatics 70(1), e102.
#### metaSPAdes
Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome research 27(5), 824-834.
#### Phyparts
Smith, S. A., Moore, M. J., Brown, J. W., Yang, Y. (2015). Analysis of phylogenomic datasets reveals conflict, concordance, and gene duplications with examples from animals and plants. BMC evolutionary biology 15(1), 1-15
#### Parallel
Tange, O. (2011). Gnu parallel-the command-line power tool. The USENIX Magazine, 36(1), 42-47.
