# UnFATE

### Languages and Packages

### Software Required
#### Trimmomatic - https://github.com/timflutre/trimmomatic

### Software Installation

### Pipeline Script Layout

##### hyb_wrap.py
This is the wrapper script designed to take user input and direct the pipeline to use amino acid targets or nucleotide targets. 

###### Input for script:
--target_enrichment_data if target enrichment sequencing method 

--whole_genome_data if data is whole genome input data 

--assembly_data if data is either a SPAdes assembly or non-SPAdes assembly


##### getNameList.py
getNameList.py writes a list of sample names to a text file namelist.txt. This text file is then used in the amino acid and nucleotide HybPiper scripts. 


### References 
Johnson, M. G., Gardner, E. M., Liu, Y., Medina, R., Goffinet, B., Shaw, A. J., Zerega, N. J. C, and Wickett, N. J. (2016). HybPiper: Extracting Coding Sequence and Introns for Phylogenetics from High-Throughput Sequencing Reads Using Target Enrichment. Applications in Plant Sciences, 4(7), 1600016. doi:10.3732/apps.1600016

