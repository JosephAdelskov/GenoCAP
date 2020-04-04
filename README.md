# GenoCAP

## Prokaryotic Genome Phylogenetic and Functional Comparative analysis pipeline

## Introduction
This project consists of a set of 3 bash scripts that when used in sequential order, utilize a variety of genomic analysis tools and databases to perform an automated, phylogenomic and functional comparative analysis of prokaryotic genomes.
These scripts automate the process of retrieving genomes from online databases (NCBI), filtering out poor quality genomes, selection of genomes at the phylogenetic level of genus, while also performing a range of comparative genome analyses.

### Key results produced by this pipeline

1. Universal marker gene phylogenetic tree of all genomes that pass the quality filtering requirements. This produces a view of the target genome(s) of interest.
2. Selection of set of genomes related to target genome(s) at similarity level specified.
2. ANIb (Average nucleotide identity (blast)) analysis of selected set of genomes.
3. Pan/Core genome analysis & statistics, including core genome phylogenetic tree, pan genome parsimony tree.
4. Genome functional annotation and link with KEGG Orthology, COG and Gene Ontology functional data.
5. Pipeline analysis report including important stats and data from the analysis.

### Get\_genomes\_filter\_checkm.sh

This first script uses the NCBI taxid files to find and retrieve all genomes (WGS) related to the taxonomic name/id, rank genus or higher specified by the user. The genomes retrieved are filtered from the analysis by assembly quailty and genome quality (Completeness, contamination and heterogeneity). Once filtered Checkm is then used to extract universal marker genes to perform a multi-gene phylogenetic analysis.

1. Retrieve all genomes (WGS) from NCBI related to the _genus_, _family_, or higher order taxonomy (not tested!) specified.
2. Genomes screened by assembly quality (N contigs, N50) and those under the specified threshold are removed from analysis.
3. Perform checkm lineage_wf to screen genomes for completeness, contamination, and heterogeneity. Filter out genomes below the thresholds specified.
4. Use checkm to get a set of universal marker genes from the remaining genomes.
5. Align marker genes and concatenate into a multi-gene protein alignment for phylogenetic analysis

#### Note on Get\_genomes\_phyla\_amphora.sh

This extra script is the predecessor to get\_genomes\_filter\_checkm.sh that used the program [Phyla Amphora](https://github.com/martinwu/Phyla_AMPHORA) to get the universal marker genes instead of checkm. Phyla Amphora could search for and retrieve universal marker genes using phylum specific predefined sets. This could take a large amount of computing time and even fail to find an acceptable number of universal marker genes (at least 5), if a large number of and/or more recent/novel genome sequences were retrieved in the first step. After looking for another solution and fully understanding the capabilities of the checkm, I found it to be a much more reliable and efficient solution for finding universally present marker genes from the set of genomes being analysed.

### Refine\_cluster\_genes.sh

The purpose of the next script in the pipeline is to narrow the analysis to a smaller set of genomes to those related to one or more target genomes designated by the user. The new set of genomes selected are then more suitable for in-depth comparative genomic analysis. Genetic relatedness of the selected genomes is determined by Average Nucleotide Identity (ANIb). The genomes are annotated using PROKKA and homologue pan-core genome comparative analaysis is performed by get_homologues. 

1. Use protein BLAST (blastp) with universal marker genes from previous script to identify genomes closely related to target genome(s).
2. Perform Average Nucleotide Identity (ANI) analysis on select genomes
3. Annotate select genomes using PROKKA annotation tool.
4. Run get_homologues.pl using genbank files produced by PROKKA.
5. Perform pan-core analysis using compare_clusters.pl for nucl, prot and synt-prot homologues.
6. Build phylogenetic trees from concatenated core-gene/protein alignments.


### Genome\_func\_parallel.sh

The role of the third and final script in this pipeline is to retrieve, associate and categorise extended gene function data from the annotations produced in the previous script. This is done by retrieving Uniprot data linking annotated genes to functional data in [KEGG](https://www.genome.jp/kegg/), [COG](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/), and [GO](http://geneontology.org/) databases. This allows the gene function data from the genomes being analysed to be organised into higher order functional categories or by Gene ontology terms.

1. Get uniprot IDs from annotation data and retrieve Uniprot metadata flatfiles
2. Extract KEGG KO, COG, and GO IDs linked to Uniprot IDs and link to annotated ORFs in tabulated database.
3. Extract and link COG IDs to biological functions and COG categories.
4. Extract and tabulate KEGG biological functions, pathways and higher hierarchical categories and link to KO IDs.
5. Retrieve and link GO IDs GO terms.


## Installation and Dependencies

### Installation
Since this pipeline works via a set of bash scripts no advanced installation instructions (i.e. compilation) is required for this pipeline alone. Once the scripts are retrieved from github the only requirements are to:

First give the scripts execution privileges.

Example (on an Ubuntu system)

```
command prompt$ chmod +x get_genomes_filter_checkm.sh
```

Second add the local directory path where these scripts reside to your shells environment PATH. If using bash as default shell on an Ubuntu system for example, adding the following line to the end of your .bashrc file located in the home directory:

Example

```
export PATH=$PATH:/directory/path/to/genocap/scripts
```

_note: use 'pwd' command to get the current directory path_



### Dependencies
The scripts were written and tested in a Ubuntu server (16.04) environment with Bash shell and required command-utilities (e.g. grep, sed, sort) installed by default. The default version of awk (i.e. mawk) was replaced with [Gawk](https://www.gnu.org/software/gawk/). To take advantage of multi-CPU/core systems parallel processing was incorporated into the workflow where possible requiring [GNU parallel](https://www.gnu.org/software/parallel/) to be installed. 

As these scripts essentially run an analysis workflow each is dependent on a different set of third-party analysis software that must be installed correctly (including their own dependencies). When each script is executed by the user it will first check whether the required dependencies are installed and can be called by command.
Below are listed the third-party dependencies including those specific for each script in the pipeline. If you do use this pipeline for any published research please acknowledge and/or cite the relevant third-party software.

### Non-specific dependencies

GNU parallel - <https://www.gnu.org/software/parallel/>

RAxML - <https://cme.h-its.org/exelixis/web/software/raxml/index.html>

Trimal - <http://trimal.cgenomics.org/>

### Get\_genomes\_filter\_checkm.sh

CheckM - <https://github.com/Ecogenomics/CheckM>

Catfasta2phyml.pl - <https://github.com/nylander/catfasta2phyml>

Fasttree - <http://microbesonline.org/fasttree/>


### Refine\_cluster\_genes.sh

NCBI BLAST suite - <https://www.ncbi.nlm.nih.gov/books/NBK52640/>

PROKKA - <http://www.vicbioinformatics.com/software.prokka.shtml>

pyani - <https://github.com/widdowquinn/pyani>

Get\_homologues - <https://github.com/eead-csic-compbio/get_homologues>


### Genome\_func\_parallel.sh

No specific dependencies




## Usage


### Get_genomes_filter_checkm.sh

The first script in the pipeline requires at minimum a valid taxonomy ID/Name of prokaryotic Genus or Family and specify an output directory.


		command_prompt$ get_genomes_filter_checkm.sh -d Clostridiaceae -o Clostridiaceae_analysis


If you have genome(s) you want to include in the analysis that are not deposited in the NCBI genome repository, you can use the '-i' option to include them via the path to the directory or specifying the specific genome file (fasta format only).

		command_prompt$ get_genomes_filter_checkm.sh -d Clostridiaceae -o Clostridiaceae_analysis -i relative/path/to/my_genomes


More optional arguments are provided including full control over the genome filtering parameters. This includes (-g) and (-c) for specifying the minimum N50 and maximum number of contigs respectively for assembly quality filtering.

Three options are provided for the phylogenetic analysis of the universal markers. Fasttree (-p fast) analysis is provided as the optimal solution considering how computationally complex the analysis could become (> 500 genomes). However, RAxML (-p raxml) method is also provided for a maximum likelihood solution if you feel you have enough computing power to resolve the problem in a reasonable time frame. A third option exists (-p none) to skip the phylogenetic analysis of the genomic universal markers, if it is not required.

The full list of optional arguments for this script are listed below and can be generated using the (-h) argument on the command line.

       	-o      Specify output directory
	   	-i      Directory containing genomes in fasta format to add to the analysis (Optional)
		-t      Number of threads to use
		-d      A valid taxonomic name genus or above (e.g. Clostridiaceae)
		-g      Minimum N50 value to keep a genome (default = 20,000)
		-c      Maximum number of contigs allowed to keep a genome (default = 200)
     	-a      Minimum genome completeness as rated by checkm to keep a genome
	 	-b      Maximum genome contamination as rated by checkm to keep a genome
		-m      Maximum genome heterogeneity as rated by checkm to keep a genome
		-n      Directory containing NCBI taxonomy_dmp files (names.dmp and nodes.dmp)
		-p      Perform phylogenetic analysis on marker genes options [none|fast|raxml]

				[none]:         No phylogenetic analysis performed
				[fast]:         Construct phylogenetic tree with fasttree
				[raxml]:        Construct phylogenetic tree with RAxML (slow!)



### Refine\_cluster\_genes.sh

Once the first script (get\_genomes\_filter\_checkm.sh) has finished its analysis with no error it is safe to run this second script in the pipeline, indicating the output directory of the previous script as the input directory. 
At least one target genome needs to be specified for the analysis, whether it be one one provided by you or one retrieved by the first script. Multiple target genomes can be specified as either comma-separated list on the command line or list in a text file with each name on a new line. 


Example: targets specified in comma-separated list

		command_prompt$ refine_cluster_genes.sh -i Clostridiaceae_analysis -g target_genome_1,target_genome_2

Example: targets in text file list

		command_prompt$ refine_cluster_genes.sh -i Clostridiaceae_analysis -g target_genomes.txt

By default the results of this script will be saved in this same directory. However this can be altered using the (-o) optional argument to specify an alternative output directory. This is useful if you want to analyse different target genomes that belong to different _genera_ that reside within the same _family_ that was analysed by the first script. This helps to cut down on redundant analysis. The (-p) optional argument can also be used to change the prefix of output files produced by this script from that used by the first script.

Example: Running two separate analyses of different genera using the same input

		command_prompt$ refine_cluster_genes.sh -i Closdridiaceae_analysis -g target_clostridium_genome -o Clostridium_analysis -p Clostridium


		command_prompt$ refine_cluster_genes.sh -i Clostridiaceae_analysis -g target_caloramator_genome -o Caloramator_analysis -p Caloramator



Full list of optional arguments for refine\_cluster\_genes.sh
       
	   	-i      output directory from get_genomes_phyla_amphora.sh
		-o      Name of output directory to save analysis inside 'indir'
		-p      Output file prefix to use, default is outdir
		-g      target genome of interest or comma-separated list of genomes
		-t      number of threads/cpus to use
		-b      minimum blast identity score from marker genes to keep genome in analysis



### Genome\_func\_parallel.sh

This script performs the final portion of the analysis by linking the genome annotation data produced by the previous script to data in the Uniprot, KEGG (KO), COG and gene ontology (GO) databases. With internet access the script should be able to automatically retrieve the Uniprot flatfiles, COG and gene ontology (GO) data required. However, to connect the annotations to the KEGG (KO) functional data, a copy of the KEGG KO mapping data has to be retrieved by the user beforehand. This is done by visiting the following URL <https://www.genome.jp/kegg-bin/get_htext?ko00001>, clicking the _'Download htext'_ at the top of the page and saving the file _(ko00001.keg)_ to your system. The location of this file then needs to be indicated when running the script using the (-k) optional argument.

		command_prompt$ genome_func_parallel.sh -k /path/to/ko00001.keg -i Clostridiaceae_analysis


Full list of optional arguments for genome\_func\_parallel.sh

     	-i      Path to output directory from REFINE_CLUSTER_GENES.sh script
	 	-t      Number of cpu threads to use
		-k      Path to KEGG orthology map file

