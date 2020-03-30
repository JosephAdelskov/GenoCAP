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
Since this pipeline works via a set of bash scripts no advanced installation instructions (i.e. compilation) is required for this pipeline alone. Once the scripts are retrieved from github the only requirement is to add the local directory path where these scripts reside to your shells environment PATH. If using bash as default shell on an Ubuntu system for example, adding the following line to the end of your .bashrc file located in the home directory:

(Example)

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
