#! /bin/bash

#: Title        : Get Genomes Filter Checkm
#: Date         : 29/10/2018
#: Author       : Joseph Adelskov
#: Version      : 1.0
#: Description  : Retrieve genomes from NCBI by taxonomy and build phylogenomic tree


usage () { printf "\n\n\tGET_GENOMES_FILTER_CHECKM\n\n\tUSAGE:\t[-o OUTDIR] [-i INDIR] [-t THREADS] [-d TAXON] [-g MIN N50] [-c MAX CONTIGS]\n\t\t[-a MIN COMPLETENESS] [-b MAX CONTAMINATION] [-m MAX HETEROGENEITY]\n\t\t[-n NCBI TAXONOMY DIR] [-h HELP]\n\n\tWhen supplied with a taxonomic name the script will find and retrieve all available genomes\n\tunder that taxon rank, filter genomes by N50, No. contigs and genome completeness and finally\n\tBuild a shared multi-gene phylogenetic tree.\n\n\t-o\tSpecify output directory\n\t-i\tDirectory containing genomes in fasta format to add to the analysis (Optional)\n\t-t\tNumber of threads to use\n\t-d\tA valid taxonomic name genus or above (e.g. Clostridiaceae)\n\t-g\tMinimum N50 value to keep a genome (default = 20,000)\n\t-c\tMaximum number of contigs allowed to keep a genome (default = 200)\n\t-a\tMinimum genome completeness as rated by checkm to keep a genome\n\t-b\tMaximum genome contamination as rated by checkm to keep a genome\n\t-m\tMaximum genome heterogeneity as rated by checkm to keep a genome\n\t-n\tDirectory containing NCBI taxonomy_dmp files (names.dmp and nodes.dmp)\n\t-p\tPerform phylogenetic analysis on marker genes options [none|fast|raxml]\n\n\t\t[none]:\t\tNo phylogenetic analysis performed\n\t\t[fast]:\t\tConstruct phylogenetic tree with fasttree\n\t\t[raxml]:\tConstruct phylogenetic tree with RAxML (slow!)\n\n\t-h\tThis help screen\n\n\n"

   exit;
	}


# exit on error (testing)

#set -e

################################################################################
## Genome assembly quality filter function
################################################################################

## genome_stat() function will filter genomes by No. contigs and N50. It will also check first if the files are gzipped and auto decompress files.

genome_stat() {
		   
	if [[ "$1" = *.gz ]]; then

		gzpi="TRUE"
		gzip -d "$1"
		genf="${1%.gz}"

	else

		genf="$1"

	fi

	fileonly="${genf##*/}"
	header="$(printf "$fileonly" | awk -F"_" '{print $1"_"$2}')"
	species="$(printf "$genf" | awk -F"/" '{gsub( /['\'';:,.]/ , ""); print $2}' | awk -F"_" '{print $1"_"$2}')" 
	ncont="$(grep -c ">" "$genf")"
	nfift="$(N50.py "$genf" | grep "N50" | awk '{print $2}')"
#	header="$(head -1 "$genf" | awk '/^>/ {gsub(">", ""); print $1}')"
	afact="$(awk -v n50="$nfift" -v cont="$ncont" 'BEGIN { printf (n50 / cont)}')"
	mkdir -p "$outdir"/sort_temp/


	if [ "$nfift" -ge "$gn50" ] && [ "$ncont" -le "$nocont" ]; then

		printf "$header\t$ncont\t$nfift\t$afact\n" >> "$outdir"/sort_temp/"$species"_temp.txt

		awk -v header="$header" '/^>/ {i++; printf(">"header"_contig%06d\n", i); next}{print}' "$genf" > "$outdir"/genome_quality_analysis/genstat_passed/"$header".fna
				
		printf ""$header" N50 = $nfift , #Contigs = $ncont , passed\n" >&2
		 
		printf "\n"$header"\t"$nfift"\t"$ncont"\tPassed" >> "$outdir"/genome_quality_analysis/"$taxname"_genstat_passed.tab

	else

		printf ""$header" N50 = $nfift , #Contigs = $ncont , failed\n" >&2
			
		printf "\n"$header"\t"$nfift"\t"$ncont"\tFailed" >> "$outdir"/genome_quality_analysis/"$taxname"_genstat_passed.tab
	fi

	if [ "$gzpi" == "TRUE" ]; then

		gzip "$genf"

	fi

}


################################################################################
## Split genomes function
################################################################################

split_genomes() {

   split_file="$(printf "$1" | awk -F"/" '{print $NF}')"

   echo "$split_file" >&2

   mkdir -p "$outdir"/genome_quality_analysis/genome_groups/"$split_file"

   while read line; do

	  # echo "copy $line to genome_groups" >&2

      cp "$outdir"/genome_quality_analysis/genstat_passed/"$line" "$outdir"/genome_quality_analysis/genome_groups/"$split_file"/

   done < "$1"

}


################################################################################
## Parallel compatible gzip function
################################################################################

par_gzip() {

	if [ "${1##*.}" == "gz" ]; then

		gzip -d $1

	else

		gzip $1

	fi

}

################################################################################
## Process checkm produced marker genes function
################################################################################

process_markers() {

	gene_id="$(printf "$1" | cut -d ' ' -f 3)"
	contig="$(printf "$1" | cut -d ' ' -f 2)"
	marker_name="$(printf "$1" | cut -d ' ' -f 4 | cut -d ';' -f 1 | sed 's/marker=//g')"

	awk -v contig="$contig" -v gene="$gene_id" 'BEGIN{RS=">";ORS="\n"}$2==contig&&$3==gene{gsub("\\*","");$2="";$3="";$4="";print ">" $1"\n"$NF}' "$outdir"/checkm_marker_analysis/"$taxname"_marker_genes.faa >> "$outdir"/checkm_marker_analysis/marker_split/"$marker_name"_tmp

}

################################################################################
## Remove duplicate marker genes function
################################################################################

remove_dupes() {

	file="${1##*/}"

	while read num acc; do

		if [ "$num" -gt "1" ]; then

			awk -v acc="$acc" 'BEGIN{RS=">"}$1 ~ acc{if(length($NF) > length(x))x=$NF}END{print ">"acc"\n"x}' "$outdir"/checkm_marker_analysis/marker_split/"$file" >> "$outdir"/checkm_marker_analysis/marker_split/"${file/_tmp/.faa}"

		else

			awk -v acc="$acc" 'BEGIN{RS=">"}$1 ~ acc{print ">"acc"\n"$NF}' "$outdir"/checkm_marker_analysis/marker_split/"$file" >> "$outdir"/checkm_marker_analysis/marker_split/"${file/_tmp/.faa}"

		fi

	done < <(awk '/^>/{gsub(">",""); print}' "$outdir"/checkm_marker_analysis/marker_split/"$file" | sort | uniq -c)

}


################################################################################
## Direct user script input to variables
################################################################################

while getopts hi:t:d:o:g:c:n:a:b:m:p: option
do
   case "${option}"
	  in
	  i) indir=${OPTARG};;
	  t) threads=${OPTARG};;
	  d) intax=${OPTARG};;
	  o) outdir=${OPTARG};;
	  g) gn50=${OPTARG};;
	  c) nocont=${OPTARG};;
	  n) ncbitax=${OPTARG};;
	  a) mincomp=${OPTARG};;
	  b) maxcont=${OPTARG};;
	  m) maxhetero=${OPTARG};;
	  p) phylopt=${OPTARG};;
	  h | *)
		 usage
		 ;;
   esac
done


# Check if required variables have a given value and/or revert to defaults

if [ -z "$intax" ]; then

   	printf "\n  Please provide a taxonomic name of genus or greater rank\n"

	usage

   	exit 1

fi


if [ -z "$outdir" ]; then

	outdir="$PWD/tax_gen_output"

	printf "\n  No output directory specified, setting output directory to $outdir\n"

fi

if [ -z "$threads" ]; then

    threads="1"
	
fi

printf "  Using $threads threads\n"

	
if [ -z "$gn50" ]; then

	gn50="20000"

fi

printf "  Minimum N50 set: $gn50\n"


if [ -z "$nocont" ]; then
	
	nocont="200"

fi

printf "  Maximum No. Contigs set: $nocont\n"

if [ -z "$mincomp" ]; then

	mincomp="90"

fi

printf "  Minimum genome clompleteness set: $mincomp\n"

if [ -z "$maxcont" ]; then

	maxcont="2"

fi

printf "  Maximum genome contamination set: $maxcont\n"

if [ -z "$maxhetero" ]; then

	maxhetero=40

fi

printf "  Maximum genome heterogeneity set: $maxhetero\n"

if [ -z "$phylopt" ]; then

	# set phylopt to default (fasttree) analysis

	printf "  Will run fasttree phylogenetic analysis as default!\n"

	phylopt="fast"

elif [ "$phylopt" == "none" ]; then

	printf "  No phylogenetic analysis of marker genes selected\n"

elif [ "$phylopt" == "fast" ]; then
	
	printf "  Will run fasttree phylogenetic analysis\n"
	
elif [ "$phylopt" != "raxml" ]; then

	printf "  Will run RAxML phylogenetic analysis\n"

else

	printf "  Please use valid option [fast|raxml|none] with -p option, exiting\n\n"

	exit

fi

# store running date/time

run_date=$(date)

date="$(date | sed 's/ /-/g')"

# export functions

export -f genome_stat
export -f split_genomes
export -f par_gzip
export -f process_markers
export -f remove_dupes

# export variables to work with functions

export outdir="$outdir"
export threads="$threads"
export nocont="$nocont"
export gn50="$gn50"


################################################################################
## Check that software dependencies are available
################################################################################

printf "\n  Checking software dependencies are installed\n"

if [ -n "$(which checkm)" ] && [ -n "$(checkm -h)" ]; then

	printf "  checkm is installed: $(which checkm)\n"

else

	printf "  checkm is not installed (correctly) !!\n"
	missdep="yes"

fi

printf "  Check for catfasta2phyml.pl script\n"

if [ -n "$(which catfasta2phyml.pl)" ]; then

	printf "  catfasta2phyml.pl is installed: $(which catfasta2phyml.pl)\n"

else

	printf "  catfasta2phyml.pl could not be found!\n"
	missdep="yes"

fi

if [ -n "$(which parallel)" ]; then

	printf "  Parallel is installed: $(which parallel)\n"

else

	printf "  Parallel is not installed (correctly?) !!\n"
	missdep="yes"

fi

if [ "$phylopt" == "fast" ]; then

	if [ -n "$(which fasttree)" ]; then

		printf "  fasttree is installed: $(which fasttree)\n"

	else 
		
		printf "  fasttree is not installed (correctly?) !!\n"
		missdep="yes"

	fi

elif [ "$phylopt" == "raxml" ]; then

	if [ -n "$(which raxmlHPC-PTHREADS)" ]; then

		printf "  RAxML is installed: $(which raxmlHPC-PTHREADS)\n"

	else

		printf "  RAxML is not installed (correctly?) !!\n"
		missdep="yes"

	fi

fi

if [ "$missdep" == "yes" ]; then

	printf "\n  One or more software dependencies cannot be called correctly\n  by this workflow script, Exiting\n\n"

	exit 1

else

	printf "  Required dependencies have been found, starting workflow\n"
fi


################################################################################
## NCBI taxonomy files
################################################################################

# Check that the given directory is correct and the required NCBI taxonomy files are present otherwise prompt user if the script should retrieve tax_dump from the NCBI ftp site

if [ -z "$ncbitax" ]; then

	if [ -f "$outdir"/NCBI_tax/nodes.dmp ] && [ -f "$outdir"/NCBI_tax/names.dmp ]; then

		printf "  NCBI files found in $outdir/NCBI_tax from previous run\n"

		ncbitax="$outdir/NCBI_tax"
		notax="0"

	else

        printf "\n  Directory path to ncbi tax dump files not supplied\n" 
        notax="1"

	fi

elif [[ ! $(ls -A "$ncbitax") ]]; then

        printf "\n  Directory path to ncbi tax dump files is empty\n"
        notax="1"

elif [ ! -f "$ncbitax"/nodes.dmp ]; then

        printf "\n  nodes.dmp file not found at $ncbitax\n"
        notax="1"

elif [ ! -f "$ncbitax"/names.dmp ]; then

        printf "\n  names.dmp file not found at $ncbitax\n"
        notax="1"
fi

if [ "$notax" = "1" ]; then

        printf "  Download NCBI tax_dump ?\n" 

        select get_taxdump in "Yes" "No"; do

                case $get_taxdump in
                        Yes )
                                mkdir -p "$outdir"/NCBI_tax
								
								if [ ! -f "$outdir"/NCBI_tax/new_taxdump.tar.gz ]; then
                                
									printf "  Retrieving NCBI new_taxdump.tar.gz\n"

									wget -q -O "$outdir"/NCBI_tax/new_taxdump.tar.gz ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz || exit 1
                                
									wget -q -O "$outdir"/NCBI_tax/new_taxdump.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5 || exit 1
                                
								fi

								while read md5code filename; do

									if [ "$(md5sum "$outdir"/NCBI_tax/"$filename" | awk '{print $1}')" = "$md5code" ]; then

                                                printf "  md5sum check passed, extracting tax dump files\n"
                                                tar -xzvf "$outdir"/NCBI_tax/new_taxdump.tar.gz -C "$outdir"/NCBI_tax
                                                ncbitax="$outdir/NCBI_tax"
												printf "\n  NCBITAX directory set to $ncbitax\n"

                                        else
											printf "\n  problem with files\n"
                                                exit 1
											
									fi

                                done < "$outdir"/NCBI_tax/new_taxdump.tar.gz.md5
								
								if [ -f "$outdir"/NCBI_tax/nodes.dmp ] && [ -f "$outdir"/NCBI_tax/names.dmp ]; then

									#printf "  Tax dump has been downloaded and extracted\n"

									break

								fi
								;;


                        No )
                                printf "\n  You have selected not to retrieve the NCBI taxonomy files, please provide a valid directory containing\nthe NCBI taxonomy nodes.dmp and names.dmp files when rerunning the script\n"

                                exit 1
								;;
						* )
								printf "\n  Please select a proper option 1 (yes) or 2 (no)\n"
								;;
                esac

        done

fi
 	
 export nodefile="$ncbitax/nodes.dmp"
 export namefile="$ncbitax/names.dmp"
 
################################################################################
## Extract taxonomy information
################################################################################

 let "iteration = 0"

## Check if taxonomy provided is name or ID and find the associated ID or name

if [ -f "$ncbitax"/names.dmp ] && [ -f "$ncbitax"/nodes.dmp ]; then


	if [[ $intax =~ ^[0-9]+$ ]]; then
	 
	    taxname="$( grep "$intax" "$namefile" | awk -v tax="$intax" -F"\t[|]\t" ' $1 == tax {print $3}')"
	
	    parenttax="$intax"
	 
	 	intaxrank="$(grep "$parenttax" "$nodefile" | awk -F "\t[|]\t" -v tax="$parenttax" '$1 == tax {print $3}')"
	
	   	printf "\n  Taxonomy ID $intax provided, known as $taxname\n"
	
	else
	 
		parenttax="$( grep "$intax" "$namefile" | awk -v tax="$intax" -F"\t[|]\t" ' $2 == tax {print $1}')"
	
		taxname="$intax"
	
	 	intaxrank="$(grep "$parenttax" "$nodefile" | awk -F "\t[|]\t" -v tax="$parenttax" '$1 == tax {print $3}')"
	 	 
	 	printf "\n  Taxonomic name $taxname provided, Taxon ID = $parenttax\n"
	 
	fi
	 
export taxname
export intaxrank

	printf "  parenttax = $parenttax\n"
	 
	## Determine the phylum of the assigned taxonomy, this will be important later for phyla amphora
	 
#	uptax="$parenttax"
	 
#	printf "  Finding phylum $taxname belongs to\n"
	 
#	until [ "$taxrank" = "class" ]; do
	 
#		taxrank="$(grep "$uptax" "$nodefile" | awk -F"\t[|]\t" -v tax="$uptax" '$1 == tax {print $3}')"
	 	
		#printf "taxrank = $taxrank\n"
		#printf "uptax = $uptax\n"
		#printf "uptaxname = $uptaxname\n"

#	 	if [ "$taxrank" = "order" ]; then
	 
#			genorder="$uptaxname"
	 
#		 	printf "Genus $taxname belongs to the Order $genorder\n"
		 
#		fi
		 
#		uptax="$(grep "$uptax" "$nodefile" | awk -F"\t[|]\t" -v tax="$uptax" '$1 == tax {print $2}')"
		 
#		uptaxname="$(grep "$uptax" "$namefile" | grep "scientific name" | awk -F"\t[|]\t" -v tax="$uptax" ' $1 == tax {print $2}')"
		    
#		if [[ $uptaxname = *"proteobacteria"* ]]; then
		 
#		 	printf "  $taxname is part of Proteobacteria, breaking at class specific $uptaxname\n"
		 
#		 	break
		 
#		fi
		 
#		let "iteration++"
		
	#	printf " $iteration..\n"
		
#	done

else

	printf "\n  Required NCBI taxonomy files not found, exiting\n\n"

	exit 1

fi

#phylum="$(grep "$uptax" "$namefile"| grep "scientific name" |  awk -F"\t[|]\t" -v tax="$uptax" '$1== tax {print $2}')"
 
#phylum="$(printf "$phylum" | awk '{print $1}')"
 
#printf "$taxname identified as belonging to the phylum $phylum\n"

#if [ -z "$phylum" ]; then

#	printf "The provided taxonomy or taxid may be invalid, exiting script\n\n\n"

#	exit 1

#fi


## make the output directory

mkdir -p "$outdir"

################################################################################
## Start script log file ## analysis_log.txt ## and add script parameters
################################################################################

dir_path="$PWD"

printf "\n\n  Log file\t$run_date\n\n  Running script get_genomes_phyla_amphora.sh\n\n\n  Script running parameters\n  --------------------------------------------------------------------------------\n  Input directory\t=\t$indir\n  Output directory\t=\t$outdir\n  Number of threads\t=\t$threads\n  Working directory\t=\t$dir_path\n\n\n  Genome selection parameters\n  --------------------------------------------------------------------------------\n  Minimum N50\t\t\t=\t$gn50 bp\n  Maximum Number of Contigs\t=\t$nocont\n  Minimum genome completeness\t=\t$mincomp %%\n  Maximum genome contamination\t=\t$maxcont %%\n  Maximum genome heterogeneity\t=\t$maxhetero %%\n\n\n  Target taxonomy\n  --------------------------------------------------------------------------------\n  Provided taxon\t=\t$intax\n  Taxonomy name\t\t=\t$taxname\n  Taxon rank\t\t=\t$intaxrank\n  Phylum\t\t=\t$phylum\n" >> "$outdir"/get_genome_log_"$date".txt


################################################################################
## Get list of genera under selected taxon
################################################################################


## Check first if the input taxonomy is rank genus otherwise find child genera

if [ ! -f "$outdir"/genus_list.txt ]; then

	if [ "$intaxrank" = "genus" ]; then

   		printf "%s\n" "$taxname" > "$outdir"/genus_list.txt

   		printf "  Input taxonomy is a genus\n"

	else

		printf "  Find all genera under $taxname\n"

		grep "$parenttax" "$nodefile" | awk -F"\t[|]\t" -v ptax="$parenttax" ' $2 == ptax {print $1"\t"$3}' > "$outdir"/tax_list.txt	# Get a list of the child taxon IDs under the parent taxon id


	## Until loop required if a higher than 'Family' rank taxon was specified. Will continue until the taxon IDs for child 'genera' is produced


	until [ ! -s "$outdir"/tax_list.txt ]; do

   		while read taxid taxr; do

			if [ "$taxr" = "genus" ]; then

				grep "$taxid" "$namefile" | awk -F"\t[|]\t" -v tax="$taxid" ' $1 == tax {print $2"\t"$NF}' | grep "name" | awk '{print $1}' | sed -e '/unclassified/d' -e '/environmental/d' -e '/Candidatus/d' >> "$outdir"/genus_list.txt

	  		elif [ "$taxr" = "no rank" ] || [ "$taxr" = "species" ]; then

				continue

	  		else

				awk -F"\t[|]\t" -v tax="$taxid" ' $2 == tax {print $1"\t"$3}' $nodefile >> "$outdir"/new_list.txt

	  		fi

		done < "$outdir"/tax_list.txt

		rm "$outdir"/tax_list.txt

		if [ -f "$outdir"/new_list.txt ]; then

			mv "$outdir"/new_list.txt "$outdir"/tax_list.txt

 		fi
	
	done

	nogenera="$(wc -l $outdir/genus_list.txt | awk '{print $1}')"

	printf "%s" "$taxname" >> "$outdir"/genus_list.txt


	printf "  $nogenera genera found\n"

	fi

else

	nogenera="$(wc -l $outdir/genus_list.txt | awk '{print $1}')"

fi

## Add no. of genera information to log file

printf "  No. of genera\t\t=\t$nogenera\n\n  The $intaxrank $taxname was found to contain $nogenera genera.\n\n  A full list of the different genera is stored in the file: genus_list.txt\n\n  PhylAmphora gene marker set to $phylum\n\n\n" >> "$outdir"/get_genome_log_"$date".txt


################################################################################
## Retrieve genomes from NCBI
################################################################################

# retrieve the genomes from ncbi using rsync. The while loop will run through the file genus_list.txt for each genus name that was found under the higher taxonomic rank specified.

mkdir -p "$outdir"/genome_quality_analysis

if [[ ! $(ls -A "$outdir"/genome_quality_analysis/retrieved_genomes/ 2>/dev/null) ]]; then

mkdir -p "$outdir"/genome_quality_analysis/retrieved_genomes

	printf "\n  Retrieving genomes from NCBI for $taxname\n"

while read genus; do

	printf "  Retrieving $genus genomes\t"

	curr_gen_no="$(ls "$outdir"/genome_quality_analysis/retrieved_genomes/ | wc -l)"

	rsync --quiet --no-motd --msgs2stderr --exclude="*_cds_from_genomic.fna.gz" --exclude="*_rna_from_genomic.fna.gz"  ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/${genus}_*/latest_assembly_versions/*/*_genomic.fna.gz "$outdir"/genome_quality_analysis/retrieved_genomes/ 2>/dev/null || printf "No genomes found!\n"

	new_gen_no="$(ls "$outdir"/genome_quality_analysis/retrieved_genomes/ | wc -l)"

	if [ "$curr_gen_no" != "$new_gen_no" ]; then
	
		printf "Total genomes = %d\n" "$(ls "$outdir"/genome_quality_analysis/retrieved_genomes/ | wc -l)"

	fi

done < "$outdir"/genus_list.txt

fi

if [[ ! $(ls -A "$outdir"/genome_quality_analysis/retrieved_genomes/ 2>/dev/null) ]]; then

	printf "\n  No genomes found, exiting script\n\n"

	exit

else
	
	downloaded_genomes="$(ls "$outdir"/genome_quality_analysis/retrieved_genomes/ | wc -l)"

fi


################################################################################
## Build info files on retrieved genomes
################################################################################

mkdir -p "$outdir"/genome_information

if [ ! -f "$outdir"/genome_information/"$taxname"_name_key.txt ] || [ ! -f "$outdir"/genome_information/"$taxname"_genome_info.tab ]; then

	if [ -d "$outdir"/genome_quality_analysis/retrieved_genomes ]; then
	

		if [[ $(ls -A "$outdir"/genome_quality_analysis/retrieved_genomes/*.gz 2>/dev/null) ]]; then

			printf "\n  Decompressing genome files\n"

			find "$outdir"/genome_quality_analysis/retrieved_genomes/ | grep -e \.gz | parallel -P "$threads" --no-notice par_gzip

		fi

		printf "  Building "$taxname"_name_key.txt file\n"

		for file in "$outdir"/genome_quality_analysis/retrieved_genomes/*.fna ; do

	    	gen_acc="$(printf "${file##*/}" | awk -F"_" '{print $1"_"$2}')"
			head -1 "$file" | awk -F"," '{print $1}' | sed -e 's/>//g' -e 's/ $//g' -e 's/^ //g' -e 's/:/-/g' -e 's/\S*scaffold\S*//g' -e 's/\S*contig\S*//g' -e 's/\(\[\|\]\|{\|}\||\)//g' -e "s/\((\|)\|'\)//g" | awk -v acc="$gen_acc" '{$1 = ""; print acc"\t"$0}' >> "$outdir"/genome_information/"$taxname"_name_key.txt

		done

		## New way to get more NCBI genome information

		# if the 'prokaryotes.txt' file is not present, retrieve it from NCBI ftp

		if [ ! -f "$outdir"/prokaryotes.txt ]; then

			printf "  Retrieving the prokaryotes.txt file from NCBI ftp server\n"

			wget -q -O "$outdir"/prokaryotes.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
	
		fi

		for file in "$outdir"/genome_quality_analysis/retrieved_genomes/*.fna ; do

			gca="$(printf "${file##*/}" | awk -F"_" '{print $1"_"$2}')"

			genome_acc+=("$gca")

		done

		if [ ! -f "$outdir"/"$taxname"_acc_list.txt ]; then

			printf "%s\n" "${genome_acc[@]}" > "$outdir"/"$taxname"_acc_list.txt

		fi

		# loop through the genome GCA accessions to gather useful information from prokaryote.txt and save to a tab delimited file

		if [ ! -f "$outdir"/genome_information/"$taxname"_genome_char.tab ] || [ ! -f "$outdir"/genome_information/"$taxname"_genome_info.tab ]; then

			printf "\n  Building genome information and general characteristics files\n"

			printf "Accession\tName\tGroup\tSubgroup\ttaxid\tSize\tGC%%\tContigs\tGenes\tProteins\n" > "$outdir"/genome_information/"$taxname"_genome_char.tab
		
			printf "Accession\tName\tBioproject Accession\tBiosample Accession\tCentre\tRelease Date\tModify Date\tStatus\tPubmed ID\n" > "$outdir"/genome_information/"$taxname"_genome_info.tab

			for i in "${!genome_acc[@]}"; do

				awk -F"\t" -v gca="${genome_acc[$i]}" '$19 ~ gca {print gca"\t"$1"\t"$5"\t"$6"\t"$2"\t"$7"\t"$8"\t"$11"\t"$12"\t"$13}' "$outdir"/prokaryotes.txt >> "$outdir"/genome_information/"$taxname"_genome_char.tab

				awk -F"\t" -v gca="${genome_acc[$i]}" '$19 ~ gca {print gca"\t"$1"\t"$3"\t"$18"\t"$17"\t"$14"\t"$15"\t"$16"\t"$22}' "$outdir"/prokaryotes.txt >> "$outdir"/genome_information/"$taxname"_genome_info.tab

			done

		fi

		printf "\n  Recompressing genome files\n"

		find "$outdir"/genome_quality_analysis/retrieved_genomes/ | grep -e \.fna  | parallel -P "$threads" --no-notice par_gzip

	fi

fi
## Add retrieved genome information to log file

printf "  Retrieval of genomes from NCBI\n  --------------------------------------------------------------------------------\n  The NCBI ftp genome archive \"ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria\" was\n  searched using the list of genera and $taxname to retrieve relevant genomes.\n\n  No. of genomes retrieved\t=\t$downloaded_genomes\n\n  Dir: "$outdir"/genome_quality_analysis/retrieved_genomes\n\n\n" >> "$outdir"/get_genome_log_"$date".txt

tail -n +2 "$outdir"/genome_information/"$taxname"_genome_info.tab | awk '{print $2}' | sort |uniq -c | awk 'BEGIN{printf "       Genomes retrieved per Genus\n     -------------------------------\n"}{printf "     %20s  =  %d\n", ($2), ($1)}END{printf "\n\n"}' >> "$outdir"/get_genome_log_"$date".txt

# testing control exit

#exit

################################################################################
## Filter genomes by assembly quality
################################################################################

extension=".fna"

if [[ ! $(ls -A "$outdir"/genome_quality_analysis/genstat_passed/ 2>/dev/null) ]]; then

		printf "Genome\tN50\tContigs\tPass/Fail\n" > "$outdir"/genome_quality_analysis/"$taxname"_genstat_passed.tab

	mkdir -p "$outdir"/genome_quality_analysis/genstat_passed

	find "$outdir"/genome_quality_analysis/retrieved_genomes/ | grep -e \\$extension -e \.gz | parallel -P "$threads" --no-notice genome_stat

fi

genstat_passed_genomes="$(ls "$outdir"/genome_quality_analysis/genstat_passed/ 2>/dev/null | wc -l 2>/dev/null)"

## Add genome quality filtering info to log file

printf "  Genome assembly quality filtering\n  --------------------------------------------------------------------------------\n  Genomes were filtered for acceptable assembly quality. To pass filtering genomes\n  must have an N50 larger than $gn50 and consist of less than $nocont contigs.\n\n  Genstat passed genomes\t=\t$genstat_passed_genomes\n\n\n" >> "$outdir"/get_genome_log_"$date".txt


################################################################################
## Filter genomes using Checkm analysis
################################################################################

## Split the genomes into groups of 500 before running checkm

mkdir -p "$outdir"/tmp_gen_list

if [ ! -d "$outdir"/genome_quality_analysis/genome_groups ]; then

	mkdir -p "$outdir"/genome_quality_analysis/genome_groups

	ls -R "$outdir"/genome_quality_analysis/genstat_passed/ | grep ".fna" > "$outdir"/genome_quality_analysis/sel_gen_list.txt

	split -l 500 -a 4 -d "$outdir"/genome_quality_analysis/sel_gen_list.txt "$outdir"/tmp_gen_list/gen_grp_

	find "$outdir"/tmp_gen_list/ | tail -n +2 | parallel -P "$threads" --no-notice split_genomes

fi

# remove tmp_gen_list dir before checkm

if [ -d "$outdir"/tmp_gen_list ]; then

	rm -r "$outdir"/tmp_gen_list

fi

## Run checkm on the genomes to check the quality

if [[ ! $(ls -A "$outdir"/quality_genomes/ 2>/dev/null) ]]; then

	for file in "$outdir"/genome_quality_analysis/genome_groups/*/ ; do

		mkdir -p "$outdir"/genome_quality_analysis/checkm_lwf

    	grp_dir="$(printf "$file" | awk -F"/" '{print $(NF-1)}')"

		if [ ! -f "$outdir"/genome_quality_analysis/checkm_lwf/"$grp_dir"/lineage_wf.tab ]; then

			mkdir -p "$outdir"/genome_quality_analysis/checkm_lwf/"$grp_dir"

    		checkm lineage_wf -t "$threads" --pplacer_threads "$threads" -x "$extension" -f "$outdir"/genome_quality_analysis/checkm_lwf/"$grp_dir"/lineage_wf.tab --tab_table "$file" "$outdir"/genome_quality_analysis/checkm_lwf/"$grp_dir" || exit 1

		fi

		qual_genome_list+=($(awk -v gen_comp="$mincomp" -v gen_cont="$maxcont" -v gen_hetero="$maxhetero" -F"\t" 'NR>1&&$2!~/k__Bacteria/&&$(NF-2)>=gen_comp&&$(NF-1)<=gen_cont&&$NF<=gen_hetero{print $1}' "$outdir"/genome_quality_analysis/checkm_lwf/"$grp_dir"/lineage_wf.tab))


	done

    mkdir -p "$outdir"/quality_genomes

	for i in ${!qual_genome_list[@]}; do

		cp "$outdir"/genome_quality_analysis/genstat_passed/"${qual_genome_list[$i]}".fna "$outdir"/quality_genomes/

	done

fi

if [ -z "$qual_genome_list" ]; then

	for file in "$outdir"/genome_quality_analysis/genome_groups/*/ ; do
		
		grp_dir="$(printf "$file" | awk -F"/" '{print $(NF-1)}')"

		qual_genome_list+=($(awk -v gen_comp="$min_comp" -v gen_cont="$maxcont" -v gen_hetero="$maxhetero" -F"\t" 'NR>1&&$2!~/k__Bacteria/&&$(NF-2)>=gen_comp&&$(NF-1)<gen_cont&&$NF<gen_hetero{print $1}' "$outdir"/genome_quality_analysis/checkm_lwf/"$grp_dir"/lineage_wf.tab))


	done

fi

# Extract genome information on only remaining quality filtered genomes
printf "Accession\tName\tBioproject Accession\tBiosample Accession\tCentre\tRelease Date\tModify Date\tStatus\tPubmed ID\n" > "$outdir"/genome_information/"$taxname"_quality_genome_info.tab

printf "Accession\tName\tGroup\tSubgroup\ttaxid\tSize\tGC%%\tContigs\tGenes\tProteins\n" > "$outdir"/genome_information/"$taxname"_quality_genome_char.tab

for i in ${!qual_genome_list[@]}; do


	grep "^"${qual_genome_list[$i]}"" "$outdir"/genome_information/"$taxname"_genome_info.tab >> "$outdir"/genome_information/"$taxname"_quality_genome_info.tab

	grep "^"${qual_genome_list[$i]}"" "$outdir"/genome_information/"$taxname"_genome_char.tab >> "$outdir"/genome_information/"$taxname"_quality_genome_char.tab

done


## Accession and Name only from genome_char.tab file to new file

if [ ! -f "$outdir"/genome_information/"$taxname"_tax_naming.txt ]; then

	for  file in "$outdir"/quality_genomes/*.fna; do

		genome_file="${file##*/}"

		grep "${genome_file/.fna/}" "$outdir"/genome_information/"$taxname"_genome_char.tab | cut -f 1,2 -s -d$'\t' | sed -e 's/\(\[\|\]\|{\|}\||\)//g' -e "s/\((\|)\|'\)//g" -e 's/:/-/g' >> "$outdir"/genome_information/"$taxname"_tax_naming.txt

	done
	
	#awk -F"\t" -v acc=

fi

## Track number of checkm filtered genomes only

no_quality_genomes="$(ls "$outdir"/quality_genomes/ 2>/dev/null | wc -l 2>/dev/null)"

## Add the user selected genomes to the quality genomes

if [ ! -z "$indir" ]; then

	cp "$indir"/*.fna "$outdir"/quality_genomes/

	printf "  Adding genomes from "$indir" to quality selected genomes\n"

fi

# Count number of remaining genomes

gencount="$(ls "$outdir"/quality_genomes/ | wc -l | awk '{print $1}')"

printf "  $gencount high quality genomes selected\n"


## Add post-checkm filtering information


printf "  Checkm genome filtering\n  --------------------------------------------------------------------------------\n  The genomes that passed assembly quality filtering were analysed using checkm to\n  determine the genome completeness, contamination and heterogeneity. A minimum\n  genome completeness of $mincomp %%, maximum contamination of $maxcont %%, and maximum heterogeneity\n  of $maxhetero %% was required to pass filtering. The genomes that passed this fitering were\n  stored in the \""$outdir"/quality_genomes/\" directory. Following filtering\n  genomes supplied via input \""$indir"\" were added to this directory.\n\n  \tCheckm passed genomes\t\t=\t$no_quality_genomes\n\n  \tInput genomes added\t\t=\t$gencount\n\n\n" >> "$outdir"/get_genome_log_"$date".txt

tail -n +2 "$outdir"/genome_information/"$taxname"_quality_genome_info.tab | awk '{print $2}' | sort |uniq -c | awk 'BEGIN{printf "       Quality Genomes per Genus\n     -------------------------------\n"}{printf "     %20s  =  %d\n", ($2), ($1)}END{printf "\n\n"}' >> "$outdir"/get_genome_log_"$date".txt


################################################################################
## Get full taxonomic nomenclature for quality genomes using NCBI taxonomy
################################################################################

#rm "$outdir"/genome_information/"$taxname"_tax_nomenclature.txt 2>/dev/null

if [ ! -f "$outdir"/genome_information/"$taxname"_tax_nomenclature.txt ]; then

	printf "  Building tax_nomenclature file\n"

	while read gen_info; do

		genometax="$(printf "$gen_info" | awk -F"\t" '{print $5}')"
		accession="$(printf "$gen_info" | awk -F"\t" '{print $1}')"

			taxstring="$( grep "$genometax" "$namefile" | grep "strain\|sp." | awk -F"\t[|]\t" -v gtax="$genometax" '$1 == gtax {print $2}')"

			if [ -z "$taxstring" ]; then

				taxstring="$( grep "$genometax" "$namefile" | grep "scientific name" | awk -F"\t[|]\t" -v gtax="$genometax" '$1 == gtax {print $2}')"

				if [ -z "$taxstring" ]; then

					taxstring="$( grep "$genometax" "$namefile" | grep "authority" | awk -F"\t[|]\t" -v gtax="$genometax" '$1 == gtax {print $2}')"

				fi

			fi

			# delete characters like parentheses '()' that will screw up newick formats when renaming the nodes

			taxstring="$(printf "$taxstring" | sed -e 's/\(\[\|\]\|{\|}\||\|(\|)\)//g' | tr -d '\n')"

			printf "$accession\t$taxstring\n" >> "$outdir"/genome_information/"$taxname"_tax_nomenclature.txt

		done < <(tail -n +2 "$outdir"/genome_information/"$taxname"_genome_char.tab)

fi

## Add genome info to log file

printf "  Extracting genome information\n  --------------------------------------------------------------------------------\n  Relevant genomic information (Metadata) was extracted for each of the remaining\n  high quality genomes. This information is stored in several files in the\n  directory: $outdir/genome_information/\n\n  Genome Naming key\t\t=\t"$taxname"_name_key.txt\n  Genome characteristics\t=\t"$taxname"_genome_char.tab\n  Genome information\t\t=\t"$taxname"_genome_info.tab\n  Genome tax naming\t\t=\t"$taxname"_tax_naming.txt\n  Taxonomic nomenclature\t=\t"$taxname"_tax_nomenclature\n\n  The tax naming file represents full genome identification as per NCBI genome\n  prokaryotes.txt file, where as taxonomic nomenclature naming is retreived from\n  NCBI taxonomy files\n\n\n" >> "$outdir"/get_genome_log_"$date".txt

################################################################################
## Checkm marker analysis
################################################################################

if [ ! -f "$outdir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta ]; then

## Use checkm to get marker genes from the quality genomes

mkdir -p "$outdir"/checkm_marker_analysis

if [ ! -f "$outdir"/checkm_marker_analysis/"$taxname"_markers.mrk ]; then

	checkm taxon_set "$intaxrank" "$taxname" "$outdir"/checkm_marker_analysis/"$taxname"_markers.mrk

fi

if [[ ! $(ls -A "$outdir"/checkm_marker_analysis/find_markers 2>/dev/null ) ]]; then

	checkm analyze -t "$threads" "$outdir"/checkm_marker_analysis/"$taxname"_markers.mrk "$outdir"/quality_genomes "$outdir"/checkm_marker_analysis/find_markers

fi

if [ ! -f "$outdir"/checkm_marker_analysis/"$taxname"_marker_genes.faa ]; then

	checkm qa -t "$threads" -o 9 -f "$outdir"/checkm_marker_analysis/"$taxname"_marker_genes.faa "$outdir"/checkm_marker_analysis/"$taxname"_markers.mrk "$outdir"/checkm_marker_analysis/find_markers

fi

# clean up marker_split output if stopped abruptly in a prevous run

if [ ! -d "$outdir"/checkm_marker_analysis/marker_split ]; then

	mkdir -p "$outdir"/checkm_marker_analysis/marker_split

else

	total_markers="$(grep -c ">" "$outdir"/checkm_marker_analysis/"$taxname"_marker_genes.faa)"
	total_tmp_markers="$(grep -c ">" "$outdir"/checkm_marker_analysis/marker_split/*_tmp 2>/dev/null | awk -F":" '{sum+=$2}END{print sum}')"

	if [ "$total_tmp_markers" != "$total_markers" ];then

		rm  "$outdir"/checkm_marker_analysis/marker_split/*.tmp

	else

		printf "  All markers accounted for in marker_tmp files\n"

	fi

fi

# split the combined marker file output from checkm into individual files for each different marker

if [[ ! $(ls -A "$outdir"/checkm_marker_analysis/marker_split/*_tmp 2>/dev/null) ]]; then

	printf "  \n\nSplitting marker genes to individual marker files\n"

	awk '/^>/{gsub(">","");print}' "$outdir"/checkm_marker_analysis/"$taxname"_marker_genes.faa | parallel -P "$threads" --no-notice process_markers

fi

# if a genome has two instances of the same marker gene pick the longest of the two or more

rm "$outdir"/checkm_marker_analysis/marker_split/*.faa 2>/dev/null

if [[ ! $(ls -A "$outdir"/checkm_marker_analysis/marker_split/*.faa 2>/dev/null) ]]; then

	printf "  Removing excess markers from the same genome\n"

	find "$outdir"/checkm_marker_analysis/marker_split | grep -e _tmp | parallel -P "$threads" --no-notice remove_dupes

fi

mkdir -p "$outdir"/checkm_marker_analysis/good_markers
mkdir -p "$outdir"/checkm_marker_analysis/aligned_markers


# select marker files with one gene for each quality genome, align and trim

if [[ ! $(ls -A "$outdir"/checkm_markers_analysis/aligned_markers/*_aln_trm.faa 2>/dev/null) ]]; then

	while read marker count; do

		if [ "$count" -eq "$gencount" ]; then

			marker_file="${marker##*/}"
		
			good_markers+=("$marker_file")

			cp "$marker" "$outdir"/checkm_marker_analysis/good_markers

			if [ ! -f "$outdir"/checkm_marker_analysis/aligned_markers/"${marker_file/.faa/_aln_trm.faa}" ]; then

				printf "  Aligning $marker_file with muscle\n"

				muscle -quiet -in "$marker" -out "$outdir"/checkm_marker_analysis/aligned_markers/"${marker_file/.faa/_aln.faa}"

				trimal -in "$outdir"/checkm_marker_analysis/aligned_markers/"${marker_file/.faa/_aln.faa}" -out "$outdir"/checkm_marker_analysis/aligned_markers/"${marker_file/.faa/_aln_trm.faa}" -gt 0.8 -st 0.001

			fi

		fi

	done < <(grep -c ">" "$outdir"/checkm_marker_analysis/marker_split/*.faa | sed 's/:/ /g')

fi

printf "  %d acceptable marker genes identified\n" "${#good_markers[@]}"

# use the script catfasta2phyml.pl to merge the alignment files into one

if [ ! -f "$outdir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta ]; then

	printf "  Concatentating the %d aligned marker genes together\n" "${#good_markers[@]}"

	catfasta2phyml.pl -f "$outdir"/checkm_marker_analysis/aligned_markers/*_trm.faa > "$outdir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta

fi

fi

## Add checkm marker analysis to log

printf "  Checkm gene marker analysis\n  --------------------------------------------------------------------------------\n  Checkm genome quality commands were used to identify and extract marker genes from\n  the quality genomes for the tax = "$taxname", rank = "$intaxrank". The markers are extracted by checkm into a\n  single fasta file by default. Therefore the file was split up into individual files\n  for each marker gene. For marker gene files which contained multiple genes from the same genome, the longest of\n  these genes was selected and the other removed. Marker gene files which 1 gene\n  sequence from each of the quality genomes were selected, aligned using muscle alignment program. Highly\n  varient regions (less than 80% genome representation) of the alignments were trimmed\n  using trimal and the trimmed alignments were then concatenated using the catfasta2phyml.pl script.\n\n" >> "$outdir"/get_genome_log_"$date".txt


grep ">" -c "$outdir"/checkm_marker_analysis/marker_split/*_tmp | awk -F":" '{sum = sum + $2; if($2 > max) max = $2; if($2<min||NR==1)min = $2}END{print "\n  \tAll Marker Statistics:\n\n  \tNo. marker genes\t=\t"NR"\n  \tTotal genes\t\t=\t"sum"\n  \tMean\t\t\t=\t"sum/NR"\n  \tMin\t\t\t=\t"min"\n  \tMax\t\t\t=\t"max"\n  \tRange\t\t\t=\t"max-min}' >> "$outdir"/get_genome_log_"$date".txt

printf "\n\n  \tGood Marker Statistics:\n\n" >> "$outdir"/get_genome_log_"$date".txt

grep ">" -c "$outdir"/checkm_marker_analysis/good_markers/*.faa | awk -F":" '{sum = sum + $2; if($2 > max) max = $2; if($2<min||NR==1)min = $2}END{print "  \tNo. marker genes\t=\t"NR"\n  \tTotal genes\t\t=\t"sum"\n  \tMean\t\t\t=\t"sum/NR"\n  \tMin\t\t\t=\t"min"\n  \tMax\t\t\t=\t"max"\n  \tRange\t\t\t=\t"max-min}' >> "$outdir"/get_genome_log_"$date".txt

awk '/^>/{i++}!/>/{if(i>1) exit; print}' "$outdir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta | awk '{seqlen += length($0)}END{print "  \tAlignment length\t=\t"seqlen" AA\n\n\n"}' >> "$outdir"/get_genome_log_"$date".txt

awk 'BEGIN{RS=">";ORS=""}{$1 =""; gsub(" ",""); print $0}' "$outdir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta | sed -e 's/-//g' -e 's/./&\n/g' | sort | uniq -c | awk 'BEGIN{print "  \tAlignment composition\n  \t---------------------\n"}{AA[$i] = $2; count[$i] = $1}{sum += $1}END{for (i in AA) {printf "           %s  =  %.3f %\n", (AA[i]), (count[i]/sum*100)}; print "\n\n\n"}' >> "$outdir"/get_genome_log_"$date".txt


################################################################################
## Marker Phylognetic analysis
################################################################################

if [ "$phylopt" == "fast" ]; then


################################################################################
## FastTree Marker phylogenetic analysis
################################################################################

	if [[ ! $( ls -A "$outdir"/checkm_marker_analysis/fasttree_analysis 2>/dev/null) ]]; then

		printf "\nPerforming phylogenetic analysis on marker alignment using FastTree\n"

		mkdir -p "$outdir"/checkm_marker_analysis/fasttree_analysis

		fasttree "$outdir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta > "$outdir"/checkm_marker_analysis/fasttree_analysis/"$taxname"_fasttree.nwk

	fi


	if [ -f "$outdir"/checkm_marker_analysis/fasttree_analysis/"$taxname"_fasttree.nwk ]; then

		mkdir -p "$outdir"/checkm_marker_analysis/fasttree_analysis/renamed_tree

		rename_nodes.sh -i "$outdir"/checkm_marker_analysis/fasttree_analysis/"$taxname"_fasttree.nwk -n "$outdir"/genome_information/"$taxname"_tax_naming.txt -o "$outdir"/checkm_marker_analysis/fasttree_analysis/renamed_tree/"$taxname"_fasttree_tname.nwk
		
		rename_nodes.sh -i "$outdir"/checkm_marker_analysis/fasttree_analysis/"$taxname"_fasttree.nwk -n "$outdir"/genome_information/"$taxname"_name_key.txt -o "$outdir"/checkm_marker_analysis/fasttree_analysis/renamed_tree/"$taxname"_fasttree_fname.nwk


	fi

## Add FastTree analysis info to log file

	printf "  FastTree phylogenetic analysis\n  --------------------------------------------------------------------------------\n  FastTree was used to perform a phylogenetic analysis on the concatenated alignment.\n  The parameters used in the phylogenetic analysis are listed below:\n\n  \tAlignment\t\t=\tAmino acid\n  \tMethod\t\t\t=\tHeuristic neighbor-joining\n  \tEvolutionary model\t=\tGamma distritubion\n  \tSubstitution matrix\t=\tJones, Taylor & Thornton (JTT)\n  \tTest of phylogeny\t=\tNone\n  \tReplicates\t\t=\tNA\n\n\n  Following the FastTree analysis the tips of the phylogenetic tree produced were\n  renamed using either "$taxname"_name_key.txt or "$taxname"_tax_naming.txt files.\n\n" >> "$outdir"/get_genome_log_"$date".txt 

	printf "  Renamed trees stored in "$outdir"/checkm_marker_analysis/fasttree_analysis/renamed_tree\n\n  "$taxname"_name_key.txt: "$taxname"_fasttree_fname.nwk\n  "$taxname"_tax_naming.txt: "$taxname"_fasttree_tname.nwk\n\n\n" >> "$outdir"/get_genome_log_"$date".txt

elif [ "$phylopt" == "raxml" ]; then

################################################################################
## run RAxML on the concatenated alignment
################################################################################

	if [[ ! $( ls -A "$outdir"/checkm_marker_analysis/raxml_analysis 2>/dev/null) ]]; then

		printf "\nPerforming RAXML maximum likelyhood phylogenetic analysis"
	
		mkdir -p "$outdir"/checkm_marker_analysis/raxml_analysis

		raxmlHPC-PTHREADS -f a -s "$outdir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta -x 127 -# 100 -m PROTGAMMAJTT -n "$taxname"_100b -p 129 -T "$threads"

		mv RAxML* "$outdir"/checkm_marker_analysis/raxml_analysis/

	fi


## Rename the nodes of the phylogenetic trees produced by RAxML using the naming lists

	if [[ ! $(ls -A "$outdir"/checkm_marker_analysis/raxml_analysis/renamed_tree 2>/dev/null) ]] && [[ $(ls -A "$outdir"/checkm_marker_analysis/raxml_analysis 2>/dev/null) ]]; then


		printf "\nRenaming nodes of ML trees\n"

		mkdir -p "$outdir"/checkm_marker_analysis/raxml_analysis/renamed_tree

		export taxname="$taxname"

		rename_nodes.sh -i "$PWD"/"$outdir"/checkm_marker_analysis/raxml_analysis/RAxML_bestTree."$taxname"_100b -n "$PWD"/"$outdir"/genome_information/"$taxname"_tax_name_key.txt -o "$PWD"/"$outdir"/checkm_marker_analysis/raxml_analysis/renamed_tree/"$taxname"_MLtree_100bs_tax_fname_no_bs_labels.nwk

		rename_nodes.sh -i "$outdir"/checkm_marker_analysis/raxml_analysis/RAxML_bipartitions."$taxname"_100b -n "$outdir"/genome_information/"$taxname"_name_key.txt -o "$outdir"/checkm_marker_analysis/raxml_analysis/renamed_tree/"$taxname"_MLtree_100bs_fnamed_with_bs_labels.nwk

		rename_nodes.sh -i "$PWD"/"$outdir"/checkm_marker_analysis/raxml_analysis/RAxML_bestTree."$taxname"_100b -n "$PWD"/"$outdir"/genome_information/"$taxname"_tax_naming.txt -o "$PWD"/"$outdir"/checkm_marker_analysis/raxml_analysis/renamed_tree/"$taxname"_MLtree_100bs_tnamed_no_bs_labels.nwk
		
		rename_nodes.sh -i "$outdir"/checkm_marker_analysis/raxml_analysis/RAxML_bipartitions."$taxname"_100b -n "$outdir"/genome_information/"$taxname"_tax_naming.txt -o "$outdir"/checkm_marker_analysis/raxml_analysis/renamed_tree/"$taxname"_MLtree_100bs_tnamed_with_bs_labels.nwk

	fi

## Add RAXML analysis info to log file

	printf "  RAXML maximum likelihood phylogenetic analysis\n  --------------------------------------------------------------------------------\n  RAXML was used to perform a phylogenetic analysis on the concatenated alignment.\n  The parameters used in the phylogenetic analysis are listed below:\n\n  \tAlignment\t\t=\tAmino acid\n  \tMethod\t\t\t=\tMaximum Likelihood\n  \tEvolutionary model\t=\tGamma distritubion\n  \tSubstitution matrix\t=\tJones, Taylor & Thornton (JTT)\n  \tTest of phylogeny\t=\tBootstrap\n  \tReplicates\t\t=\t100\n\n\n  Following the RAXML analysis the tips of the phylogenetic tree produced were\n  renamed using either "$taxname"_name_key.txt or "$taxname"_tax_naming.txt files.\n\n" >> "$outdir"/get_genome_log_"$date".txt 

	printf "  Renamed trees stored in "$outdir"/checkm_marker_analysis/raxml_analysis/renamed_tree\n\n  ML tree no bootstrap renamed "$taxname"_name_key.txt: "$taxname"_MLtree_100bs_fnamed_no_bs_labels.nwk\n  ML tree with bootstrap renamed "$taxname"_name_key.txt: "$taxname"_MLtree_100bs_fnamed_with_bs_labels.nwk\n  ML tree no bootstrap renamed "$taxname"_tax_naming.txt: "$taxname"_MLtree_100bs_tnamed_no_bs_labels.nwk\n  ML tree with bootstrap renamed "$taxname"_tax_naming.txt: "$taxname"_MLtree_100bs_tnamed_with_bs_labels.nwk\n\n" >> "$outdir"/get_genome_log_"$date".txt

else

	printf "\nNo phylogenetic analysis on markers performed as requested\n"

fi

cat "$outdir"/get_genome_log_"$date".txt > "$outdir"/"$taxname"_analysis_report.txt


exit


