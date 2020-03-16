#! /bin/bash

#: Title        : Get Genomes Phyla Amphora
#: Date         : 05/05/2018
#: Author       : Joseph Adelskov
#: Version      : 1.0
#: Description  : Retrieve genomes from NCBI by taxonomy and build phylogenomic tree


usage () { printf "\n\n\tGET_GENOMES_PHYLA_AMPHORA\n\n\tUSAGE:\t[-o OUTDIR] [-i INDIR] [-t THREADS] [-d TAXON] [-g MIN N50] [-c MAX CONTIGS]\n\t\t[-a MIN COMPLETENESS] [-b MAX CONTAMINATION] [-m MAX HETEROGENEITY]\n\t\t[-n NCBI TAXONOMY DIR] [-q MARKER DATABASE ] [-h HELP]\n\n\tWhen supplied with a taxonomic name the script will find and retrieve all available genomes\n\tunder that taxon rank, filter genomes by N50, No. contigs and genome completeness and finally\n\tBuild a shared multi-gene phylogenetic tree.\n\n\t-o\tSpecify output directory\n\t-i\tDirectory containing genomes in fasta format to add to the analysis (Optional)\n\t-t\tNumber of threads to use\n\t-d\tA valid taxonomic name genus or above (e.g. Clostridiaceae)\n\t-g\tMinimum N50 value to keep a genome (default = 20,000)\n\t-c\tMaximum number of contigs allowed to keep a genome (default = 200)\n\t-a\tMinimum genome completeness as rated by checkm to keep a genome\n\t-b\tMaximum genome contamination as rated by checkm to keep a genome\n\t-m\tMaximum genome heterogeneity as rated by checkm to keep a genome\n\t-n\tDirectory containing NCBI taxonomy_dmp files (names.dmp and nodes.dmp)\n\t-q\tPath to Phyla Amphora marker database\n\t-h\tThis help screen\n\n\n"

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

		   fileonly="$(printf "$1" | awk -F"/" '{print $2}')"

		   species="$(printf "$genf" | awk -F"/" '{gsub( /['\'';:,.]/ , ""); print $2}' | awk -F"_" '{print $1"_"$2}')" 

		  	ncont="$(grep -c ">" "$genf")"
		   
		  	nfift="$(N50.py "$genf" | grep "N50" | awk '{print $2}')"

		   	header="$(head -1 "$genf" | awk '/^>/ {gsub(">", ""); print $1}')"

			afact="$(awk -v n50="$nfift" -v cont="$ncont" 'BEGIN { printf (n50 / cont)}')"


			mkdir -p "$outdir"/sort_temp/


			if [ "$nfift" -ge "$gn50" ] && [ "$ncont" -le "$nocont" ]; then

				printf "$header\t$ncont\t$nfift\t$afact\n" >> "$outdir"/sort_temp/"$species"_temp.txt

			  	awk -v header="$header" '/^>/ {i++; printf(">"header"_contig%06d\n", i); next}{print}' "$genf" > "$outdir"/genome_quality_analysis/genstat_passed/"$header".fna
				
				echo "$(basename $1) N50 = $nfift , #Contigs = $ncont , passed" >&2
		 
		 	else

				echo "$(basename $1) N50 = $nfift , #Contigs = $ncont , failed" >&2
			
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
## Prodigal function
################################################################################

get_prot_seqs() {

   genome_file="$( printf "$1" | awk -F"/" '{print $NF}')"
   genome_name="${genome_file%.fna}"

   mkdir -p "$outdir"/prodigal_output/prod_gff
   mkdir -p "$outdir"/prodigal_output/gene_prot_seqs
   mkdir -p "$outdir"/prodigal_output/gene_nucl_seqs

       prodigal -i "$1" -a "$outdir"/prodigal_output/gene_prot_seqs/"$genome_name".tmp -d "$outdir"/prodigal_output/gene_nucl_seqs/"$genome_name".tmp -o "$outdir"/prodigal_output/prod_gff/"$genome_name".gff -q

       awk '/^>/ {print $1; next}{print}' "$outdir"/prodigal_output/gene_prot_seqs/"$genome_name".tmp > "$outdir"/prodigal_output/gene_prot_seqs/"$genome_name".faa
       awk '/^>/ {print $1; next}{print}' "$outdir"/prodigal_output/gene_nucl_seqs/"$genome_name".tmp > "$outdir"/prodigal_output/gene_nucl_seqs/"$genome_name".ffn

       rm "$outdir"/prodigal_output/gene_prot_seqs/"$genome_name".tmp
       rm "$outdir"/prodigal_output/gene_nucl_seqs/"$genome_name".tmp

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
## Direct user script input to variables
################################################################################

while getopts hi:t:d:o:g:c:n:a:b:m:q: option
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
	  q) pard=${OPTARG};;
	  h | *)
		 usage
		 ;;
   esac
done


# Check if required variables have a given value and/or revert to defaults

if [ -z "$intax" ]; then

   	echo "Please provide a taxonomic name of genus or greater rank"

	usage

   	exit 1

fi


if [ -z "$outdir" ]; then

	outdir="$PWD/tax_gen_output"

	echo "no output directory specified, setting output directory to $outdir"

fi

if [ -z "$threads" ]; then

    threads="1"
	
fi

printf "\n\nUsing $threads threads"

	
if [ -z "$gn50" ]; then

	gn50="20000"

fi

printf "Minimum N50 set: $gn50\n"


if [ -z "$nocont" ]; then
	
	nocont="200"

fi

printf "Maximum No. Contigs set: $nocont\n"

if [ -z "$mincomp" ]; then

	mincomp="85"

fi

printf "Minimum genome clompleteness set: $mincomp\n"

if [ -z "$maxcont" ]; then

	maxcont="2"

fi

printf "Maximum genome contamination set: $maxcont\n"

if [ -z "$maxhetero" ]; then

	maxhetero=40

fi

printf "Maximum genome heterogeneity set: $maxhetero\n"

if [ -z "$pard" ]; then

	pard="/home/joseph/biotools/programs/Phyla_AMPHORA/Marker"

fi


# store running date/time

run_date=$(date)

date="$(date | sed 's/ /-/g')"

# export functions

export -f genome_stat
export -f split_genomes
export -f get_prot_seqs
export -f par_gzip

# export variables to work with functions

export outdir="$outdir"
export threads="$threads"
export nocont="$nocont"
export gn50="$gn50"


################################################################################
## Check that software dependencies are available
################################################################################

printf "\n\nChecking software dependencies are installed\n"

if [ -n "$(which checkm)" ] && [ -n "$(checkm -h)" ]; then

	printf "\ncheckm is installed: $(which checkm)\n"

else

	printf "\ncheckm is not installed (correctly) !!\n"
	missdep="yes"

fi

printf "\nCheck Phyla AMPHORA scripts\n"

if [ -n "$(which MarkerScanner.pl)" ]; then

	printf "MarkerScanner.pl is installed: $(which MarkerScanner.pl)\n"

else

	printf "MarkerScanner.pl is not installed (correctly) !!\n"
	missdep="yes"

fi

if [ -n "$(which MarkerAlignTrim.pl)" ]; then

	printf "MarkerAlignTrim.pl is installed: $(which MarkerAlignTrim.pl)\n"

else

	printf "MarkerAlignTrim.pl is not installed (correctly) !!\n"
	missdep="yes"

fi

if [ -n "$(which prodigal)" ]; then

	printf "\nProdigal is installed: $(which prodigal)\n"

else

	printf "\nProdigal is not installed (correctly) !!\n"
	missdep="yes"

fi

if [ -n "$(which parallel)" ]; then

	printf "\nParallel is installed: $(which parallel)\n"

else

	printf "\nParallel is not installed (correctly) !!\n"
	missdep="yes"

fi

if [ "$missdep" == "yes" ]; then

	printf "\nOne or more software dependencies cannot be called correctly by this workflow script, Exiting\n\n"

	exit 1

else

	printf "\nRequired dependencies have been found, starting workflow\n"
fi


################################################################################
## NCBI taxonomy files
################################################################################

# Check that the given directory is correct and the required NCBI taxonomy files are present otherwise prompt user if the script should retrieve tax_dump from the NCBI ftp site

if [ -z "$ncbitax" ]; then

	if [ -f "$outdir"/NCBI_tax/nodes.dmp ] && [ -f "$outdir"/NCBI_tax/names.dmp ]; then

		printf "NCBI files found in $outdir/NCBI_tax from previous run\n"

		ncbitax="$outdir/NCBI_tax"
		notax="0"

	else

        printf "\n\nDirectory path to ncbi tax dump files not supplied\n" 
        notax="1"

	fi

elif [[ ! $(ls -A "$ncbitax") ]]; then

        printf "\n\nDirectory path to ncbi tax dump files is empty\n"
        notax="1"

elif [ ! -f "$ncbitax"/nodes.dmp ]; then

        printf "\n\nnodes.dmp file not found at $ncbitax\n"
        notax="1"

elif [ ! -f "$ncbitax"/names.dmp ]; then

        printf "\n\nnames.dmp file not found at $ncbitax\n"
        notax="1"
fi

if [ "$notax" = "1" ]; then

        printf "Download NCBI tax_dump ?\n" 

        select get_taxdump in "Yes" "No"; do

                case $get_taxdump in
                        Yes )
                                mkdir -p "$outdir"/NCBI_tax
								
								if [ ! -f "$outdir"/NCBI_tax/new_taxdump.tar.gz ]; then
                                
									wget -O "$outdir"/NCBI_tax/new_taxdump.tar.gz ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz || exit 1
                                
									wget -O "$outdir"/NCBI_tax/new_taxdump.tar.gz.md5 ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5 || exit 1
                                
								fi

								while read md5code filename; do

									if [ "$(md5sum "$outdir"/NCBI_tax/"$filename" | awk '{print $1}')" = "$md5code" ]; then

                                                printf "md5sum check passed, extracting tax dump files"
                                                tar -xzvf "$outdir"/NCBI_tax/new_taxdump.tar.gz -C "$outdir"/NCBI_tax
                                                ncbitax="$outdir/NCBI_tax"
												printf "\nNCBITAX directory set to $ncbitax\n"

                                        else
											printf "\n problem with files\n"
                                                exit 1
											
									fi

                                done < "$outdir"/NCBI_tax/new_taxdump.tar.gz.md5
								
								if [ -f "$outdir"/NCBI_tax/nodes.dmp ] && [ -f "$outdir"/NCBI_tax/names.dmp ]; then

									printf "\ntax dump has been downloaded and extracted\n"

									break

								fi
								;;


                        No )
                                printf "\n\nYou have selected not to retrieve the NCBI taxonomy files, please provide a valid directory containing\nthe NCBI taxonomy nodes.dmp and names.dmp files when rerunning the script\n"

                                exit 1
								;;
						* )
								printf "\nPlease select a proper option 1 (yes) or 2 (no)\n"
								;;
                esac

        done

fi
 	
echo "point A cleared" 

# ncbitax="/biodata/database/NCBI-taxonomy/taxdump"

 export nodefile="$ncbitax/nodes.dmp"
 export namefile="$ncbitax/names.dmp"
 #export gbtax="$ncbitax/nucl_gb.accession2taxid"
 #export wgstax="$ncbitax/nucl_wgs.accession2taxid"
 
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
	
	   	printf "\n\n\nTaxonomy ID $intax provided, known as $taxname\n"
	
	else
	 
		parenttax="$( grep "$intax" "$namefile" | awk -v tax="$intax" -F"\t[|]\t" ' $2 == tax {print $1}')"
	
		taxname="$intax"
	
	 	intaxrank="$(grep "$parenttax" "$nodefile" | awk -F "\t[|]\t" -v tax="$parenttax" '$1 == tax {print $3}')"
	 	 
	 	printf "\n\n\nTaxonomic name $taxname provided, Taxon ID = $parenttax\n"
	 
	fi
	 
	printf "\nparenttax = $parenttax\n"
	 
	## Determine the phylum of the assigned taxonomy, this will be important later for phyla amphora
	 
	uptax="$parenttax"
	 
	printf "\n\nFinding phylum $taxname belongs to\n\n"
	 
	until [ "$taxrank" = "class" ]; do
	 
		taxrank="$(grep "$uptax" "$nodefile" | awk -F"\t[|]\t" -v tax="$uptax" '$1 == tax {print $3}')"
	 	
		#printf "taxrank = $taxrank\n"
		#printf "uptax = $uptax\n"
		#printf "uptaxname = $uptaxname\n"

	 	if [ "$taxrank" = "order" ]; then
	 
			genorder="$uptaxname"
	 
		 	printf "Genus $taxname belongs to the Order $genorder\n"
		 
		fi
		 
		uptax="$(grep "$uptax" "$nodefile" | awk -F"\t[|]\t" -v tax="$uptax" '$1 == tax {print $2}')"
		 
		uptaxname="$(grep "$uptax" "$namefile" | grep "scientific name" | awk -F"\t[|]\t" -v tax="$uptax" ' $1 == tax {print $2}')"
		    
		if [[ $uptaxname = *"proteobacteria"* ]]; then
		 
		 	printf "\n$taxname is part of Proteobacteria, breaking at class specific $uptaxname \n\n"
		 
		 	break
		 
		fi
		 
		let "iteration++"
		
		printf " $iteration..\n"
		
	done

else

	printf "\nRequired NCBI taxonomy files not found, exiting\n\n"

	exit 1

fi
 
		# print vars for error checking
		#printf "taxrank = $taxrank\n"
		#printf "uptax = $uptax\n"

phylum="$(grep "$uptax" "$namefile"| grep "scientific name" |  awk -F"\t[|]\t" -v tax="$uptax" '$1== tax {print $2}')"
 
phylum="$(printf "$phylum" | awk '{print $1}')"
 
printf "\n\n$taxname identified as belonging to the phylum $phylum\n\n\n"

if [ -z "$phylum" ]; then

	printf "The provided taxonomy or taxid may be invalid, exiting script\n\n\n"

	exit 1

fi

#echo "point B cleared"


################################################################################
## Set correct Phyla-amphora marker gene set based on the identified phylum
################################################################################

phylum_list=(phylum Alphaproteobacteria Betaproteobacteria Gammaproteobacteria Deltaproteobacteria Epsilonproteobacteria Acidobacteria Actinobacteria Aquificae Bacteroidetes Chlamydiae Chlorobi Chloroflexi Cyanobacteria Deinococcus-Thermus Firmicutes Fusobacteria Planctomycetes Spirochaetes Tenericutes Thermotogae Verrucomicrobia)

for i in "${!phylum_list[@]}"; do

	if [[ "${phylum_list[$i]}" = "$phylum" ]]; then

		pamsp="$i"
		
		if [[ "$pamsp" = "21" ]]; then	# under phyla amphora Chlamydiae and Verrucomicrobia use the same marker set

			pamsp="10"

		fi

	fi

done

if [ -z "$pamsp" ]; then

	echo "Phyla-amphora phylum could not be set correctly, exiting"

	exit

fi



printf "\npamsp set to $pamsp\n"	# print pamsp value to screen for error checking


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

   		printf "input taxonomy is a genus\n\n"

	else

		printf "\nFind all genera under $taxname\n\n"

		grep "$parenttax" "$nodefile" | awk -F"\t[|]\t" -v ptax="$parenttax" ' $2 == ptax {print $1"\t"$3}' > "$outdir"/tax_list.txt	# Get a list of the child taxon IDs under the parent taxon id


	## Until loop required if a higher than 'Family' rank taxon was specified. Will continue until the taxon IDs for child 'genera' is produced

	let "iteration = 0"

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

		let "iteration=iteration+1"

		printf "$iteration.."
	
	done

	nogenera="$(wc -l $outdir/genus_list.txt | awk '{print $1}')"

	printf "%s" "$taxname" >> "$outdir"/genus_list.txt


	printf "\n$nogenera genera found\n\n\n"

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

while read genus; do

	printf "  Retrieving genomes from NCBI for $genus genus\t"

	rsync --quiet --no-motd --msgs2stderr --exclude="*_cds_from_genomic.fna.gz" --exclude="*_rna_from_genomic.fna.gz"  ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/${genus}_*/latest_assembly_versions/*/*_genomic.fna.gz "$outdir"/genome_quality_analysis/retrieved_genomes/ 2>/dev/null || printf "No genomes found!\n"

	printf"Total genomes = %d\n" "$(ls "$outdir"/genome_quality_analysis/retrieved_genomes/ | wc -l)"

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

		printf "\n  Building "$taxname"_name_key.txt file\n"

		for file in "$outdir"/genome_quality_analysis/retrieved_genomes/*.fna ; do

			head -1 "$file" | awk -F"," '{print $1}' | sed -e 's/>//g' -e 's/ $//g' -e 's/^ //g' -e 's/:/-/g' -e 's/\S*scaffold\S*//g' -e 's/\S*contig\S*//g' -e 's/\(\[\|\]\|{\|}\||\)//g' -e "s/\((\|)\|'\)//g" | awk '{acc = $1; $1 = ""; print acc"\t"$0}' >> "$outdir"/genome_information/"$taxname"_name_key.txt

		done

		## New way to get more NCBI genome information

		# if the 'prokaryotes.txt' file is not present, retrieve it from NCBI ftp

		if [ ! -f "$outdir"/prokaryotes.txt ]; then

			printf "\n  Retrieving the prokaryotes.txt file from NCBI ftp server\n\n"

			wget -O "$outdir"/prokaryotes.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
	
		fi

		declare -A genome_acc_gca
	
		printf "\n  Linking Accession and Genbank accessions for genomes\n"

		for file in "$outdir"/genome_quality_analysis/retrieved_genomes/*.fna ; do

			gca="$(printf "${file##*/}" | awk -F"_" '{print $1"_"$2}')"

			genome_acc="$(head -1 "$file" | awk '{gsub(">",""); print $1}')"

			genome_acc_gca[$genome_acc]=$gca
		
		done

		if [ ! -f "$outdir"/"$taxname"_acc_gca.txt ]; then

			for i in "${!genome_acc_gca[@]}"; do

				printf "$i\t${genome_acc_gca[$i]}\n" >> "$outdir"/"$taxname"_acc_gca.txt

			done

		fi

		# loop through the genome GCA accessions to gather useful information from prokaryote.txt and save to a tab delimited file

		if [ ! -f "$outdir"/genome_information/"$taxname"_genome_char.tab ] || [ ! -f "$outdir"/genome_information/"$taxname"_genome_info.tab ]; then

			printf "\n  Building genome information and general characteristics files\n"

			printf "Accession\tName\tGroup\tSubgroup\ttaxid\tSize\tGC%%\tContigs\tGenes\tProteins\n" > "$outdir"/genome_information/"$taxname"_genome_char.tab
		
			printf "Accession\tName\tBioproject Accession\tBiosample Accession\tCentre\tRelease Date\tModify Date\tStatus\tPubmed ID\n" > "$outdir"/genome_information/"$taxname"_genome_info.tab

			for i in "${!genome_acc_gca[@]}"; do

				awk -F"\t" -v gca="${genome_acc_gca[$i]}" -v acc="$i" '$19 ~ gca {print acc"\t"$1"\t"$5"\t"$6"\t"$2"\t"$7"\t"$8"\t"$11"\t"$12"\t"$13}' "$outdir"/prokaryotes.txt >> "$outdir"/genome_information/"$taxname"_genome_char.tab

				awk -F"\t" -v gca="${genome_acc_gca[$i]}" -v acc="$i" '$19 ~ gca {print acc"\t"$1"\t"$3"\t"$18"\t"$17"\t"$14"\t"$15"\t"$16"\t"$22}' "$outdir"/prokaryotes.txt >> "$outdir"/genome_information/"$taxname"_genome_info.tab

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

		qual_genome_list+=($(awk -v gen_comp="$min_comp" -v gen_cont="$maxcont" -v gen_hetero="$maxhetero" -F"\t" 'NR>1&&$2!~/k__Bacteria/&&$(NF-2)>=gen_comp&&$(NF-1)<gen_cont&&$NF<gen_hetero{print $1}' "$outdir"/genome_quality_analysis/checkm_lwf/"$grp_dir"/lineage_wf.tab))


	done

    mkdir -p "$outdir"/quality_genomes

	# Extract genome information on only remaining quality filtered genomes
	printf "Accession\tName\tBioproject Accession\tBiosample Accession\tCentre\tRelease Date\tModify Date\tStatus\tPubmed ID\n" > "$outdir"/genome_information/"$taxname"_quality_genome_info.tab

	printf "Accession\tName\tGroup\tSubgroup\ttaxid\tSize\tGC%%\tContigs\tGenes\tProteins\n" > "$outdir"/genome_information/"$taxname"_quality_genome_char.tab

	for i in ${!qual_genome_list[@]}; do

		cp "$outdir"/genome_quality_analysis/genstat_passed/"${qual_genome_list[$i]}".fna "$outdir"/quality_genomes/

	grep "^"${qual_genome_list[$i]}"" "$outdir"/genome_information/"$taxname"_genome_info.tab >> "$outdir"/genome_information/"$taxname"_quality_genome_info.tab

	grep "^"${qual_genome_list[$i]}"" "$outdir"/genome_information/"$taxname"_genome_char.tab >> "$outdir"/genome_information/"$taxname"_quality_genome_char.tab

	done


fi




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

	printf "\nAdding genomes from "$indir" to quality selected genomes\n"

fi

# Count number of remaining genomes

gencount="$(ls "$outdir"/quality_genomes/ | wc -l | awk '{print $1}')"

printf "\n$gencount high quality genomes selected\n"


## Add post-checkm filtering information


printf "  Checkm genome filtering\n  --------------------------------------------------------------------------------\n  The genomes that passed assembly quality filtering were analysed using checkm to\n  determine the genome completeness, contamination and heterogeneity. A minimum\n  genome completeness of $mincomp %%, maximum contamination of $maxcont %%, and maximum heterogeneity\n  of $maxhetero %% was required to pass filtering. The genomes that passed this fitering were\n  stored in the \""$outdir"/quality_genomes/\" directory. Following filtering\n  genomes supplied via input \""$indir"\" were added to this directory.\n\n  \tCheckm passed genomes\t\t=\t$no_quality_genomes\n\n  \tInput genomes added\t\t=\t$gencount\n\n\n" >> "$outdir"/get_genome_log_"$date".txt

tail -n +2 "$outdir"/genome_information/"$taxname"_quality_genome_info.tab | awk '{print $2}' | sort |uniq -c | awk 'BEGIN{printf "       Quality Genomes per Genus\n     -------------------------------\n"}{printf "     %20s  =  %d\n", ($2), ($1)}END{printf "\n\n"}' >> "$outdir"/get_genome_log_"$date".txt


################################################################################
## Get full taxonomic nomenclature for quality genomes using NCBI taxonomy
################################################################################

# This method is old and slow because it has to query very large accession-to-taxid files that are filed with every NCBI entry and takes grep several seconds to check each one for a matching taxid this has been replaced above with a new method using the prokaryote.txt file which also gives access to more genome metadata and characteristics

# if [ ! -f "$outdir"/"$taxname"_tax_nomenclature.txt ]; then
# 
# 	for file in "$outdir"/quality_genomes/*.fna; do
# 
# 		gfile="${file##*/}"
# 
# 		genomeacc="$( printf "${gfile/.fna}" | awk -F"." '{print $1}')"
# 
# 		printf "\ngenomeacc=$genomeacc, "
# 
# 		genometax="$( grep "\<$genomeacc\>" "$gbtax" | awk -F"\t" -v gacc="$genomeacc" '$1 == gacc {print $3}')"
# 
# 		
# 		if [ -z "$genometax" ]; then
# 
# 			genometax="$( grep "\<$genomeacc\>" "$wgstax" | awk -F"\t" -v gacc="$genomeacc" '$1 == gacc {print $3}')"
# 
# 			#printf "$genometax\t"
# 
# 			if [ -z "$genometax" ]; then
# 
# 				printf "\nNo taxonomic name found for $file\n"
# 
# 				continue
# 
# 			fi
# 
# 		fi
# 
# 		printf "genometax=$genometax, "
# 

rm "$outdir"/genome_information/"$taxname"_tax_nomenclature.txt 2>/dev/null

if [ ! -f "$outdir"/genome_information/"$taxname"_tax_nomenclature.txt ]; then

	printf "\nBuilding tax_nomenclature file\n"

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
## Predict protein coding ORFs in genomes using prodigal
################################################################################

if [[ ! $(ls -A "$outdir"/prodigal_output/gene_prot_seqs/ 2>/dev/null) ]]; then

	printf "\nPredicting ORFs in genomes\n"

	find "$outdir"/quality_genomes/ | grep -e \\$extension | parallel -P "$threads" --no-notice get_prot_seqs

fi

## Add prodigal ORF info to log

printf "  Predicting open reading frames (ORF)\n  --------------------------------------------------------------------------------\n  Open reading frames were found for each of the remaining high quality genomes\n  using prodigal ORF prediction  software. ORF prediction statistics are listed\n  below:\n\n" >> "$outdir"/get_genome_log_"$date".txt

grep -c ">" "$outdir"/prodigal_output/gene_prot_seqs/*.faa | awk -F":" '{sum = sum + $2; if ($2 > max) max = $2; if(min>$2||NR==1)min = $2}END{print "  \tNo. genomes\t\t\t=\t"NR"\n  \tTotal ORFs identified\t\t=\t"sum"\n  \tMean ORFs per genome\t\t=\t"sum/NR"\n  \tMin ORF\t\t\t\t=\t"min"\n  \tMax ORF\t\t\t\t=\t"max"\n  \tRange\t\t\t\t=\t"max-min"\n\n"}' >> "$outdir"/get_genome_log_"$date".txt

printf "  The ORFs were merged into one file for Phyla AMPHORA marker scanning\n\n\n" >> "$outdir"/get_genome_log_"$date".txt


## Concatenate the protein ORFs into a single file

mkdir -p "$outdir"/phyla_amphora_analysis

if [ ! -f "$outdir"/phyla_amphora_analysis/"$taxname"_combined.faa ]; then

	printf "\nCombining protein ORFs\n\n"

	cat "$outdir"/prodigal_output/gene_prot_seqs/*.faa > "$outdir"/phyla_amphora_analysis/"$taxname"_combined.faa

fi


################################################################################
## Run Marker scanner to find marker genes
################################################################################

mkdir -p "$outdir"/phyla_amphora_analysis/marker_scan

if [[ ! $(ls -A "$outdir"/phyla_amphora_analysis/marker_scan 2>/dev/null ) ]]; then

	printf "\nFinding marker genes"

	MarkerScanner.pl -Phylum "$pamsp" -ReferenceDirectory "$pard" "$outdir"/phyla_amphora_analysis/"$taxname"_combined.faa

	mv *.pep "$outdir"/phyla_amphora_analysis/marker_scan/

fi

## Select the marker genes with 1 sequences from each genome


if [[ ! $(ls -A "$outdir"/phyla_amphora_analysis/good_markers/ 2>/dev/null ) ]]; then

	mkdir -p "$outdir"/phyla_amphora_analysis/good_markers

	printf "\nSelect good marker genes"

	for file in "$outdir"/phyla_amphora_analysis/marker_scan/*.pep; do

		if [[ $(awk -F"_" '/^>/ {print $(NF-2)}' $file | sort | uniq -c | awk '$1 == "1" {count++} END {print count}') == $gencount ]]; then
			cp "$file" "$outdir"/phyla_amphora_analysis/good_markers/
		fi

	done

	mgenecount="$(ls "$outdir"/phyla_amphora_analysis/good_markers/ | wc -l | awk '{print $1}')"

	if [ "$mgenecount" -lt "5" ]; then

		printf "Less than 5 compliant marker genes have been found for this set of genomes\nWhich will result in poor phylogenetic analysis, exiting script\n"

		exit 1

	else

		printf "\n$mgenecount compliant marker genes found"

	fi

	if [[ ! $( find *.pep 2>/dev/null) ]]; then

		cp "$outdir"/phyla_amphora_analysis/good_markers/*.pep $PWD/

	fi

fi

## Add Phyla AMPHORA marker scan to log

printf "  Phyla AMPHORA analysis\n  --------------------------------------------------------------------------------\n  The predicted ORFs from all the genomes combined were scanned using Phyla\n  AMPHORA Markerscanner.pl against the $phylum marker gene set. The sets of\n  marker genes generated were then checked for conformity for phylogenetic analysis.\n  Conformity required that for each marker gene, a single gene sequence from each\n  remaining genome analysed was present in the set. If a marker gene set did not\n  follow this conformity it was removed from the remaining analysis.\n\n" >> "$outdir"/get_genome_log_"$date".txt


grep ">" -c "$outdir"/phyla_amphora_analysis/marker_scan/*.pep | awk -F":" '{sum = sum + $2; if($2 > max) max = $2; if($2<min||NR==1)min = $2}END{print "\n  \tAll Marker Statistics:\n\n  \tNo. marker genes\t=\t"NR"\n  \tTotal genes\t\t=\t"sum"\n  \tMean\t\t\t=\t"sum/NR"\n  \tMin\t\t\t=\t"min"\n  \tMax\t\t\t=\t"max"\n  \tRange\t\t\t=\t"max-min}' >> "$outdir"/get_genome_log_"$date".txt

printf "\n\n  \tGood Marker Statistics:\n\n" >> "$outdir"/get_genome_log_"$date".txt

grep ">" -c "$outdir"/phyla_amphora_analysis/good_markers/*.pep | awk -F":" '{sum = sum + $2; if($2 > max) max = $2; if($2<min||NR==1)min = $2}END{print "  \tNo. marker genes\t=\t"NR"\n  \tTotal genes\t\t=\t"sum"\n  \tMean\t\t\t=\t"sum/NR"\n  \tMin\t\t\t=\t"min"\n  \tMax\t\t\t=\t"max"\n  \tRange\t\t\t=\t"max-min}' >> "$outdir"/get_genome_log_"$date".txt


################################################################################
## Run MarkerAlignTrim.pl script to align the good marker genes
################################################################################

if [[ ! $(ls -A "$outdir"/phyla_amphora_analysis/marker_alignments/ 2>/dev/null) ]]; then

	if [[ ! $( find *.aln 2>/dev/null) ]] && [[ $(find *.pep 2>/dev/null) ]]; then

		printf "\nAligning compliant marker genes"

		MarkerAlignTrim.pl -Trim -ReferenceDirectory "$pard" -OutputFormat fasta

	fi

	# output alignments are cleaned up and moved to marker_alignments directory

	mkdir -p "$outdir"/phyla_amphora_analysis/marker_alignments

	for file in *.aln; do

		awk -F"_" '/^>/ {print $1"_"$(NF-2); next}{print}' "$file" > "$outdir"/phyla_amphora_analysis/marker_alignments/"$file"

	done

fi

## check that alignment files have been correctly moved then clean up the working directory

if [[ $( ls -A "$outdir"/phyla_amphora_analysis/marker_alignments/ 2>/dev/null) ]]; then


	if [[  $(find *.pep 2>/dev/null) ]]; then

		rm "$PWD"/*.pep

	fi

	if [[ $(find *.aln 2>/dev/null) ]]; then

		rm "$PWD"/*.aln

	fi

	if [[ $(find *.mask 2>/dev/null) ]]; then

		rm "$PWD"/*.mask

	fi

else

	printf "\nSomething has gone wrong. Exiting"

	exit

fi


## use the script catfasta2phyml.pl to merge the alignment files into one

if [ ! -f "$outdir"/phyla_amphora_analysis/"$taxname"_comb_aln.fasta ]; then

	printf "\nConcatentating the aligned marker genes together"

	catfasta2phyml.pl -f "$outdir"/phyla_amphora_analysis/marker_alignments/*.aln > "$outdir"/phyla_amphora_analysis/"$taxname"_comb_aln.fasta

fi

## Add Marker alignment info to log file

printf "\n\n  The Good marker genes were aligned using MarkerAlignTrim.pl and the resulting\n  alignments were concatenated using catfasta2phyml.pl. The resulting concatenated\n  alignment was then ready for phylogenetic analysis.\n\n\n" >> "$outdir"/get_genome_log_"$date".txt

awk '/^>/{i++}!/>/{if(i>1) exit; print}' "$outdir"/phyla_amphora_analysis/"$taxname"_comb_aln.fasta | awk '{seqlen += length($0)}END{print "  \tAlignment length\t=\t"seqlen" AA\n\n\n"}' >> "$outdir"/get_genome_log_"$date".txt

awk 'BEGIN{RS=">";ORS=""}{$1 =""; gsub(" ",""); print $0}' "$outdir"/phyla_amphora_analysis/"$taxname"_comb_aln.fasta | sed -e 's/-//g' -e 's/./&\n/g' | sort | uniq -c | awk 'BEGIN{print "  \tAlignment composition\n  \t---------------------\n"}{AA[$i] = $2; count[$i] = $1}{sum += $1}END{for (i in AA) {printf "           %s  =  %.3f %\n", (AA[i]), (count[i]/sum*100)}; print "\n\n\n"}' >> "$outdir"/get_genome_log_"$date".txt


################################################################################
## run RAxML on the concatenated alignment
################################################################################

if [[ ! $( ls -A "$outdir"/phyla_amphora_analysis/raxml_analysis 2>/dev/null) ]]; then

	printf "\nPerforming RAXML maximum likelyhood phylogenetic analysis"

	mkdir -p "$outdir"/phyla_amphora_analysis/raxml_analysis

	raxmlHPC-PTHREADS -f a -s "$outdir"/phyla_amphora_analysis/"$taxname"_comb_aln.fasta -x 127 -# 100 -m PROTGAMMAJTT -n "$taxname"_100b -p 129 -T "$threads"

	mv RAxML* "$outdir"/phyla_amphora_analysis/raxml_analysis/

fi

## Rename the nodes of the phylogenetic trees produced by RAxML using the naming lists

if [[ ! $(ls -A "$outdir"/phyla_amphora_analysis/raxml_analysis/renamed_tree 2>/dev/null) ]] && [[ $(ls -A "$outdir"/phyla_amphora_analysis/raxml_analysis 2>/dev/null) ]]; then


		printf "\nRenaming nodes of ML trees\n"

		mkdir -p "$outdir"/phyla_amphora_analysis/raxml_analysis/renamed_tree

		export taxname="$taxname"

		rename_nodes.sh -i "$PWD"/"$outdir"/phyla_amphora_analysis/raxml_analysis/RAxML_bestTree."$taxname"_100b -n "$PWD"/"$outdir"/genome_information/"$taxname"_name_key.txt -o "$PWD"/"$outdir"/phyla_amphora_analysis/raxml_analysis/renamed_tree/"$taxname"_MLtree_100bs_fnamed_no_bs_labels.nwk

		rename_nodes.sh -i "$outdir"/phyla_amphora_analysis/raxml_analysis/RAxML_bipartitions."$taxname"_100b -n "$outdir"/genome_information/"$taxname"_name_key.txt -o "$outdir"/phyla_amphora_analysis/raxml_analysis/renamed_tree/"$taxname"_MLtree_100bs_fnamed_with_bs_labels.nwk

		rename_nodes.sh -i "$PWD"/"$outdir"/phyla_amphora_analysis/raxml_analysis/RAxML_bestTree."$taxname"_100b -n "$PWD"/"$outdir"/genome_information/"$taxname"_tax_naming.txt -o "$PWD"/"$outdir"/phyla_amphora_analysis/raxml_analysis/renamed_tree/"$taxname"_MLtree_100bs_tnamed_no_bs_labels.nwk
		
		rename_nodes.sh -i "$outdir"/phyla_amphora_analysis/raxml_analysis/RAxML_bipartitions."$taxname"_100b -n "$outdir"/genome_information/"$taxname"_tax_naming.txt -o "$outdir"/phyla_amphora_analysis/raxml_analysis/renamed_tree/"$taxname"_MLtree_100bs_tnamed_with_bs_labels.nwk


fi

## Add RAXML analysis info to log file

printf "  RAXML maximum likelihood phylogenetic analysis\n  --------------------------------------------------------------------------------\n  RAXML was used to perform a phylogenetic analysis on the concatenated alignment.\n  The parameters used in the phylogenetic analysis are listed below:\n\n  \tAlignment\t\t=\tAmino acid\n  \tMethod\t\t\t=\tMaximum Likelihood\n  \tEvolutionary model\t=\tGamma distritubion\n  \tSubstitution matrix\t=\tJones, Taylor & Thornton (JTT)\n  \tTest of phylogeny\t=\tBootstrap\n  \tReplicates\t\t=\t100\n\n\n  Following the RAXML analysis the tips of the phylogenetic tree produced were\n  renamed using either "$taxname"_name_key.txt or "$taxname"_tax_naming.txt files.\n\n" >> "$outdir"/get_genome_log_"$date".txt 

printf "  Renamed trees stored in "$outdir"/phyla_amphora_analysis/raxml_analysis/renamed_tree\n\n  ML tree no bootstrap renamed "$taxname"_name_key.txt: "$taxname"_MLtree_100bs_fnamed_no_bs_labels.nwk\n  ML tree with bootstrap renamed "$taxname"_name_key.txt: "$taxname"_MLtree_100bs_fnamed_with_bs_labels.nwk\n  ML tree no bootstrap renamed "$taxname"_tax_naming.txt: "$taxname"_MLtree_100bs_tnamed_no_bs_labels.nwk\n  ML tree with bootstrap renamed "$taxname"_tax_naming.txt: "$taxname"_MLtree_100bs_tnamed_with_bs_labels.nwk\n\n" >> "$outdir"/get_genome_log_"$date".txt


cat "$outdir"/get_genome_log_"$date".txt >> "$outdir"/"$taxname"_analysis_report.txt


################################################################################
## Final cleanup of output directory
################################################################################

if [[ $( ls -A "$outdir"/phyla_amphora_analysis/raxml_analysis/renamed_tree 2> /dev/null) ]]; then

if [ -f "$outdir"/"$taxname"_combined.faa ]; then

	mv "$outdir"/"$taxname"_combined.faa "$outdir"/phyla_amphora_analysis/

fi

if [ -f "$outdir"/"$taxanme"_comb_aln.fasta ]; then

	mv "$outdir"/"$taxname"_comb_aln.fasta "$outdir"/phyla_amphora_analysis/
	
fi

if [ -f "$outdir"/"$taxname"_comb_aln.fasta.reduced ]; then

	rm "$outdir"/"$taxname"_comb_aln.fasta.reduced

fi

	rm -r "$outdir"/tmp_gen_list 2>/dev/null

else

	printf "\nRenamed trees not produced\n"

	exit 1

fi



exit


		 







