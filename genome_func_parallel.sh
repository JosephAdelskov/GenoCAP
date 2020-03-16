#! /bin/bash


#: Title            : Genome func parallel
#: Date             : 10/07/2018
#: Author           : Joseph Adelskov
#: Version          : 1.0
#: Description		: Continues on from refine_cluster_genes.sh script. Extract functional UniProt information from PROKKA annotations. links genes with COG, KO and GO terms via UniProt IDs.

## Testing parallelizaiton versions of COMP_GENOME_FUNC_ANALYSIS
usage() { printf "\nCOMP_GENOME_FUNC_ANALYSIS.sh : Continues from REFINE_CLUSTER_GENES.sh by performing\nfunctional annotation analysis of core-pan gene clusters from the refined set of genomes\n\n\t-i\tPath to output directory from REFINE_CLUSTER_GENES.sh script\n\t-t\tNumber of cpu threads to use\n\t-k\tPath to KEGG orthology map file\n\n"

	exit;
}


################################################################################
## Get Uniprot flatfile function
################################################################################

get_uniprot_ff() {

					if [ ! -f "$outdir"/uniprot_info/"$1".txt ]; then

						wget -q -O "$outdir"/uniprot_info/"$1".txt www.uniprot.org/uniprot/"$1".txt

						printf "$1.txt retrieved\n" >&2

					fi

}

################################################################################
## Get gene function IDs function
################################################################################

get_func_ids() {

	filename="${1%.*}"

   if [ ! -f "$outdir"/genome_func_ids/"$filename"_func_ids.tab ]; then

#	printf "\nExtracting CDS function identifiers for $filename\n"      
#	printf "$gene_count/$gene_no"

	if [ -f "$outdir"/prokka_annotation/"$filename".gff ]; then

		printf "Extracting function IDs for %s\n" "$filename" >&2

		printf "Gene_id\tUniprot_id\tKO_id\tCOG_id\tGO_id\n" > "$outdir"/genome_func_ids/"$filename"_func_ids.tab

		IFS=$'\t'

		while read cont_id tool type start stop dot direction zero annotation; do

			if [ "$type" = "CDS" ]; then

		   	#let "gene_count++"

		   	#printf "\rParsing $filename, Total = $gene_count/$gene_no" >&2


				gene_id="$(printf "%s" "$annotation" | awk -F";" '{gsub("ID=", ""); print $1}')"

				uprot_id="$(printf "%s" "$annotation" | awk -F":" '{for(i=1;i<=NF;i++) if($i == "UniProtKB") print $(i+1)}' | awk -F";" '{print $1}')"


				if [ -z "$uprot_id" ]; then

					hypoprot="TRUE"

					uprot_id="Hypothetical protein"

					printf "$gene_id\tHypothetical protein\n" >> "$outdir"/genome_func_ids/"$filename"_func_ids.tab
				
				 else

					if [ ! -f "$outdir"/uniprot_info/"$uprot_id".txt ]; then

						wget -q -O "$outdir"/uniprot_info/"$uprot_id".txt www.uniprot.org/uniprot/"$uprot_id".txt

					fi

					unset KO_id

					while read ko_line; do

						KO_id+=($ko_line)

					done < <(grep "^DR" "$outdir"/uniprot_info/"$uprot_id".txt | awk '$2 ~/KO;/{gsub(";",""); print $3}')
				
					unset KO_list

					if [ -n "$KO_id" ]; then

						KO_list="$( printf "%s," ${KO_id[@]})"

					fi

					unset COG_id

					while read cog_line; do

						COG_id+=($cog_line)

					done < <(grep "^DR" "$outdir"/uniprot_info/"$uprot_id".txt | awk -F";" '$2 ~ /COG */ {gsub(" ",""); print $2}')

					unset COG_list

					if [ -n "$COG_id" ]; then
					 
						COG_list="$( printf "%s," ${COG_id[@]})"

		 			fi

					unset GO_id

					# use loop to extract all GO ids from uniprot flatfile to array

					while read go_line; do

						GO_id+=($go_line)

					done < <(grep "^DR" "$outdir"/uniprot_info/"$uprot_id".txt | awk -F";" '$2 ~ /GO:[0-9]{7}/ {gsub(" ",""); print $2}')

					unset GO_list

					if [ -n "$GO_id" ]; then

				   		GO_list="$( printf "%s," ${GO_id[@]})"

					fi

					printf "%s\t%s\t%s\t%s\t%s\n" "$gene_id" "$uprot_id" "${KO_list%%,}" "${COG_list%%,}" "${GO_list%%,}" >> "$outdir"/genome_func_ids/"$filename"_func_ids.tab
					
				fi	
					
			fi

		done < <( grep "^gnl|" "$outdir"/prokka_annotation/"$filename".gff)

	fi

fi

}


################################################################################
## process script input options
################################################################################

while getopts hi:t:k: option 
do	
	case "${option}"
		in
		i) indir=${OPTARG};;
		t) threads=${OPTARG};;
	 	k) kegg_file=${OPTARG};;
		h | *)
			usage ;;
	esac
done

if [ -z "$indir" ]; then

	usage

fi

if [ -z "$threads" ]; then

	threads="1"

fi

## Prokka analysis is performed on the selected genomes in refine_cluster_genes.sh
## Prepare an array of Prokka .gff files

unset infile

for file in "$indir"/prokka_annotation/*.gff; do

	infile+=($(printf "${file##*/}"))
	
done

# prepare array of genome id only

for i in "${!infile[@]}"; do

	filename[$i]="${infile[$i]%.*}"

done

mkdir -p "$indir"/genome_func_ids
mkdir -p "$indir"/uniprot_info

# set outdir as indir
outdir="$indir"

# set taxname

taxname="${indir##*/}"

# Export functions and variables

export outdir
export -f get_func_ids
export -f get_uniprot_ff

## Get date and time

date="$(date | sed 's/ /-/g')"


################################################################################
## Check that software dependencies are available
################################################################################

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
## Extract protein functional annotation information from .gff file
################################################################################


## alternative way to get Uniprot ID flat files faster ( in parallel )

# get a non-redundant list of uprot ids

while read uprot_id; do

	if [ ! -f "$outdir"/uniprot_info/"$uprot_id".txt ]; then

		uprot_ids+=($uprot_id)

	fi

done < <(awk -F"\t" '{print $9}' "$indir"/prokka_annotation/*.gff | awk -F":" '{for(i=1;i<=NF;i++) if($i == "UniProtKB") print $(i+1)}' | awk -F";" '{print $1}' | sort | uniq)


printf "\n${#uprot_ids[@]} missing UniProt flat files.\n\n"


# get all required uniprot flat files in parallel

if [ -n  "$uprot_ids" ]; then


printf "Retrieving missing UniProt flatfiles for IDs from www.uniprot.org\n"

parallel -P "$threads" get_uniprot_ff ::: "${uprot_ids[@]}"

else

	if [[ ! "$(ls -A "$outdir"/uniprot_info/ 2>/dev/null)" ]]; then
		
		printf "\nNo UniProt IDs found, something may be wrong\n"
		exit

	fi

fi

# Build func_id files from genome .gff files in parallel
	
gene_no="$( awk -F"\t" '$3 ~ "CDS" {print}' "$outdir"/prokka_annotation/*.gff | wc -l)"
	
let "gene_count = 0"

export gene_no
export gene_count

parallel -P "$threads" get_func_ids ::: "${infile[@]}"

## Copy genome func ID files to new directory

if [[ ! "$(ls -A "$outdir"/func_id_files/ 2>/dev/null)" ]]; then

	mkdir -p "$outdir"/func_id_files

	cp "$outdir"/genome_func_ids/*.tab "$outdir"/func_id_files/

fi

## Start adding to analysis_report.txt

printf "\n  ================================================================================\n\n  Running script genome_func_parallel.sh\t$date\n\n\n  Script running parameters\n  --------------------------------------------------------------------------------\n  \tInput directory\t\t=\t"$indir"\n  \tKEGG file\t\t=\t"$kegg_file"\n  \tNumber of threads\t=\t"$threads"\n\n" >> "$outdir"/func_log_"$date".txt

## Add genome func ID extraction to analysis report

printf "\n  Building functional ID profiles from PROKKA annotations\n  --------------------------------------------------------------------------------\n  Protein coding sequence (CDS) annotation information was retrieved from the\n  PROKKA annotations produced by the previous script. The information includes\n  linked Uniprot Accession IDs for which flatfiles were retrieved from:\n  www.uniprot.org. The Uniprot flatfiles were searched for linked KEGG orthology\n  (KO), Gene Ontology (GO), and Clusters of Orthologous sequences (COG) IDs. For\n  each genome analysed a tab delimited file containing linked IDs for each CDS was\n  produced. These func_id.tab files have been stored in the directory:\n\n  "$outdir"/func_id_files\n\n\n" >> "$outdir"/func_log_"$date".txt

printf "  \tProkka annotation statistics\n  \t--------------------------------------------------\n  \tGenomes\t\t\t\t=%'10d\n" "$(ls "$outdir"/prokka_annotation/*.faa | wc -l)" >> "$outdir"/func_log_"$date".txt

printf "  \tTotal CDS\t\t\t=%'10d\n" "$(tail -n +2 "$outdir"/genome_func_ids/*func_ids.tab | wc -l)" >> "$outdir"/func_log_"$date".txt

printf "  \tHypothetical CDS\t\t=%'10d\n" "$(tail -n +2 "$outdir"/genome_func_ids/*func_ids.tab | awk -F"\t" '$2 == "Hypothetical protein"{count++}END{print count}')" >> "$outdir"/func_log_"$date".txt

printf "  \tCDS with Uniprot ID\t\t=%'10d\n" "$(tail -n +2 "$outdir"/genome_func_ids/*func_ids.tab | awk -F"\t" '$2 != "Hypothetical protein"{count++}END{print count}')" >> "$outdir"/func_log_"$date".txt

printf "  \tUniprot/Hypothetical CDS\t=     %.3f\n" "$(tail -n +2 "$outdir"/genome_func_ids/*func_ids.tab | awk -F"\t" '$2 == "Hypothetical protein"{hypo++}$2 != "Hypothetical protein"{uniprot++}END{print uniprot/hypo}')" >> "$outdir"/func_log_"$date".txt

printf "  \tCDS with KO ID(S)\t\t=%'10d\n" "$(tail -n +2 "$outdir"/genome_func_ids/*func_ids.tab | awk -F"\t" '$3 != ""{count++}END{print count}')" >> "$outdir"/func_log_"$date".txt

printf "  \tCDS with COG ID(S)\t\t=%'10d\n" "$(tail -n +2 "$outdir"/genome_func_ids/*func_ids.tab | awk -F"\t" '$4 != ""{count++}END{print count}')" >> "$outdir"/func_log_"$date".txt

printf "  \tCDS with GO ID(S)\t\t=%'10d\n" "$(tail -n +2 "$outdir"/genome_func_ids/*func_ids.tab | awk -F"\t" '$5 != ""{count++}END{print count}')" >> "$outdir"/func_log_"$date".txt

printf "\n\n  Unique Uniprot flatfiles retrieved =%'6d\n\n" "$(ls "$outdir"/uniprot_info/ | wc -l)" >> "$outdir"/func_log_"$date".txt

printf "\n\n  Genome CDS Function Statistics:\n\n  %15s |%9s |%9s |%9s |%7s |%7s |%7s\n  -----------------------------------------------------------------------------\n" Type Total mean STD min max range >> "$outdir"/func_log_"$date".txt

printf "  %15s |%9d |%9.2f |%9.2f |%7d |%7d |%7d\n" $(awk -F"\t" '$2!=""&&$1!="Gene_id"{gsub("_[0-9]{5}",""); print $1}' "$outdir"/genome_func_ids/*.tab | uniq -c | awk '{sum+=$1; if ($1 > max) max = $1; if($1<min||NR==1) min = $1; val[NR]=$1;size=NR}END{mean=(sum/NR);for(x in val){diff+=(val[x]-(sum/size))^2};std=sqrt(diff/(size-1));print "Total\t"sum"\t"mean"\t"std"\t"min"\t"max"\t"max-min}') >> "$outdir"/func_log_"$date".txt

printf "  %15s |%9d |%9.2f |%9.2f |%7d |%7d |%7d\n" $(awk -F"\t" '$2=="Hypothetical protein"&&$1!="Gene_id"{gsub("_[0-9]{5}",""); print $1}' "$outdir"/genome_func_ids/*.tab | uniq -c | awk '{sum+=$1; if ($1 > max) max = $1; if($1<min||NR==1) min = $1;val[NR]=$1;size=NR}END{mean=(sum/NR);for(x in val){diff+=(val[x]-mean)^2};std=sqrt(diff/(size-1));print "Hypothetical\t"sum"\t"mean"\t"std"\t"min"\t"max"\t"max-min}') >> "$outdir"/func_log_"$date".txt

printf "  %15s |%9d |%9.2f |%9.2f |%7d |%7d |%7d\n" $(awk -F"\t" '$2!="Hypothetical protein"&&$1!="Gene_id"{gsub("_[0-9]{5}",""); print $1}' "$outdir"/genome_func_ids/*.tab | uniq -c | awk '{sum+=$1; if ($1 > max) max = $1; if($1<min||NR==1) min = $1;val[NR]=$1;size=NR}END{mean=(sum/size);for(x in val){diff+=(val[x]-mean)^2};std=sqrt(diff/(size-1));print "Uniprot\t"sum"\t"mean"\t"std"\t"min"\t"max"\t"max-min}') >> "$outdir"/func_log_"$date".txt

printf "  %15s |%9d |%9.2f |%9.2f |%7d |%7d |%7d\n" $(awk -F"\t" '$3!=""&&$1!="Gene_id"{gsub("_[0-9]{5}",""); print $1}' "$outdir"/genome_func_ids/*.tab | uniq -c | awk '{sum+=$1; if ($1 > max) max = $1; if($1<min||NR==1) min = $1;val[NR]=$1;size=NR}END{mean=(sum/size);for(x in val){diff+=(val[x]-mean)^2};std=sqrt(diff/(size-1));print "KO\t"sum"\t"mean"\t"std"\t"min"\t"max"\t"max-min}') >> "$outdir"/func_log_"$date".txt

printf "  %15s |%9d |%9.2f |%9.2f |%7d |%7d |%7d\n" $(awk -F"\t" '$4!=""&&$1!="Gene_id"{gsub("_[0-9]{5}",""); print $1}' "$outdir"/genome_func_ids/*.tab | uniq -c | awk '{sum+=$1; if ($1 > max) max = $1; if($1<min||NR==1) min = $1;val[NR]=$i;size=NR}END{mean=(sum/size);for(x in val){diff+=(val[x]-mean)^2};std=sqrt(diff/(size-1));print "COG\t"sum"\t"mean"\t"std"\t"min"\t"max"\t"max-min}') >> "$outdir"/func_log_"$date".txt

printf "  %15s |%9d |%9.2f |%9.2f |%7d |%7d |%7d\n" $(awk -F"\t" '$5!=""&&$1!="Gene_id"{gsub("_[0-9]{5}",""); print $1}' "$outdir"/genome_func_ids/*.tab | uniq -c | awk '{sum+=$1; if ($1 > max) max = $1; if($1<min||NR==1) min = $1;val[NR]=$1;size=NR}END{mean=(sum/size);for(x in val){diff+=(val[x]-mean)^2};std=sqrt(diff/(size-1));print "GO\t"sum"\t"mean"\t"std"\t"min"\t"max"\t"max-min}') >> "$outdir"/func_log_"$date".txt

printf "\n\n\n" >> "$outdir"/func_log_"$date".txt

################################################################################
## Extract function IDs for homologue gene clusters
################################################################################

# CDS sequences (ORFs) were clustered by get_homologues in the previous script for phylogenomic anlaysis.
# Build a functional profile for each CDS cluster using the genome_func_id.tab files produced above.

# first search if the get_homologues output exists

if [ -d "$outdir"/cluster_analysis ]; then

	mkdir -p "$outdir"/cluster_analysis/core_pan_func

	if [ ! -f "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab ]; then

		printf "Cluster\tUniProt IDs\tKO IDs\tCOG IDs\tGO IDs\n" > "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab

		for file in "$outdir"/cluster_analysis/prot_clst_intersection/*.faa; do

			file_only="${file##*/}"

			# get a list of the genes in the current cluster

			unset cluster_ids

			while read gene_id; do

				cluster_ids+=("$gene_id")

			done < <(grep ">" "$file" | awk -F"|" '{gsub(">ID:Prokka:",""); print $1}')

			printf "\nCluster $file_only has ${#cluster_ids[@]} associated genes\n"

	
			# clear variables from last loop
	
			unset uniprot_ids
			unset ko_ids
			unset cog_ids
			unset go_ids
	
			unset uniprot_list
			unset ko_list
			unset cog_list
			unset go_list

			# retrieve Uniprot, KO, COG, GO IDs for each gene in the cluster

			for i in ${!cluster_ids[@]}; do

				genome_id="${cluster_ids[$i]%_*}"
	
				# get Uniprot IDs
				uniprot_ids+=("$(grep "^${cluster_ids[$i]}" "$outdir"/func_id_files/"$genome_id"_func_ids.tab | awk -F"\t" '$2 != "Hypothetical protein"{print $2}')")

				# get KO IDs
				ko_ids+=("$(grep "^${cluster_ids[$i]}" "$outdir"/func_id_files/"$genome_id"_func_ids.tab | awk -F"\t" '{print $3}')")

				# get COG IDs
				cog_ids+=("$(grep "^${cluster_ids[$i]}" "$outdir"/func_id_files/"$genome_id"_func_ids.tab | awk -F"\t" '{print $4}')")
		
				# get GO IDs
				go_ids+=("$(grep "^${cluster_ids[$i]}" "$outdir"/func_id_files/"$genome_id"_func_ids.tab | awk -F"\t" '{print $5}')")

			done

			# remove duplicate IDs

			uniprot_ids=($(printf "%s\n" "${uniprot_ids[@]}" | sed 's/,/\n/g' | sort | uniq))
			ko_ids=($(printf "%s\n" "${ko_ids[@]}" | sed 's/,/\n/g' | sort | uniq))
			cog_ids=($(printf "%s\n" "${cog_ids[@]}" | sed 's/,/\n/g' | sort | uniq))
			go_ids=($(printf "%s\n" "${go_ids[@]}" | sed 's/,/\n/g' | sort | uniq))

			# sort non-redundant IDs into lists
	
			uniprot_list="$(printf "%s," "${uniprot_ids[@]}" | sed 's/,$//g')"
			ko_list="$(printf "%s," "${ko_ids[@]}" | sed 's/,$//g')"
			cog_list="$(printf "%s," "${cog_ids[@]}" | sed 's/,$//g')"
			go_list="$(printf "%s," "${go_ids[@]}" | sed 's/,$//g')"
	
			if [ -z "$uniprot_list" ]; then
		
				uniprot_list="Hypothetical protein"

			fi

			printf "%s\t%s\t%s\t%s\t%s\n" "${file_only/.faa/}" "$uniprot_list" "$ko_list" "$cog_list" "$go_list" >> "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab

		done

	fi

	## Use the cluster function IDs file to build a function ID file for core, soft-core, shell, and cloud genomes

	for file in "$outdir"/cluster_analysis/prot_clst_intersection/*_list.txt; do

		safefile="${file/\%/}"

		outfile="$(printf "${safefile/_list.txt/}" | awk -F"__" '{print $NF"_func_ids.tab"}')"

		if [ ! -f "$outdir"/func_id_files/"$outfile" ]; then

			total_clusters="$(wc -l "$file" | awk '{print $1}')"

			let "cluster_count = 0"


			printf "\n$outfile\t$cluster_count/$total_clusters"


			printf "Cluster\tUniProt ID\tKO ID\tCOG ID\tGO ID\n" > "$outdir"/cluster_analysis/core_pan_func/"$outfile"

			while read cluster; do

				let "cluster_count++"

				printf "\r$outfile\t$cluster_count/$total_clusters"

				grep "\<${cluster/.faa/}\>" "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab >> "$outdir"/cluster_analysis/core_pan_func/"$outfile"

			done < "$file"

		else

			printf "\n$outfile already exists in "$outdir"/func_id_files directory"

		fi

	done

	## copy core-pan genome func_id files to func_id_files

	cp "$outdir"/cluster_analysis/core_pan_func/*.tab "$outdir"/func_id_files/

	# add the core-pan cluster func files to the infile and filename arrays


	for file in "$indir"/cluster_analysis/core_pan_func/*.tab; do

		infile+=($(printf "${file##*/}"))
		filename+=($(printf "${file##*/}" | sed 's/_func_ids.tab//g'))

	done

fi

## Add Homologue func ID extraction to analysis report

printf "  Homologue cluster functional IDs\n  --------------------------------------------------------------------------------\n  Functional annotation for Homologue clusters produced by get_homologues.pl was\n  also determined. This was done by retrieving linked functional IDs (Uniprot, KO,\n  COG, GO) of the genes in each cluster from the (genome)_func_ids.tab files\n  produced previously. The resulting set of functional IDs for each cluster were\n  saved into new func_ids.tab files stored in the directory:\n\n  "$outdir"/cluster_analysis/core_pan_func \n\n\n"  >> "$outdir"/func_log_"$date".txt

printf "  \tProt Homologue annotation statistics\n  \t----------------------------------------------\n" >> "$outdir"/func_log_"$date".txt

printf "  \tTotal Homologues\t\t=%'10d\n" "$(awk 'END{print NR-1}' "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab)">> "$outdir"/func_log_"$date".txt

printf "  \tHypothetical Homologues\t\t=%'10d\n" "$(awk -F"\t" '$2 == "Hypothetical protein"{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with Uniprot ID(s)\t=%'10d\n" "$(awk -F"\t" '$2 != "Hypothetical protein"{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with KO ID(s)\t=%'10d\n" "$(awk -F"\t" '$3 != ""{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with COG ID(s)\t=%'10d\n" "$(awk -F"\t" '$4 != ""{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with GO ID(s)\t=%'10d\n" "$(awk -F"\t" '$5 != ""{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/prot-homologue_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "\n\n\n  The functional IDs for each class: core, softcore, shell and cloud genomes were\n  also determined.\n\n" >> "$outdir"/func_log_"$date".txt

printf "  \tCore Genome\n  \t----------------------------------------------\n"  >> "$outdir"/func_log_"$date".txt

printf "  \tTotal Homologues\t\t=%'10d\n" "$(awk 'END{print NR-1}' "$outdir"/cluster_analysis/core_pan_func/core_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHypothetical Homologues\t\t=%'10d\n" "$(awk -F"\t" '$2 == "Hypothetical protein"{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/core_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with Uniprot ID(s)\t=%'10d\n" "$(awk -F"\t" '$2 != "Hypothetical protein"{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/core_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with KO ID(s)\t=%'10d\n" "$(awk -F"\t" '$3 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/core_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with COG ID(s)\t=%'10d\n" "$(awk -F"\t" '$4 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/core_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with GO ID(s)\t=%'10d\n" "$(awk -F"\t" '$5 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/core_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "\n\n  \tSoftcore Genome\n  \t----------------------------------------------\n" >> "$outdir"/func_log_"$date".txt 

printf "  \tTotal Homologues\t\t=%'10d\n" "$(awk 'END{print NR-1}' "$outdir"/cluster_analysis/core_pan_func/softcore_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHypothetical Homologues\t\t=%'10d\n" "$(awk -F"\t" '$2 == "Hypothetical protein"{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/softcore_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with Uniprot ID(s)\t=%'10d\n" "$(awk -F"\t" '$2 != "Hypothetical protein"{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/softcore_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with KO ID(s)\t=%'10d\n" "$(awk -F"\t" '$3 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/softcore_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with COG ID(s)\t=%'10d\n" "$(awk -F"\t" '$4 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/softcore_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with GO ID(s)\t=%'10d\n" "$(awk -F"\t" '$5 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/softcore_func_ids.tab)" >> "$outdir"/func_log_"$date".txt


printf "\n\n  \tShell Genome\n  \t----------------------------------------------\n" >> "$outdir"/func_log_"$date".txt 

printf "  \tTotal Homologues\t\t=%'10d\n" "$(awk 'END{print NR-1}' "$outdir"/cluster_analysis/core_pan_func/shell_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHypothetical Homologues\t\t=%'10d\n" "$(awk -F"\t" '$2 == "Hypothetical protein"{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/shell_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with Uniprot ID(s)\t=%'10d\n" "$(awk -F"\t" '$2 != "Hypothetical protein"{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/shell_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with KO ID(s)\t=%'10d\n" "$(awk -F"\t" '$3 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/shell_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with COG ID(s)\t=%'10d\n" "$(awk -F"\t" '$4 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/shell_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with GO ID(s)\t=%'10d\n" "$(awk -F"\t" '$5 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/shell_func_ids.tab)" >> "$outdir"/func_log_"$date".txt



printf "\n\n  \tCloud Genome\n  \t----------------------------------------------\n" >> "$outdir"/func_log_"$date".txt 

printf "  \tTotal Homologues\t\t=%'10d\n" "$(awk 'END{print NR-1}' "$outdir"/cluster_analysis/core_pan_func/cloud_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHypothetical Homologues\t\t=%'10d\n" "$(awk -F"\t" '$2 == "Hypothetical protein"{count++}END{print count}' "$outdir"/cluster_analysis/core_pan_func/cloud_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with Uniprot ID(s)\t=%'10d\n" "$(awk -F"\t" '$2 != "Hypothetical protein"{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/cloud_func_ids.tab)" >> "$outdir"/func_log_"$date".txt

printf "  \tHomologues with KO ID(s)\t=%'10d\n" "$(awk -F"\t" '$3 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/cloud_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with COG ID(s)\t=%'10d\n" "$(awk -F"\t" '$4 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/cloud_func_ids.tab)" >> "$outdir"/func_log_"$date".txt
printf "  \tHomologues with GO ID(s)\t=%'10d\n" "$(awk -F"\t" '$5 != ""{count++}END{print count-1}' "$outdir"/cluster_analysis/core_pan_func/cloud_func_ids.tab)" >> "$outdir"/func_log_"$date".txt


################################################################################
## Analyse genome & homologue annotation statistics
################################################################################

# build a list of the genome ids 

	for i in ${!filename[@]}; do

		genome_id="$genome_id	${filename[$i]}"

	done

# Get genome annotation statistics

printf "\n\nTabulating genome annotation statistics "

printf "Genome\tTotal Genes\trRNA\ttRNA\tTotal CDS\tCDS w/ Uprot\tCDS w/ KO\tCDS w/ COG\t CDS w/ GO\n" > "$outdir"/Genome_annotation_stats.tab

for i in ${!infile[@]}; do

#	filename="${infile[$i]/_func_ids.tab/}"

	if [ -f "$outdir"/prokka_annotation/"${filename[$i]}".txt ]; then

		total_genes="$(grep "^gene" "$outdir"/prokka_annotation/"${filename[$i]}".txt | awk '{print $2}')"
		total_trna="$(grep "^tRNA" "$outdir"/prokka_annotation/"${filename[$i]}".txt | awk '{print $2}')"
		total_rrna="$(grep "^rRNA" "$outdir"/prokka_annotation/"${filename[$i]}".txt | awk '{print $2}')"

		if [ -z "$total_rrna" ]; then

			total_rrna="0"

		fi

		if [ -z "$total_trna" ]; then

			total_trna="0"

		fi

	else

		total_genes="NA"
		total_trna="NA"
		total_rrna="NA"

	fi

	total_cds="$(tail -n +2 "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab | wc -l)"
	uprot_anno_count="$(awk -F"\t" '$2!="Hypothetical protein"{print}' "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab | tail -n +2 | wc -l)"
	ko_anno_count="$(awk -F"\t" '$3!=""{print}' "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab | tail -n +2 | wc -l)"
	cog_anno_count="$(awk -F"\t" '$4!=""{print}' "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab | tail -n +2| wc -l)"
	go_anno_count="$(awk -F"\t" '$5!=""{print}' "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab | tail -n +2 | wc -l)"


	tsv_genes="$tsv_genes	$total_genes"
	tsv_trna="$tsv_trna	$total_trna"
	tsv_rrna="$tsv_rrna	$total_rrna"
	tsv_cds="$tsv_cds	$total_cds"
	tsv_uprot="$tsv_uprot	$uprot_anno_count"
	tsv_ko="$tsv_ko	$ko_anno_count"
	tsv_cog="$tsv_cog	$cog_anno_count"
	tsv_go="$tsv_go	$go_anno_count"

	printf "${filename[$i]}\t$total_genes\t$total_rrna\t$total_trna\t$total_cds\t$uprot_anno_count\t$ko_anno_count\t$cog_anno_count\t$go_anno_count\n" >> "$outdir"/Genome_annotation_stats.tab

	printf "."

done

# build a transposed table of annotation statistics

printf "Stat$genome_id\nTotal genes$tsv_genes\nrRNA$tsv_rrna\ntRNA$tsv_trna\nCDS$tsv_cds\nCDS w/ UPROT$tsv_uprot\nCDS w/ KO$tsv_ko\nCDS w/ COG$tsv_cog\nCDS w/ GO$tsv_go" > "$outdir"/Genome_annotation_stats_transverse.tsv

printf "DONE\n"

## Add annotation stat info to analysis report

printf "\n\n\n  Genome annotation statistics table\n  --------------------------------------------------------------------------------\n  The PROKKA annotation characteristics and functional ID statistic for the genomes\n  and pancore genome groups were tabulated. The table and transposed table have\n  been saved the the files listed below:\n\n  "$outdir"/Genome_annotation_stats.tab\n  "$outdir"/Genome_annotation_stats_transversed.tsv\n\n\n" >> "$outdir"/func_log_"$date".txt


################################################################################
## Retrieve and Prepare COG and KO functional information
################################################################################

# make separate output directories for each functional analysis ( make output cleaner )

mkdir -p "$outdir"/COG_analysis
mkdir -p "$outdir"/KEGG_analysis
mkdir -p "$outdir"/GO_analysis

# Retreive the COG names file if it doesn't already exist

if [ ! -f "$outdir"/COG_analysis/cog_function_names.tab ]; then

	wget -q ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab -O "$outdir"/COG_analysis/cog_function_names.tab

fi

declare -A cog_cat
declare -A cog_func

printf "\nLoading COG functions into memory\n\n"

while read id cat func; do

	cog_cat[$id]="$cat"
	cog_func[$id]="$func"

done < <(tail -n +2  "$outdir"/COG_analysis/cog_function_names.tab)

## read through kegg orthology file and assign values to associative arrays

if [ ! -f "$outdir"/KEGG_analysis/Kegg_KO_func_map.tab ]; then

   	if [ -f "$kegg_file" ]; then
	
	   KO_number="$(grep -c "^D" "$kegg_file")"
   	
	   let "KO_count = 0"
   	
	   printf "\nConverting Kegg Orthology map to tab delimited file\n"
	   printf "$KO_count/$KO_number"
   
		declare -A kegg_cat
   		declare -A kegg_subcat
   		declare -A kegg_pathway
   		declare -A kegg_short
   		declare -A kegg_name
   		declare -A kegg_ec
   		printf "KO_ID\tCategory\tSub_category\tPathway\tFunction\tshort_name\tEC_ID\n" > "$outdir"/KEGG_analysis/Kegg_KO_func_map.tab

  		 while read part_a part_b; do
				
	  		if [[ $part_a == A* ]]; then

				cat_name="$(printf "$part_b" | sed 's/^ //g')"

	  		elif [[ $part_a == B ]]; then

				subcat_name="$(printf "$part_b" | awk '{out=""; for(x=2;x<=NF;x++){out=out" "$x}; print out}' | sed 's/^ //g')"
		
	  		elif [[ $part_a == C ]]; then

				pathway_name="$(printf "$part_b" | awk '{out=""; for(x=2;x<=NF;x++){out=out" "$x}; print out}' | sed 's/^ //g')"

	  		elif [[ $part_a == D ]]; then

				let "KO_count++"
		 		printf "\r$KO_count/$KO_number"

				kegg_KO="$(printf "$part_b" | awk '{print $1}')"

		 		short_name="$(printf "$part_b" | awk -F";" '{print $1}' | awk '{out="";for(x=2;x<=NF;x++){out=out$x}; print out}')"

		 		func_name="$(printf "$part_b" | awk -F";" '{print $2}' | awk -F"[" '{print $1}' | sed 's/^ //g')"

		 		ec_id="$(printf "$part_b" | awk -F"EC:" '{gsub("]",""); gsub(" ",","); print $2}')"


		 		kegg_cat[$kegg_KO]="$cat_name"
		 		kegg_subcat[$kegg_KO]="$subcat_name"
		 		kegg_pathway[$kegg_KO]="$pathway_name"
		 		kegg_short[$kegg_KO]="$short_name"
		 		kegg_name[$kegg_KO]="$func_name"
		 		kegg_ec[$kegg_KO]="$ec_id"

		 		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$kegg_KO" "${kegg_cat[$kegg_KO]}" "${kegg_subcat[$kegg_KO]}" "${kegg_pathway[$kegg_KO]}" "${kegg_name[$kegg_KO]}" "${kegg_short[$kegg_KO]}" "${kegg_ec[$kegg_KO]}" >> "$outdir"/KEGG_analysis/Kegg_KO_func_map.tab

	  		fi

   		done < "$kegg_file"

   		printf "\n\n"

	 else

		printf "\nNo Kegg Orthology map file detected, will skip KO function assignments\n"

	 fi

else

	printf "KEGG_KO_func_map.tab file found in $outdir/KEGG_analysis/\n\n"

fi

################################################################################
## Link genes to KO and COG function data
################################################################################

for i in ${!infile[@]}; do

#	filename="${infile[$i]%.*}"

	if [ ! -f "$outdir"/COG_analysis/"${filename[$i]}"_cog_function.tab ] || [ ! -f "$outdir"/KEGG_analysis/"${filename[$i]}"_KO_function.tab ]; then

# start with KO annotations

	printf "\nRetrieving KO and COG data for genome %s\n" ${filename[$i]}

	printf "Gene_id\tCOG_category\tCOG_function\n" > "$outdir"/COG_analysis/"${filename[$i]}"_cog_function.tab

	printf "Gene_ID\tKO_ID\tCategory\tSub_category\tPathway\tFunction\tShort_name\tEC_ID\n" > "$outdir"/KEGG_analysis/"${filename[$i]}"_KO_function.tab


	IFS='|'

	#gene_number="$(grep -c "^${filename[$i]}" "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab)"

	#if [ "$gene_number" = 0 ]; then

		gene_number="$(tail -n +2 "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab | wc -l)"

	#fi

	let "gene_count = 0"

	printf "$gene_count/$gene_number"

	while read gene_id uprot_id KO_id COG_id GO_id; do

		let "gene_count++"
		printf "\r$gene_count/$gene_number"
		
		unset ko_array

		if [[ $KO_id == *,* ]]; then

			IFS=',' read -r -a ko_array <<< "$KO_id"

			IFS='|'

		else

			ko_array+=("$KO_id")

		fi

		for x in ${!ko_array[@]}; do

			ko_id="${ko_array[$x]}"

			if [ -n "$ko_id" ]; then

				while read ko_data; do

					printf "$gene_id\t$ko_data\n" >> "$outdir"/KEGG_analysis/"${filename[$i]}"_KO_function.tab

				done < <(grep "^$ko_id" "$outdir"/KEGG_analysis/Kegg_KO_func_map.tab)

			else

				printf "$gene_id\n" >> "$outdir"/KEGG_analysis/"${filename[$i]}"_KO_function.tab

			fi

		done

# Complete by processing COG annotations

		unset cog_array

		if [[ $COG_id == *,* ]]; then

			IFS=',' read -r -a cog_array <<< "$COG_id"

			IFS='|'

		else

			cog_array+=("$COG_id")

		fi

		unset gene_cog_cat
		unset gene_cog_func
		unset cog_cat_list
		unset cog_func_list

		for x in ${!cog_array[@]}; do

				cog_id="${cog_array[$x]}"

			if [ -n "$cog_id" ]; then

				gene_cog_cat+="${cog_cat[$cog_id]}"
				gene_cog_func+="[${cog_func[$cog_id]}],"

				unset cog_cat_list

				if [ -n "$gene_cog_cat" ]; then

					cog_cat_list="$( printf "%s" ${gene_cog_cat[@]})"

				fi

				if [ -n "$gene_cog_func" ]; then

					cog_func_list="$( printf "%s" ${gene_cog_func[@]})"

				##	printf "$cog_func_list\n"

				fi
			
			fi

		done

		cog_cat_list="$(echo $cog_cat_list | sed 's/\.*/ /g;s/^ //g' | awk '{for(i=1;i<=NF;i++){if(i in a == 0)a[$i]}for(i in a) b = i b}END{print b}')" # Remove duplicate cog function letters assigned to each gene

		printf "%s\t%s\t%s\n" $gene_id ${cog_cat_list%%,} ${cog_func_list%%,}  >> "$outdir"/COG_analysis/"${filename[$i]}"_cog_function.tab

	done < <(sed 's/\t/|/g'  "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab | tail -n +2)

fi

done

################################################################################
## Compute COG category totals for genomes
################################################################################

# retrieve COG function list

if [ ! -f "$outdir"/COG_analysis/COG_func_list.txt ]; then

	printf "\nRetrieving list of COG function categories\n"

	wget -q -O "$outdir"/COG_analysis/COG_func_list.txt ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab

fi

# for each genome count number of genes per COG category

printf "\nCompiling COG function counts "

IFS=' '

for i in ${!infile[@]}; do

#	filename="${infile[$i]%.*}"

	unset gen_cog_cat

	while read category; do

		gen_cog_cat+=( $category )

	done < <(tail -n +2  "$outdir"/COG_analysis/"${filename[$i]}"_cog_function.tab | awk -F"\t" '$2!=""{print $2}' | sed 's/\.*/ /g;s/^ //g')

	printf "%s\n" ${gen_cog_cat[@]} | sort | uniq -c | awk '{print $2"\t"$1}'  > "$outdir"/COG_analysis/"${filename[$i]}"_COG_cat_count.tab

	printf "."

done

printf "DONE\n"

# Use the COG function (category) list to build a table of all the genomes
printf "\nBuilding COG function table "
printf "COG\tFunction$genome_id\n" > "$outdir"/COG_analysis/COG_category_table.tab
printf "Genome\tCOG\tFunction\tcount\n" > "$outdir"/COG_analysis/COG_genome_function_counts.tab

IFS=$'\t'
while read cog_cat cat_name; do

	unset all_cog_counts

	for i in ${!infile[@]}; do

	#	filename="${infile[$i]%.*}"
	
		cog_count="$(awk -v cog_cat="$cog_cat" -F"\t" '$1==cog_cat{print $2}' "$outdir"/COG_analysis/"${filename[$i]}"_COG_cat_count.tab)"

		if [ -z "$cog_count" ]; then

			cog_count="0"

		fi

		all_cog_counts="$all_cog_counts	$cog_count"

		printf "${filename[$i]}\t$cog_cat\t$cat_name($cog_cat)\t$cog_count\n" >> "$outdir"/COG_analysis/COG_genome_function_counts.tab

	done

	printf "$cog_cat\t$cat_name$all_cog_counts\n" >> "$outdir"/COG_analysis/COG_category_table.tab

	printf "."


done < <(tail -n +2 "$outdir"/COG_analysis/COG_func_list.txt)

printf "DONE\n\n"

# split table between genomes and homologue sets
printf "Splitting COG category counts genomes/pancore genome\n\n"

awk -F"\t" '{for(i=1;i<=(NF-6);i++) printf $i"\t"; print $i}' "$outdir"/COG_analysis/COG_category_table.tab > "$outdir"/COG_analysis/COG_genome_category_table.tab

printf "COG category table genomes only: "$outdir"/COG_analysis/COG_genome_category_table.tab\n\n"

awk -F"\t" '{printf "%s\t%s\t", $1, $2; for(i=(NF-4);i<=NF;i++) printf $i"\t"; print $i}' "$outdir"/COG_analysis/COG_category_table.tab > "$outdir"/COG_analysis/COG_pancore_category_table.tab

printf "COG category table pancore genome: "$outdir"/COG_analysis/COG_pancore_category_table.tab\n\n"

## Add COG analysis info to analysis report
printf "\n\n\n  Genome COG category analysis\n  --------------------------------------------------------------------------------\n  The PROKKA annotated genes that were assigned COG IDs via UniProt flatfiles were\n  linked to their respective COG categorical and functional information. The COG\n  functional annotation data was obtained from the cognames2003-2014.tab file\n  retrieved from the NCBI ftp server and saved under the name \"cog_function_names.tab\".\n  The Total number of different COGs as well as the number of uniq COGs identified\n  in the set of genomes is show under \"COG ID statistics\" below. Note that a single\n  COG ID can be assigned to multiple categories and a single gene could be assigned\n  to more than one COG ID. COGs are sorted into functional categories represented\n  by both a one letter code (A-Z) and a more descriptive name. A statistical table\n  of COG category representation in genomes as well as a table represented COG\n  category counts in pangenome composition are shown below.\n\n\n" >> "$outdir"/func_log_"$date".txt

printf "  \tCOG ID statistics\n  \t-------------------\n" >> "$outdir"/func_log_"$date".txt

Total_cogs="$(awk 'NR>1{count++}END{print count}' "$outdir"/COG_analysis/cog_function_names.tab)"

uniq_cogs="$(awk '$3~/^COG/{print $3}' "$outdir"/uniprot_info/*.txt | sort | uniq | wc -l)"
cog_rep=$(bc -l <<< "scale=4;$uniq_cogs / $Total_cogs * 100")

printf "  \tTotal unique COG IDs\t\t=%'9d\n" $Total_cogs >> "$outdir"/func_log_"$date".txt

printf "  \tUnique COG IDs in genomes\t=%'9d\n" $uniq_cogs >> "$outdir"/func_log_"$date".txt

printf "  \tUnique COG representation\t=%'9.2f %%\t(Total COGs / Genome COGs)\n\n\n" $cog_rep >> "$outdir"/func_log_"$date".txt

printf "  COG analysis files ("$outdir"/COG_analysis/):\n\n  cog_function_names.tab\t\t->\tCOG names linked to IDs & Categories\n  COG_func_list.tab\t\t\t->\tList of COG categories\n  COG_genome_function_counts.tab\t->\tGenes per category per genome long format\n  COG_category_table.tab\t\t->\tTable of COG category counts per genome\n  COG_genome_category_table.tab\t\t->\tGenome only category table\n  COG_pancore_category_table.tab\t->\tPan-core genome only category table\n  'GENOME'_COG_function.tab\t\t->\tCOG information linked to individual genes\n  'GENOME'_COG_cat_count.tab\t\t->\tNo. genes per COG category\n\n"  >> "$outdir"/func_log_"$date".txt



printf "  Genome COG category statistics\n  ------------------------------\n\n" >> "$outdir"/func_log_"$date".txt  

awk -F"\t" 'NR>1{for(i=3;i<=NF;i++)total[NR]+=$i; printf total[NR]"\t"; print}' "$outdir"/COG_analysis/COG_genome_category_table.tab | sort -k 1 -nr | awk -F"\t" 'BEGIN{printf "%5s |%8s |%8s |%7s |%7s |%5s |%5s |%6s | %s\n  ----------------------------------------------------------------------------------------------------\n", "COG", "Total", "Mean", "STD", "SEM", "Min", "Max", "Range", "Function"}{for(i=4;i<=NF;i++){sum[NR]+=$i;mean[NR]=(sum[NR]/(NF-3));if($i>max[NR])max[NR]=$i;if($i<min[NR]||i==4)min[NR]=$i};val[i]=$i;for(x in val){diff[NR]+=($x-(mean[NR]))^2}; std[NR]=sqrt(diff[NR]/(NF-4));sem[NR]=(std[NR]/sqrt(NF-3)); printf "%5s |%8d |%8.2f |%7.2f |%7.2f |%5d |%5d |%6d | %s\n", $2, sum[NR], (sum[NR]/(NF-3)), std[NR], sem[NR], min[NR], max[NR], (max[NR]-min[NR]), $3}' >> "$outdir"/func_log_"$date".txt

# COG pancore genome 
printf "\n\n  COG pancore genome\n  ------------------\n\n" >> "$outdir"/func_log_"$date".txt

tail -n +2 "$outdir"/COG_analysis/COG_pancore_category_table.tab | sort -t $'\t' -k5 -nr | awk -F"\t" 'BEGIN{printf "%7s |%7s |%9s |%7s |%7s |%7s | %s\n  -----------------------------------------------------------------------------------------\n", "COG", "Core", "Softcore", "Shell", "Cloud", "Total", "Function"}{printf "%7s |%7d |%9d |%7d |%7d |%7d | %s\n", $1, $4, $7, $6, $3, $5, $2}' >> "$outdir"/func_log_"$date".txt


################################################################################
## Compute KO Sub-category totals for genomes
################################################################################

# First for each genome get counts for each sub category and save to file
for i in ${!infile[@]}; do

	awk -F"\t" 'NR>1&&$4!=""{ print $4}' "$outdir"/KEGG_analysis/"${filename[$i]}"_KO_function.tab | sort | uniq -c | awk '{for(i=2;i<=NF;i++){subcat[NR] = subcat[NR]" "$i}; print subcat[NR]"\t"$1}' | sed 's/^ //g' > "$outdir"/KEGG_analysis/"${filename[$i]}"_KO_subcat_count.tab

done

# Next produce a list of the subcategories to use as the table framework
awk -F"\t" 'NR>1{print $3}' "$outdir"/KEGG_analysis/Kegg_KO_func_map.tab | sort | uniq > "$outdir"/KEGG_analysis/KEGG_KO_subcat_list.txt


printf "\nBuilding KEGG subcategory table "

# Write out header for KEGG subcategory table
printf "Subcategory$genome_id\n" > "$outdir"/KEGG_analysis/KEGG_subcat_table.tab
printf "Genome\tSubcategory\tcount\n" > "$outdir"/KEGG_analysis/KEGG_genome_function_counts.tab

# Write out values (counts) to subcategory table
IFS=$'\t'
while read subcat; do

	unset all_ko_counts

	for i in ${!infile[@]}; do

	#	filename="${infile[$i]%.*}"
	
		ko_count="$(awk -v subcat="$subcat" -F"\t" '$1==subcat{print $2}' "$outdir"/KEGG_analysis/"${filename[$i]}"_KO_subcat_count.tab)"

		if [ -z "$ko_count" ]; then

			ko_count="0"

		fi

		all_ko_counts="$all_ko_counts	$ko_count"

		printf "${filename[$i]}\t$subcat\t$ko_count\n" >> "$outdir"/KEGG_analysis/KEGG_genome_function_counts.tab

	done

	printf "$subcat$all_ko_counts\n" >> "$outdir"/KEGG_analysis/KEGG_subcat_table.tab

	printf "."


done < "$outdir"/KEGG_analysis/KEGG_KO_subcat_list.txt


printf "DONE\n\n"

printf "\nALL KEGG subcat counts in file: "$outdir"/KEGG_analysis/KEGG_subcat_table.tab\n\n"


# split table between genomes and homologue sets

printf "Splitting KEGG subcat tabe at genomes/pancore genome\n\n"

awk -F"\t" '{for(i=1;i<=(NF-6);i++) printf $i"\t"; print $i}' "$outdir"/KEGG_analysis/KEGG_subcat_table.tab > "$outdir"/KEGG_analysis/KEGG_genome_subcat_table.tab

printf "KEGG subcat genomes only: "$outdir"/KEGG_analysis/KEGG_genome_subcat_table.tab\n\n"

awk -F"\t" '{printf "%s\t", $1; for(i=(NF-4);i<=NF;i++) printf $i"\t"; print $i}' "$outdir"/KEGG_analysis/KEGG_subcat_table.tab > "$outdir"/KEGG_analysis/KEGG_pancore_subcat_table.tab

printf "KEGG subcat pancore only: "$outdir"/KEGG_analysis/KEGG_pancore_subcat_table.tab\n\n"


## Add KEGG analysis info to analysis report

printf "\n\n\n  Genome KEGG annotation analysis\n  --------------------------------------------------------------------------------\n" >> "$outdir"/func_log_"$date".txt

printf "  The KO IDs that were linked to genes and homologue clusters previously were linked\n  to KEGG (Kyoto Encyclopedia of Genes and Genomes) functional annotation information.\n  This allows KO linked genes to be identified at functional category, subcategory\n  and even individual metabolic pathway levels.\n\n" >> "$outdir"/func_log_"$date".txt

printf "  \tKEGG orthology statistics\n  \t-------------------------\n" >> "$outdir"/func_log_"$date".txt 

Total_kos="$(awk 'NR>1{count++}END{print count}' "$outdir"/KEGG_analysis/Kegg_KO_func_map.tab)"

uniq_kos="$(awk '$2~/^KO;/{print $3}' "$outdir"/uniprot_info/*.txt | sort | uniq | wc -l)"
ko_rep=$(bc -l <<< "scale=4;$uniq_kos / $Total_kos * 100")

printf "  \tTotal unique KO IDs\t\t=%'9d\n" "$Total_kos" >> "$outdir"/func_log_"$date".txt

printf "  \tNo. unique KO IDs in genomes\t=%'9d\n" "$uniq_kos" >> "$outdir"/func_log_"$date".txt

printf "  \tGlobal KO representation\t=%'9.2f %%\t(Genome / Total KO IDs)\n\n\n" "$ko_rep" >> "$outdir"/func_log_"$date".txt

printf "  KEGG analysis files ("$outdir"/KEGG_analysis/):\n\n  KEGG_KO_func_map.tab\t\t\t->\ttab delim KEGG KO mapping file\n  KEGG_genome_function_counts.tab\t->\tKEGG subcategory counts per genome\n  KEGG_subcat_table.tab\t\t\t->\tTable of subcategory counts per genome\n  KEGG_genome_subcat_table.tab\t\t->\tGenome only subcat table\n  KEGG_pancore_subcat_table.tab\t\t->\tPan-core genome only subcat table\n  'GENOME'_KO_function.tab\t\t->\tKEGG information linked to individual genes\n  'GENOME'_KO_subcat_count.tab\t\t->\tNo. KO IDs per subcategory\n\n" >> "$outdir"/func_log_"$date".txt  

printf "  Genome KEGG subcategory statistics\n  ----------------------------------\n\n" >> "$outdir"/func_log_"$date".txt

awk -F"\t" 'NR>1{for(i=2;i<=NF;i++)total[NR]+=$i; printf total[NR]"\t"; print}' "$outdir"/KEGG_analysis/KEGG_genome_subcat_table.tab | sort -t $'\t' -k 1 -nr | awk -F"\t" 'BEGIN{printf "  %7s |%8s |%8s |%7s |%5s |%5s |%6s | %s\n  --------------------------------------------------------------------------------\n", "Total", "Mean", "STD", "SEM", "Min", "Max", "Range", "Subcategory"}{for(i=3;i<=NF;i++){sum[NR]+=$i;if($i>max[NR])max[NR]=$i;if($i<min[NR]||i==3)min[NR]=$i;val[i]=$i};category=$2;mean[NR]=(sum[NR]/(NF-2));for(x in val){diff[NR]+=($x-mean[NR])^2}std[NR]=sqrt(diff[NR]/(NF-3));ste[NR]=(std[NR]/sqrt(NF-2))}{printf "  %7d |%8.2f |%8.2f |%7.2f |%5d |%5d |%6d | %s\n", (sum[NR]), (mean[NR]), (std[NR]), (ste[NR]), (min[NR]), (max[NR]), (max[NR]-min[NR]), category}' >> "$outdir"/func_log_"$date".txt 

printf "\n\n  KEGG pancore genome\n  -------------------\n\n" >> "$outdir"/func_log_"$date".txt

tail -n +2 "$outdir"/KEGG_analysis/KEGG_pancore_subcat_table.tab | sort -t $'\t' -k5 -nr | awk -F"\t" 'BEGIN{printf "%7s |%9s |%7s |%7s |%7s | %s\n  -----------------------------------------------------------------------------------------\n", "Core", "Softcore", "Shell", "Cloud", "Total", "Function"}{printf "%7d |%9d |%7d |%7d |%7d | %s\n", $3, $6, $5, $2, $4, $1}' >> "$outdir"/func_log_"$date".txt

printf "\n\n" >> "$outdir"/func_log_"$date".txt


################################################################################
## Build krona graphs from Kegg KO annotations
################################################################################

mkdir -p "$outdir"/KEGG_analysis/krona_KO_graphs

for i in ${!infile[@]}; do

#	filename="${infile[$i]%.*}"

	printf "Building krona graph for $filename from KO annotation\n"

	awk -F"\t" '$2 != ""{gsub("^ ","",$4); print $3"|\t"$4"|\t"$5}' "$outdir"/KEGG_analysis/"${filename[$i]}"_KO_function.tab | tail -n +2 | sort | uniq -c | awk -F" " '{gsub("$","\t",$1); print $0}' | sed 's/| /\t/g' > "$outdir"/KEGG_analysis/krona_KO_graphs/"${filename[$i]}"_KO_krona.txt

	ktImportText -o "$outdir"/KEGG_analysis/krona_KO_graphs/"${filename[$i]}"_KO_krona.html "$outdir"/KEGG_analysis/krona_KO_graphs/"${filename[$i]}"_KO_krona.txt

done

printf "\n\n"

## ADD KEGG Krona info to Analysis report

printf "  Genome KEGG Orthology Krona Graphs\n  --------------------------------------------------------------------------------\n  The hierarchical nature of KEGGs KO functional annotation makes it possible to\n  produce Krona graph representations. Therefore it was possible to build a Krona\n  graph representing the KEGG functional distribution for each genome and pan-core\n  category. The krona graphs were saved to the directory:\n\n  "$outdir"/KEGG_analysis/krona_KO_graphs/:\n\n  'Genome'_KO_krona.html\t->\tKrona graph\n  'Genome'_KO_krona.txt\t\t->\tFile used to build krona graph\n" >> "$outdir"/func_log_"$date".txt 


################################################################################
## Retrieve and prepare GO ID information
################################################################################

## Retrieve Gene Ontology (GO) id information (go.obo) convert to a delimited file

if [ ! -f "$outdir"/GO_analysis/go.obo ]; then

	printf "\nRetrieving Gene Ontology Go ID annotation data\n"

	wget -q -O "$outdir"/GO_analysis/go.obo http://purl.obolibrary.org/obo/go.obo
	
fi

if [ -f "$outdir"/GO_analysis/go.obo ]; then
	
	if [ ! -f "$outdir"/GO_analysis/GO_term_annotation.tab ]; then

		awk -F"\n" 'BEGIN{RS="\n\n"}$1~/\[Term\]/{for(i=2;i<=NF;i++){if($i ~ /^id:/)goid[NR]=substr($i,5);if($i ~ /^name:/)name[NR]=substr($i,7);if($i ~ /^namespace:/)ns[NR]=substr($i,11);if($i ~ /def:/)def[NR]=substr($i,5)}}END{for(i in name) print goid[i]"\t"name[i]"\t"ns[i]"\t"def[i]}' "$outdir"/GO_analysis/go.obo > "$outdir"/GO_analysis/GO_term_annotation.tab

		# awk 'BEGIN{RS="\n\n"}{print}' "$outdir"/GO_analysis/go.obo | sed 's/\[Term\]/,,,,,/g' | awk -F"\n" 'BEGIN{RS=",,,,,"}$2~/id:/{gsub("id: |name:|namespace:","");print $2"\t"$4"\t"$3}' > "$outdir"/GO_analysis/GO_ids_anno.tab
		# awk -F"\n" 'BEGIN{RS="\n\n"}$1~/\[Term\]/{gsub("name: |namespace: |id: |def: ","");print $2"\t"$4"\t"$3"\t"$5}' "$outdir"/GO_analysis/go.obo > "$outdir"/GO_analysis/GO_term_annotation.tab

	fi

else

	printf "\nCould not retrieve 'go.obo' file from http://purl.obolibrary.org/obo/go.obo\n"

fi

printf "\n"

## load GO id data into memory

printf "Load Go data into memory\t"

declare -A go_namespace
declare -A go_name
#declare -A go_def

IFS=$'\t'
while read id namespace name def; do

	go_namespace[$id]="$namespace"
	go_name[$id]="$name"
#	go_def[$id]="$def"

done < "$outdir"/GO_analysis/GO_term_annotation.tab

printf "${#go_name[@]} GO ids parsed\n\n"

################################################################################
## Analyse GO IDs from genomes
################################################################################

## Compile GO ids from genomes into table format

if [ -f "$outdir"/GO_analysis/all_go_ids.txt ]; then

	rm "$outdir"/GO_analysis/all_go_ids.txt

fi

for i in ${!infile[@]}; do

#	filename="${infile[$i]%.*}"

	printf "Count GO ids from ${filename[$i]}\n"

	awk -F"\t" '$5 != ""{print $5}' "$outdir"/func_id_files/"${filename[$i]}"_func_ids.tab | sed 's/,/\n/g' | sort | uniq -c |sort -nr | awk '{print $2"\t"$1}' > "$outdir"/GO_analysis/"${filename[$i]}"_go_count.txt

	awk -F"\t" '$1~/GO:/{print $1}' "$outdir"/GO_analysis/"${filename[$i]}"_go_count.txt >> "$outdir"/GO_analysis/all_go_ids.txt

	printf "GO_id\tCount\tNamespace\tName\n" > "$outdir"/GO_analysis/"${filename[$i]}"_go_annotated.tab
	
	while read go_id count; do

		#   Get GO namespace and name using associative arrays

		printf "$go_id\t$count\t${go_namespace[$go_id]}\t${go_name[$go_id]}\n" >> "$outdir"/GO_analysis/"${filename[$i]}"_go_annotated.tab

	done < "$outdir"/GO_analysis/"${filename[$i]}"_go_count.txt
	
done

# Build a non-redundant list of GO ids from all genomes and produce a GO ID table of all genome counts

# the non-redundant list of go-ids can be loaded directly into an array

if [ ! -f "$outdir"/GO_analysis/GO_term_count_table.tab ]; then

	process_count="0"

	unset go_id_list

	printf "\nBuilding table of GO_id counts from all genomes\n"

	while read go_id; do

		go_id_list+=($go_id)

	done < <(sort $outdir/GO_analysis/all_go_ids.txt | uniq)

	printf "${#go_id_list[@]} Uniq GO IDs found\n"

	printf "%i/%i" "$process_count" "${#go_id_list[@]}"

	# print header for file GO_term_count_table.tab
	printf "GO_id\tName\tNamespace$genome_id\n" > "$outdir"/GO_analysis/GO_term_count_table.tab

	# process via the list of GO ids

	for i in ${!go_id_list[@]}; do

		let "process_count++"

		printf "\r$process_count/${#go_id_list[@]}"

		go_id="${go_id_list[$i]}"

		unset gen_go_count

		total_count=0

		for x in ${!infile[@]}; do

#			filename="${infile[$x]%%.*}"

			go_count="$(grep "$go_id" "$outdir"/GO_analysis/"${filename[$x]}"_go_count.txt | awk -F"\t" '{print $2}')"
		
			if [ -z "$go_count" ]; then

				go_count="0"

			fi

			gen_go_count="$gen_go_count	$go_count"

		#total_count=$(( total_count + go_count ))

		done

		#go_count_list=($gen_go_count)

		# go_avg_std="$(printf "$gen_go_count" | awk '{for(i=1;i<=NF;i++){val[i] = $i; sum += $i}}{avg=sum/NF}END{for(i in val){diff+=($i-avg)^2} print sum"\t"avg"\t"sqrt(diff/(NF-1))}')"

		printf "${go_id_list[$i]}\t${go_namespace[$go_id]}\t${go_name[$go_id]}$gen_go_count\n" >> "$outdir"/GO_analysis/GO_term_count_table.tab

	done


printf "\nGO id count table completed\n\n"

fi

awk -F"\t" '{for(i=1;i<=(NF-6);i++) printf $i"\t"; print $i}' "$outdir"/GO_analysis/GO_term_count_table.tab > "$outdir"/GO_analysis/GO_term_genome_count_table.tab 


awk -F"\t" '{printf "%s\t%s\t%s\t", $1, $2, $3; for(i=(NF-4);i<=NF;i++) printf $i"\t"; print $i}' "$outdir"/GO_analysis/GO_term_count_table.tab > "$outdir"/GO_analysis/GO_term_pancore_count_table.tab

## Add GO analysis info to analysis report

printf "\n\n\n  Genome GO (Gene Ontology) analysis\n  --------------------------------------------------------------------------------\n  Go terms were retrieved for genes linked to UniProt IDs via PROKKA annotation.\n  The Gene Ontology project assigns different biological classes, functions and\n  concepts to GO Terms (IDs). The GO terms are linked to annotated biological\n  sequences (genes/proteins) with experimental evidence to support the terms\n  association. In general a single gene sequence from an annotated database\n  (e.g. Uni/Swiss prot) will in general be assigned several GO terms to represent\n  different components of its biological/functional characteristics.\n\n\n" >> "$outdir"/func_log_"$date".txt 


printf "  \tGO Term Representation\n  \t---------------------\n" >> "$outdir"/func_log_"$date".txt

Total_gos="$(awk '{count++}END{print count}' "$outdir"/GO_analysis/GO_term_annotation.tab)"
uniq_gos="$(sort "$outdir"/GO_analysis/all_go_ids.txt | uniq | wc -l)"

go_rep=$(bc -l <<< "scale=4;$uniq_gos / $Total_gos * 100")

printf "  \tTotal unique GO Terms\t\t=%'9d\n" "$Total_gos" >> "$outdir"/func_log_"$date".txt

printf "  \tUnique GO Terms from genomes\t=%'9d\n" "$uniq_gos" >> "$outdir"/func_log_"$date".txt

printf "  \tUnique GO term representation\t=%'9.2f %%\t(Genome / Total GO terms)\n\n\n" "$go_rep" >> "$outdir"/func_log_"$date".txt

printf "  GO analysis files ("$outdir"/GO_analysis/):\n\n  go.obo\t\t\t\t->\tGO term information retrieved from Gene Ontology Consortium\n  GO_term_annotation\t\t\t->\tTab delim GO term info extracted from go.obo\n  all_go_ids.txt\t\t\t->\tList of unique GO IDs from each genome\n  GO_term_count_table\t\t\t->\tTable of GO term counts per genome\n  GO_term_genome_count_table.tab\t->\tGenome only GO term counts\n  GO_term_pancore_count_table.tab\t->\tPan-core genome only GO term counts\n  'GENOME'_go_annotated.tab\t\t->\tGO term counts, name & namespace per genome\n  'GENOME'_go_count.txt\t\t\t->\tNo. genes per GO term\n\n" >> "$outdir"/func_log_"$date".txt  

printf "  Genome top 50 GO term statistics\n  ----------------------------------\n\n" >> "$outdir"/func_log_"$date".txt

awk -F"\t" 'NR>1{for(i=7;i<=NF;i++){count[NR]+=$i;if($i>max[NR])max[NR]=$i;if($i<min[NR]||i==7)min[NR]=$i;val[i]=$i};for(x in val){diff[NR]+=($x-(count[NR]/(NF-3)))^2}; std[NR]=sqrt(diff[NR]/(NF-4));sem[NR]=(std[NR]/sqrt(NF-3)); print $1"\t"$3"\t"(count[NR])"\t"(count[NR]/(NF-3))"\t"(std[NR])"\t"(sem[NR])"\t"(min[NR])"\t"(max[NR])"\t"(max[NR]-min[NR])"\t"$2}' "$outdir"/GO_analysis/GO_term_genome_count_table.tab | sort -k3 -nr | head -n 50 | awk -F"\t" 'BEGIN{printf "%12s |%20s |%7s |%7s |%7s |%7s |%5s |%5s |%6s | %s\n  ----------------------------------------------------------------------------------------------------\n",  "Go ID", "Namespace" , "Total" , "Mean" ,"STD", "SEM" , "Min" , "Max", "Range" , "Name"}{printf "%12s |%20s |%7s |%7.2f |%7.2f |%7.2f |%5s |%5s |%6s | %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >> "$outdir"/func_log_"$date".txt


printf "\n\n\n Pancore genome top 50 GO Terms\n  ----------------------------------\n\n" >> "$outdir"/func_log_"$date".txt

tail -n +2 "$outdir"/GO_analysis/GO_term_pancore_count_table.tab | sort -t $'\t' -k5 -nr | head -n 51 | awk -F"\t" 'BEGIN{printf "%12s |%20s |%5s |%9s |%7s |%7s |%7s | %s\n  ----------------------------------------------------------------------------------------------------\n", "GO ID", "Namespace", "Core", "Softcore", "Shell", "Cloud", "Total", "Name"}{printf "%12s |%20s |%5s |%9s |%7s |%7s |%7s | %s\n", $1, $3, $5, $8, $7, $4, $6, $2}' >> "$outdir"/func_log_"$date".txt

printf "\n\n\n  End of Report\n\n\n" >> "$outdir"/func_log_"$date".txt

cat "$outdir"/func_log_"$date".txt >> "$outdir"/"$taxname"_analysis_report.txt


printf "\n\nAnalysis Complete\n\n"
exit




