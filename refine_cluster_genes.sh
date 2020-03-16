#! /usr/bin/env bash

#: Title			: Refine Cluster Genes
#: Date				: 18/05/2018
#: Author			: Joseph Adelskov
#: Version			: 1.0
#: Description		: Continues on from the Get Genomes Phyla Amphora script. Given the name of the genome of interest, selects a group of highly similar genomes via marker gene blast proceeded by Average Nucleotide Identity analysis. Core/Pan genomic analysis is then conducted on the selection of genomes using get_homologues.


usage () { printf "\n\nREFINE_CLUSTER_GENES\n\n\tUSAGE:\t[-I INDIR] [-T THREADS] [-G GENOME]\n\nFrom the output directory from get_genomes_phyla_amphora.sh this script will select groups of closely related genomes to genome(s) of interest. Core/Pan genomic analysis will then be performed on this set of genomes using get_homologues.\n\n\t-i\toutput directory from get_genomes_phyla_amphora.sh\n\t-o\tName of output directory to save analysis inside 'indir'\n\t-p\tOutput file prefix to use, default is outdir\n\t-g\ttarget genome of interest or comma-separated list of genomes\n\t-t\tnumber of threads/cpus to use\n\t-b\tminimum blast identity score from marker genes to keep genome in analysis\n\n\n"

	exit;
}

################################################################################
## Direct user script input to variables
################################################################################

while getopts hi:p:t:g:b:o: option
do
	case "${option}"
		in
		i) indir=${OPTARG};;
		o) outdir=${OPTARG};;
		p) prefix=${OPTARG};;
		g) gtarget=${OPTARG};;
		t) threads=${OPTARG};;
		b) minblast=${OPTARG};;
		h | *)
			usage
			;;
	esac
done


# Ensure required variables have a value

if [ -z "$indir" ]; then

	echo "No Input directory specified, please provide the output directory of get_genomes_phyla_amphora.sh via the '-i' option"

	usage

fi

if [ -z "$outdir" ]; then

	outdir="$indir"

fi

if [ -z "$gtarget" ]; then

	echo "No central genome(s) specified, exiting"

	usage

fi

if [ -z "$threads" ]; then

	threads="1"

	echo "Using a single cpu thread for analysis"

fi

if [ -z "$taxname" ]; then

	printf "\nPlease provide the taxonomic name used in the last script via the -n option"
	taxname="$(ls "$indir" | awk -F"_" '/_acc_list.txt/ {print $1}')"

fi

	printf "\ntaxname set to "$taxname"\n"

if [ -z "$prefix" ]; then

	prefix="${outdir##*/}"

fi

printf "\nOutfile prefix set to "$prefix"\n"

if [ -z "$minblast" ]; then

	minblast="default"
	
	printf "\nMinimum blast identity score not set, genomes will be\nselected based on the average blast idenity score + or - 1 standard deviation\nif greater than  or less than 100 quality genomes respectively"

fi

# Determine if $gtarget is a comma-sep list, file or single variable

if [ -r "$gtarget" ]; then

	while read line; do

		gtlist+=($line)

	done < "$gtarget"


elif [[ $gtarget =~ "" ]]; then

	IFS=',' read -a gtlist <<< $gtarget

else

	gtlist[0]=$gtarget


fi

# get current time and date of executing script

date=$(date | sed 's/ /-/g')


################################################################################
## Check that script dependencies are available and working
################################################################################

printf "\n\n Checking software dependencies are installed and working\n"

# check blastp is installed and working
if [ -n "$(which blastp)" ] && [ -n "$(blastp -help)" ]; then

	printf "\nBLASTp is installed $(which blastp)\n"

else

	printf "\nBLASTp is not installed (correctly) !!\n"
	missdep="yes"

fi

# check average_nucleotide_identity.py (pyani) is installed and working

if [ -n "$(which average_nucleotide_identity.py)" ] && [ -n "$(average_nucleotide_identity.py -h)" ]; then

	printf "\naverage_nucleotide_identity.py is installed $(which average_nucleotide_identity.py)\n"

else

	printf "\naverage_nucleotide_identity.py is not installed (correctly) !!\n"
	missdep="yes"

fi

# check PROKKA is installed and working

if [ -n "$(which prokka)" ]; then

	printf "\nPROKKA is installed $(which prokka)\n"

else

	printf "\nPROKKA is not installed (correctly) !!\n"
	missdep="yes"

fi

# check get_homologues.pl and other scripts are installed

printf "\nChecking if Get homologues suit is installed\n"

if [ -n "$(which get_homologues.pl)" ]; then

	printf "\nget_homologues.pl is installed: $(which get_homologues.pl)\n"

else

	printf "\nget_homologues.pl is not installed (correctly) !!\n"
	missdep="yes"

fi

if [ -n "$(which compare_clusters.pl)" ]; then

	printf "\ncompare_clusters.pl is available: $(which compare_clusters.pl)\n"

else

	printf "\ncompare_clusters.pl is not available !!\n"
	missdep="yes"

fi

if [ -n "$(which parse_pangenome_matrix.pl)" ]; then

	printf "\nparse_pangenome_matrix.pl is available: $(which parse_pangenome_matrix.pl)\n"

else

	printf "\nparse_pangenome_matrix.pl is not available !!\n"
	missdep="yes"

fi

# check trimal and readal are installed

if [ -n "$(which trimal)" ]; then

	printf "\ntrimal is installed $(which trimal)\n"

else

	printf "\ntrimal is not installed (correctly) !!\n"
	missdep="yes"

fi

if [ -n "$(which readal)" ]; then

	printf "\nreadal is installed $(which readal)\n"

else

	printf "\nreadal is not installed (correctly) !!\n"
	missdep="yes"

fi

# check RAXML is instaled and working

if [ -n "$(which raxmlHPC-PTHREADS)" ]; then

	printf "\nRAXML is installed $(which raxmlHPC-PTHREADS)\n"

else

	printf "\nRAXML is not installed (correctly) !!\n"
	missdep="yes"

fi

if [ "$missdep" == "yes" ]; then

	printf "\nOne or more software dependencies cannot be called correctly by this workflow script, Exiting\n\n"

	exit 1

else

	printf "\nRequired dependencies have been found, starting workflow\n"

fi


################################################################################
## Start Refine_log.txt contributions from script
################################################################################

mkdir -p "$outdir"

printf "\n\n  ================================================================================\n\n  Running script Refine_cluster_genes.sh | $date  \n\n\n  Script running parameters\n  --------------------------------------------------------------------------------\n  \tInput directory\t\t=\t"$indir"\n  \tTaxonomy name\t\t=\t"$taxname"\n  \tNumber of threads\t=\t"$threads"\n\n\n" >> "$outdir"/refine_log_"$date".txt

## Add information on BLAST genome selection to analysis log

printf "  Genome selection criteria\n  --------------------------------------------------------------------------------\n  The first task of this script is to select genomes with sufficient similarity to\n  the user specified set of target genomes. To do this BLAST analysis will be used\n  to compare the marker genes of the specified genomes against the genomes retrieved\n  and filtered in the previous script \"Get_genomes_phyla_amphora.sh\".\n\n  \tTarget Genomes\n  \t--------------\n" >> "$outdir"/refine_log_"$date".txt
 

## List Target genomes in analysis log

for i in ${!gtlist[@]}; do

	let "count = $i + 1"

	printf "  \t"$count".\t"${gtlist[$i]}"\n" >> "$outdir"/refine_log_"$date".txt

done

printf "\n\n  The minimum BLAST identity value can be specified by the user. However if a value\n  is not specified the default setting will be used.\n\n  Default: mean(+/- 1 std) {Greater/less than 100 quality genomes}\n\n   \tMinimum BLAST identity\t=\t"$minblast"\n\n" >> "$outdir"/refine_log_"$date".txt

################################################################################
## Marker BLAST genome selection
################################################################################

# Make new directory for the Marker gene blast analysis

mkdir -p "$outdir"/marker_blast

# check if checkm or phyla amphora was used to get marker alignment

if [ -f "$indir"/phyla_amphora_analysis/"$taxname"_comb_aln.fasta ]; then

# Remove the gaps (-)s from the combined alignment and save the result in the marker_blast folder

sed 's/-//g' "$indir"/phyla_amphora_analysis/"$taxname"_comb_aln.fasta > "$outdir"/marker_blast/"$taxname"_markers_comb.fasta

elif [ -f "$indir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta ]; then

sed 's/-//g' "$indir"/checkm_marker_analysis/"$taxname"_comb_aln.fasta > "$outdir"/marker_blast/"$taxname"_markers_comb.fasta

else

	printf "\nNO marker alignment file found, exiting!\n\n"

	exit

fi


# Build blast database from the marker gene alignments produced from phyla amphora

if [[ ! $(ls -A "$outdir"/marker_blast/marker_comb* 2>/dev/null) ]]; then

	makeblastdb -dbtype prot -in "$outdir"/marker_blast/"$taxname"_markers_comb.fasta -title marker_comb -out "$outdir"/marker_blast/marker_comb

fi

# Extract the target genome combined alignments from the total to use as the blast queries which will save time

if [ -f "$outdir"/marker_blast/"$prefix"_genome_queries.fasta ]; then

	rm "$outdir"/marker_blast/"$prefix"_genome_queries.fasta

fi

let "iteration = 0"

for i in ${!gtlist[@]}; do

	let "iteration++"

	echo "input genome $iteration = ${gtlist[$i]}"

	awk -v genome="${gtlist[$i]}" 'BEGIN{RS=">"} $1 ~ genome {print ">"$0}' "$outdir"/marker_blast/"$taxname"_markers_comb.fasta >> "$outdir"/marker_blast/"$prefix"_genome_queries.fasta

done

# Select closely related genomes to target genomes

	# Run blastp total alignments against query alignments

if [ ! -f "$outdir"/marker_blast/"$prefix"_blast_out.txt ]; then

	printf "\nRunning BLASTP on genome queries\n"

	blastp -db "$outdir"/marker_blast/marker_comb -query "$outdir"/marker_blast/"$prefix"_genome_queries.fasta -max_hsps 1 -outfmt '6 qseqid sseqid pident evalue qlen length' -out "$outdir"/marker_blast/"$prefix"_blast_out.txt -num_threads "$threads"

fi

# From the blast results select the more similar matches >= 85% identity, atleast 80% total alignment


#if [[ ! $(ls -A "$outdir"/selected_genomes/ 2>/dev/null) ]]; then

	if [ -n "$genome_list" ]; then

		unset -v genome_list

	fi


	qual_gen_cnt="$(ls "$indir"/quality_genomes/ | wc -l)"

	printf "\n\n%20s%10s\n  -------------------------------\n" "Target" "Genomes"
	for i in ${!gtlist[@]}; do

		# Determine mean and standard deviation of blast identity using awk redirect it through while read
		while read mean std; do

			# print important variables to terminal (for error checking)
			#printf "\n${gtlist[$i]} $mean $std\n"

			if [ "$minblast" != "default" ]; then

				genome_list+=($(awk -v gen="${gtlist[$i]}" -v minb="$minblast" '$1==gen && $3 >= minb {print $2}' "$outdir"/marker_blast/"$prefix"_blast_out.txt | sort | uniq ))
				
			elif [ "$qual_gen_cnt" -ge "100" ]; then


				genome_list+=($(awk -v gen="${gtlist[$i]}" -v mean="$mean" -v std="$std" 'BEGIN{cutoff=mean+std}$1 ~ gen && $3 >= cutoff {print $2}' "$outdir"/marker_blast/"$prefix"_blast_out.txt | sort | uniq ))

			else
			

				genome_list+=($(awk -v gen="${gtlist[$i]}" -v mean="$mean" -v std="$std" 'BEGIN{cutoff=mean-std}$1 ~ gen && $3 >= cutoff {print $2}' "$outdir"/marker_blast/"$prefix"_blast_out.txt | sort | uniq ))
			# clean up array by removing redundant values

			fi

			genome_list=($( printf "%s\n" "${genome_list[@]}" | sort | uniq))

		done < <(awk -v gen="${gtlist[$i]}" '$1 ~ gen && $2 !~ gen {print}' "$outdir"/marker_blast/"$prefix"_blast_out.txt | awk '{sum+=$3; sumsq+=$3*$3; nbr=NR} END {print sum/nbr" "sqrt(sumsq/nbr-(sum/nbr)^2)}')


			target_count[$i]="${#genome_list[@]}"

			printf "%20s%10d\n" "${gtlist[$i]}" "${#genome_list[@]}"

	done


	## Add BLAST results to analysis log

	printf "  BLAST selection results:\n\n  \tQuery Genome\tTotal Genomes Selected\n  \t--------------------------------------\n" >> "$outdir"/refine_log_"$date".txt

	for i in ${!gtlist[@]}; do

		printf "  \t%12s\t            %10d\n" "${gtlist[$i]}"  "${target_count[$i]}" >> "$outdir"/refine_log_"$date".txt

	done


	# Copy the selected genome files to a new directory called selected_genomes

	mkdir -p "$outdir"/selected_genomes
	

	for i in ${!genome_list[@]}; do

		cp "$indir"/quality_genomes/"${genome_list[$i]}".fna "$outdir"/selected_genomes/

	done

#fi

## Add contig statistics to analysis log

printf "\n\n\n  \tSelected genome contig statistics\n  \t-----------------------------------------------\n" >> "$outdir"/refine_log_"$date".txt

grep -c ">" "$outdir"/selected_genomes/*.fna | awk -F":" '{sum = sum + $2; if ($2 > max) max = $2; if(min>$2||NR==1)min = $2}END{print "  \tNo. genomes\t\t\t=\t"NR"\n  \tTotal contigs\t\t\t=\t"sum"\n  \tMean contigs\t\t\t=\t"sum/NR"\n  \tMin contigs\t\t\t=\t"min"\n  \tMax contigs\t\t\t=\t"max"\n  \tRange\t\t\t\t=\t"max-min"\n\n"}' >> "$outdir"/refine_log_"$date".txt


################################################################################
## Extract genome information for select genomes only
################################################################################

mkdir -p "$outdir"/genome_information

head -n 1 "$indir"/genome_information/"$taxname"_genome_char.tab > "$outdir"/genome_information/"$prefix"_selected_genome_char.tab
head -n 1 "$indir"/genome_information/"$taxname"_genome_info.tab > "$outdir"/genome_information/"$prefix"_selected_genome_info.tab

for file in "$outdir"/selected_genomes/*.fna; do

	filename="${file##*/}"
	gen_name="${filename/.fna/}"

	grep "^$gen_name" "$indir"/genome_information/"$taxname"_genome_char.tab >> "$outdir"/genome_information/"$prefix"_selected_genome_char.tab
	grep "^$gen_name" "$indir"/genome_information/"$taxname"_genome_info.tab >> "$outdir"/genome_information/"$prefix"_selected_genome_info.tab

	grep "^$gen_name" "$indir"/genome_information/"$taxname"_name_key.txt >> "$outdir"/genome_information/"$prefix"_name_key.txt
	grep "^$gen_name" "$indir"/genome_information/"$taxname"_tax_naming.txt >> "$outdir"/genome_information/"$prefix"_tax_naming.txt

done


################################################################################
## Determine average nucleotide identity (ANI) of the selected group of genomes
################################################################################

if [ -d "$outdir"/ANI_out ]; then

	if [[ ! $(ls -A "$outdir"/ANI_out) ]]; then

		rmdir "$outdir"/ANI_out

		average_nucleotide_identity.py -i "$outdir"/selected_genomes -o "$outdir"/ANI_out -m ANIb -s 100 -g --workers "$threads" || { echo "Cannot run ANI" ; exit 1; }

	fi

else

	printf "\nCalculating average nucleotide identity\n"

	average_nucleotide_identity.py -i "$outdir"/selected_genomes -o "$outdir"/ANI_out -m ANIb -s 100 -g --workers "$threads" || { echo "Cannot run ANI" ; exit 1; }
	
fi

## Combined ANI pID, alignment length and alignment coverage into a single file using awk

if [ -f "$outdir"/ANI_out/ANIb_percentage_identity.tab ] && [ ! -f "$outdir"/"$prefix"_ANI_long_format.tab ]; then

printf "Genome 1\tGenome 2\tANIb\tCoverage\tAlignment length\n" > "$outdir"/"$prefix"_ANI_long_format.tab

while read gen1 gen2 ani; do

	coverage=$(awk -v gen1="$gen1" -v gen2="$gen2" 'NR==1 { for(i=1; i<=NF; i++) header[i]=$i ; next } $1 ~ gen1 {for(x=1;x<=NF;x++) {if(header[x] == gen2) print ($(x+1)*100)}}' "$outdir"/ANI_out/ANIb_alignment_coverage.tab )

	alength=$(awk -v gen1="$gen1" -v gen2="$gen2" 'NR==1 { for(i=1; i<=NF; i++) header[i]=$i ; next } $1 ~ gen1 {for(x=1;x<=NF;x++) {if(header[x] == gen2) print $(x+1)}}' "$outdir"/ANI_out/ANIb_alignment_lengths.tab)
	
	printf "%s\t%s\t%s\t%s\t%s\n" $gen1 $gen2 $ani $coverage $alength >> "$outdir"/"$prefix"_ANI_long_format.tab

done < <(awk 'NR==1{for(i=1; i<=NF; i++) header[i]=$i ; next} {for(x=2; x<=NF; x++) {if ( $x > 0.80 ) print $1"\t"header[(x-1)]"\t"($x*100)} next}' "$outdir"/ANI_out/ANIb_percentage_identity.tab)

fi

## Add ANI info to analysis log

printf "\n  Average Nucleotide Identity analysis (ANI)\n  --------------------------------------------------------------------------------\n  ANI analysis was performed on the whole set of selected genomes using the python\n  pyani package. The analysis was conducted using ANIb (BLAST comparison) method\n  with a fragment size of 100 bp.\n\n  ANI results stored in dir: "$outdir"/ANI_out\n\n  The results including ANI, coverage and alignment length have been converted\n  to tabular format and stored in the file:\n\n  "$outdir"/"$prefix"_ANI_long_format.tab\n\n\n" >> "$outdir"/refine_log_"$date".txt


#####################################################################################
## Annotate genomes with PROKKA to produce genbank format files for get_homologues
#####################################################################################

mkdir -p "$outdir"/prokka_annotation

gen_count="$(ls "$outdir"/selected_genomes | wc -l )"

if [ "$(ls "$outdir"/prokka_annotation/*.gbk 2>/dev/null | wc -l)" -lt "$gen_count" ]; then

# if the array genome_list does not exist due to running the script again after ANI analysis was performed then reassign the array using the selected_genomes directory

	if [[ ! "${genome_list[@]}" ]]; then

		if [ -d "$outdir"/selected_genomes ]; then

			if [[ $( ls -A "$outdir"/selected_genomes 2>/dev/null) ]]; then

				for file in "$outdir"/selected_genomes/*.fna; do

					genome_list+=($(printf "${file##*/}" | sed 's/.fna//'))

				done

			else

				printf "\ndirectory selected_genomes is empty, exiting\n" && exit 1

			fi

		else

			printf "\nDirectory selected_genomes does not exist, exiting\n" && exit 1

		fi

	fi

	printf "\nRunning PROKKA annotations of selected genomes\n"

	for i in ${!genome_list[@]}; do

		# read in values from $prefix_tax_nomenclature.txt and set as values for the prokka command

		while read accession genus species strain; do

			prokaccession="$accession"
			prokgenus="$genus"
			prokspecies="$species"
			prokstrain="$strain"

			if [ -z "$prokaccession" ]; then

				prokaccession="${genome_list[$i]}"
				prokgenus="$(grep "${genome_list[$i]}" "$outdir"/genome_information/"$prefix"_tax_nomenclature.txt | awk '{print $2}' )"
				prokspecies="$(grep "${genome_list[$i]}" "$outdir"/genome_information/"$prefix"_tax_nomenclature.txt | awk '{print $3}' )"

			fi

			if [ -z "$prokstrain" ]; then

				prokstrain="$prokaccession"

			fi

		done < <(awk -F"=" -v gen="${genome_list[$i]}" '$1 ~ gen {print $1}' "$outdir"/genome_information/"$prefix"_tax_naming.txt)

		prokka --outdir "$outdir"/prokka_annotation  --prefix "$prokaccession" --genus "$prokgenus" --species "$prokspecies" --strain "$prokstrain" --locustag "$prokaccession" --compliant --force --cpu "$threads" "$outdir"/selected_genomes/${genome_list[$i]}.fna

	done

fi

# copy the .gbk files to a new directory for get_homologues clustering

mkdir -p "$outdir"/get_homologues_analysis/gh_gbk

if [[ ! $(ls -A "$indir"/gh_gbk/ 2>/dev/null) ]]; then
	
	cp "$outdir"/prokka_annotation/*.gbk "$outdir"/get_homologues_analysis/gh_gbk/

fi

## Add Prokka annotation info to analysis log

printf "  PROKKA annotation\n  --------------------------------------------------------------------------------\n  To prepare for Get_homologues each of the selected genomes was analysed and\n  annotated using PROKKA prokaryotic genome annotation tool.  \n\n  \tProkka annotation statistics\n  \t-----------------------------------------" >> "$outdir"/refine_log_"$date".txt

grep ">" -c "$outdir"/prokka_annotation/*.ffn | awk -F":" '{sum = sum + $2; if($2 > max) max = $2; if($2<min||NR==1)min = $2}END{print "\n  \tGenomes\t\t\t=\t"NR"\n  \tTotal genes\t\t=\t"sum"\n  \tMean\t\t\t=\t"sum/NR"\n  \tMin\t\t\t=\t"min"\n  \tMax\t\t\t=\t"max"\n  \tRange\t\t\t=\t"max-min}' >> "$outdir"/refine_log_"$date".txt


################################################################################
## Perform pan-core genome analysis with get_homologues
################################################################################

# make directory to store all get_homologues analysis

mkdir -p "$outdir"/get_homologues_analysis

if [ ! -d "$outdir"/get_homologues_analysis/gh_gbk_homologues ]; then

	printf "\nRunning Get Homologues\n"

	get_homologues.pl -d "$outdir"/get_homologues_analysis/gh_gbk -n "$threads" -t 0 -M -c -z -A || { echo "Cannot run get_homologues.pl" ; exit 1; }

	# move the get_homologues output into the in/output directory

	mv gh_gbk_homologues "$outdir"/get_homologues_analysis/

fi

# Rename the output from get_homologues before running the compare clusters script

if [ ! -d "$outdir"/get_homologues_analysis/gh_gbk_homologues/gh_algOMCL_out ]; then

	printf "\nRenaming output directory\n"

	mv "$outdir"/get_homologues_analysis/gh_gbk_homologues/*_algOMCL_e0_ "$outdir"/get_homologues_analysis/gh_gbk_homologues/gh_algOMCL_out

	# When renaming the cluster output directory, the .cluster_list file should also be given the same prefix

	mv "$outdir"/get_homologues_analysis/gh_gbk_homologues/*.cluster_list "$outdir"/get_homologues_analysis/gh_gbk_homologues/gh_algOMCL_out.cluster_list

fi

# use get_homologues to cluster intergenic regions

if [ ! -d "$outdir"/get_homologues_analysis/gh_gbk_intergenic ]; then
	
	printf "\nPerforming Intergenic sequence clustersing\n"

	get_homologues.pl -d "$outdir"/get_homologues_analysis/gh_gbk -n "$threads" -g  || { echo "Cannot run get_homologues.pl" ; exit 1; }

	mv gh_gbk_homologues "$outdir"/get_homologues_analysis/gh_gbk_intergenic

fi

## Testing run get_homologues.pl with Pfam domain prediction enabled (takes alot more computer time)

#if [[ ! $(ls -A "$indir"/gh_gbk_pfam 2>/dev/null ) ]]; then

#	get_homologues.pl -d "$indir"/gh_gbk -n "$threads" -c || { echo "Cannot run get_homologues.pl with pfam domain prediction" ; exit 1; }

#	mv gh_gbk_homologues "$indir"/gh_gbk_pfam

#fi


# Run compare_clusters.pl on gene and intergenic clusters produced via get_homologues

mkdir -p "$outdir"/cluster_analysis

if [[ ! "$(ls -A "$outdir"/cluster_analysis 2>/dev/null)" ]]; then

	printf "\n\nRunning compare_clusters.pl for protein gene homologue clusters\n\n"

	compare_clusters.pl -o "$outdir"/cluster_analysis/prot_clst_intersection -m -T -d "$outdir"/get_homologues_analysis/gh_gbk_homologues/gh_algOMCL_out

	printf "\n\nRunning compare_clusters.pl for nucleotide gene homologue clusters\n\n"

	compare_clusters.pl -o "$outdir"/cluster_analysis/nucl_clst_intersection -m -n -T -d "$outdir"/get_homologues_analysis/gh_gbk_homologues/gh_algOMCL_out

	printf "\n\nRunning compare_clusters.pl for syntenic (protein) gene homologue clusters\n\n"

	compare_clusters.pl -o "$outdir"/cluster_analysis/synt_clst_intersection -m -s -T -d "$outdir"/get_homologues_analysis/gh_gbk_homologues/gh_algOMCL_out

fi

## Analyse Pan-Core Genome ; Parse pangenome matrix

mkdir -p "$outdir"/cluster_analysis/pancore_genome/gh_prot
mkdir -p "$outdir"/cluster_analysis/pancore_genome/gh_nucl
mkdir -p "$outdir"/cluster_analysis/pancore_genome/gh_synt

if [[ ! "$(ls -A "$outdir"/cluster_analysis/pancore_genome/gh_prot/ 2>/dev/null)" ]]; then

	parse_pangenome_matrix.pl -m "$outdir"/cluster_analysis/prot_clst_intersection/pangenome_matrix_t0.tab -s

	cp "$outdir"/cluster_analysis/prot_clst_intersection/pangenome_matrix_t0* "$outdir"/cluster_analysis/pancore_genome/gh_prot/

fi

if [[ ! "$(ls -A "$outdir"/cluster_analysis/pancore_genome/gh_nucl/ 2>/dev/null)" ]]; then

	parse_pangenome_matrix.pl -m "$outdir"/cluster_analysis/nucl_clst_intersection/pangenome_matrix_t0.tab -s

	cp "$outdir"/cluster_analysis/nucl_clst_intersection/pangenome_matrix_t0* "$outdir"/cluster_analysis/pancore_genome/gh_nucl/

fi

if [[ ! "$(ls -A "$outdir"/cluster_analysis/pancore_genome/gh_synt/ 2>/dev/null)" ]]; then

	parse_pangenome_matrix.pl -m "$outdir"/cluster_analysis/synt_clst_intersection/pangenome_matrix_t0_s.tab -s

	cp "$outdir"/cluster_analysis/synt_clst_intersection/pangenome_matrix_t0_s* "$outdir"/cluster_analysis/pancore_genome/gh_synt/

fi

## Add info on cluster analysis to analysis log

printf "\n\n\n  Get_homologues Pan/core genome analysis\n  --------------------------------------------------------------------------------\n  The genbank annotation files produced by PROKKA analysis were used as input for\n  'get_homologues.pl' to cluster homologous genes into clusters. Using genbank\n  (.gbk) files as input also allowed for homologous intergenic regions between\n  genes to be clustered together. Following homologue clustering the\n  compare_clusters.pl script was then use to produce the pangenomic matrix and\n  Phylip parsimony analysis for either nucl, prot or syntenic(prot) clusters.\n  Sytenic clusters are defined by get_homologues as clusters where all homologous\n  genes in the cluster have atleast one similar neibouring gene.\n\n" >> "$outdir"/refine_log_"$date".txt

printf "  Using the pangenome matricies for each set of clusters enome compositional\n  analysis was then performed via the parse_pangenome_matrix.pl script.\n  The pangenome composition consists of core (1 gene per genome), softcore\n  (gene in 95% of genomes), shell (shared between 2 or more genomes), and\n  cloud (strain specific genes).\n\n" >> "$outdir"/refine_log_"$date".txt

# print out stats for each set of clusters

printf "  Standard Homologue Clusters\n  ---------------------------\n\n  \tCluster stats\n  \t------------------------------------------------\n" >> "$outdir"/refine_log_"$date".txt

find "$outdir"/cluster_analysis/prot_clst_intersection/ -type f -name '*.faa' | xargs -P 8 -I {} grep -c "^>" {} | awk '{sum += $1;if($1 > max)max=$1;if($1<min||NR==1)min=$1}END{print "  \tHomologue clusters\t\t=\t"NR"\n  \tTotal genes\t\t\t=\t"sum"\n  \tMaximum genes per cluster\t=\t"max"\n  \tMinimum genes per cluster\t=\t"min"\n  \tMean genes per cluster\t\t=\t"sum/NR}' >> "$outdir"/refine_log_"$date".txt

printf "\n  \tPangenome composition\n  \t------------------------------------------------\n  \tCore\t\t\t\t=\t%d\n" "$(grep ".faa" -c  "$outdir"/cluster_analysis/pancore_genome/gh_prot/pangenome_matrix_t0__core_list.txt)" >> "$outdir"/refine_log_"$date".txt

printf "  \tSoftcore\t\t\t=\t%d\n" "$(grep ".faa" -c  "$outdir"/cluster_analysis/pancore_genome/gh_prot/pangenome_matrix_t0__softcore_list.txt)" >> "$outdir"/refine_log_"$date".txt

printf "  \tShell\t\t\t\t=\t%d\n" "$(grep ".faa" -c  "$outdir"/cluster_analysis/pancore_genome/gh_prot/pangenome_matrix_t0__shell_list.txt)" >> "$outdir"/refine_log_"$date".txt

printf "  \tCloud\t\t\t\t=\t%d\n" "$(grep ".faa" -c  "$outdir"/cluster_analysis/pancore_genome/gh_prot/pangenome_matrix_t0__cloud_list.txt)" >> "$outdir"/refine_log_"$date".txt

printf "\n\n  Syntenic Homologue Clusters\n  ---------------------------\n\n  \tCluster stats\n  \t------------------------------------------------\n" >> "$outdir"/refine_log_"$date".txt

find "$outdir"/cluster_analysis/synt_clst_intersection/ -type f -name '*.faa' | xargs -P 8 -I {} grep -c "^>" {} | awk '{sum += $1;if($1 > max)max=$1;if($1<min||NR==1)min=$1}END{print "  \tHomologue clusters\t\t=\t"NR"\n  \tTotal genes\t\t\t=\t"sum"\n  \tMaximum genes per cluster\t=\t"max"\n  \tMinimum genes per cluster\t=\t"min"\n  \tMean genes per cluster\t\t=\t"sum/NR}' >> "$outdir"/refine_log_"$date".txt


printf "\n  \tPangenome composition\n  \t------------------------------------------------\n  \tCore\t\t\t\t=\t%d\n" "$(grep ".faa" -c  "$outdir"/cluster_analysis/pancore_genome/gh_synt/pangenome_matrix_t0_s__core_list.txt)" >> "$outdir"/refine_log_"$date".txt

printf "  \tSoftcore\t\t\t=\t%d\n" "$(grep ".faa" -c  "$outdir"/cluster_analysis/pancore_genome/gh_synt/pangenome_matrix_t0_s__softcore_list.txt)" >> "$outdir"/refine_log_"$date".txt

printf "  \tShell\t\t\t\t=\t%d\n" "$(grep ".faa" -c  "$outdir"/cluster_analysis/pancore_genome/gh_synt/pangenome_matrix_t0_s__shell_list.txt)" >> "$outdir"/refine_log_"$date".txt

printf "  \tCloud\t\t\t\t=\t%d\n" "$(grep ".faa" -c  "$outdir"/cluster_analysis/pancore_genome/gh_synt/pangenome_matrix_t0_s__cloud_list.txt)" >> "$outdir"/refine_log_"$date".txt


################################################################################
## Rename tips of pangenome parsimony trees produced by compare_clusters.pl
################################################################################

mkdir -p "$outdir"/cluster_analysis/pangenome_pars_trees

if [[ ! "$(ls -A "$outdir"/cluster_analysis/pangenome_pars_trees 2>/dev/null)" ]]; then


	printf "\nRenaming nodes of pangenomic parsimony trees produced from compare clusters\n"

	# load trees into variables

	prot_pg_pars_tree="$(sed -e 's/_/./g' -e 's/GCA./GCA_/g' "$outdir"/cluster_analysis/prot_clst_intersection/pangenome_matrix_t0.phylip.ph | sed 's/.gbk...//g' | awk '{print $0";"}')"

	nucl_pg_pars_tree="$(sed -e 's/_/./g' -e 's/GCA./GCA_/g' "$outdir"/cluster_analysis/nucl_clst_intersection/pangenome_matrix_t0.phylip.ph | sed 's/.gbk...//g' | awk '{print $0";"}')"

	synt_pg_pars_tree="$(sed -e 's/_/./g' -e 's/GCA./GCA_/g' "$outdir"/cluster_analysis/synt_clst_intersection/pangenome_matrix_t0_s.phylip.ph | sed 's/.gbk...//g' | awk '{print $0";"}')"

	# rename tree nodes via variable substituion using taxonomy nomenclature file

	while read acc nomen; do

		nomen="$(printf "$nomen" | sed 's/ /_/g')"

		prot_pg_pars_tree="${prot_pg_pars_tree/"$acc"/"$nomen"_"$acc"}"
	
		nucl_pg_pars_tree="${nucl_pg_pars_tree/"$acc"/"$nomen"_"$acc"}"

		synt_pg_pars_tree="${synt_pg_pars_tree/"$acc"/"$nomen"_"$acc"}"

	done < "$outdir"/genome_information/"$prefix"_tax_naming.txt

	# remove any characters that are invalid in newick format

	prot_pg_pars_tree="$(printf "$prot_pg_pars_tree" | sed  's/\(\[\|\]\|{\|}\||\)//g')" 

	nucl_pg_pars_tree="$(printf "$nucl_pg_pars_tree" | sed  's/\(\[\|\]\|{\|}\||\)//g')" 

	synt_pg_pars_tree="$(printf "$synt_pg_pars_tree" | sed  's/\(\[\|\]\|{\|}\||\)//g')" 

	printf "%s\n" "$prot_pg_pars_tree" > "$outdir"/cluster_analysis/pangenome_pars_trees/"$prefix"_prot_pangenome_parsimony_tree.nwk

	printf "%s\n" "$nucl_pg_pars_tree" > "$outdir"/cluster_analysis/pangenome_pars_trees/"$prefix"_nucl_pangenome_parsimony_tree.nwk
	
	printf "%s\n" "$synt_pg_pars_tree" > "$outdir"/cluster_analysis/pangenome_pars_trees/"$prefix"_synt_pangenome_parsimony_tree.nwk

fi

## Add info for pangenome trees to analysis log

printf "\n\n\n  The pangenome parsimony trees produced using the compare_clusters.pl script\n  processed to rename the phylogenetic tips to full genome names/identifiers\n  and have been stored in the directory:\n\n  "$outdir"/cluster_analysis/pangenome_pars_trees/\n\n" >> "$outdir"/refine_log_"$date".txt

################################################################################
## Build phylogenomic tree from core genes/proteins
################################################################################

mkdir -p "$outdir"/cluster_analysis/core_genome/gh_prot_omcl
mkdir -p "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl

# get the core genome gene clusters using the core core_list.txt produced earlier
# use awk commands to clean up gene files and remove duplicate genes to the same genome in the files

if [[ ! $(ls -A "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/ 2>/dev/null) ]]; then

	printf "\nCleaning up core genes (proteins) for phylogenomic tree analysis\n"

	let "core_prot_inc=0"

	while read core_gene; do

		let "core_prot_inc++"

		awk -F"|" '/^>/ {gsub("ID:Prokka:",""); gsub("_[0-9]{5} "," "); print $1; next}1' "$outdir"/cluster_analysis/prot_clst_intersection/"$core_gene" > "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"$core_gene"_tmp

		while read num acc; do

			if [ "$num" -gt "1" ]; then

				awk -v acc="$acc" 'BEGIN{RS=">"; y = "0"; x = "0"}$1 ~ acc{x = length($NF); {if(x > y) y = x}}$1 ~ acc{if(length($NF) == y && z == "0") z = "1"; gene=$1"\n"$NF} END{print ">"gene}' "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"$core_gene"_tmp >> "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"$core_gene"

			else

				awk -v acc="$acc" 'BEGIN{RS=">"}$1 ~ acc{print ">"$1"\n"$NF}' "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"$core_gene"_tmp >> "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"$core_gene"

			fi

		done < <(awk -F"|" '/^>/ {gsub(">",""); print $1}' "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"$core_gene"_tmp | sort | uniq -c)

		rm "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"$core_gene"_tmp

		printf "\nAligning core gene $core_prot_inc ($core_gene) using muscle"

		muscle -quiet -in "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"$core_gene" -out "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/"${core_gene/.faa/_aligned.faa}"

	
	done < "$outdir"/cluster_analysis/pancore_genome/gh_prot/pangenome_matrix_t0__core_list.txt

	printf "\nConcatenating aligned genes (proteins) into single alignment\n"

	# use catfasta2phyml.pl to concatenate the aligned core genes into a single alignment
	catfasta2phyml.pl -f "$outdir"/cluster_analysis/core_genome/gh_prot_omcl/*_aligned.faa > "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aligned.faa

	# use trimal to remove poorly represented regions of the alignment

	trimal -in "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aligned.faa -out "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aln_trimed.faa -gt 0.8 -st 0.001

	readal -in "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aln_trimed.fna -info
	readal -in "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aln_trimed.fna -info > "$outdir"/cluster_analysis/core_genome/gh_prot_core_alignment_info.txt
fi

if [[ ! $(ls -A "$outdir"/cluster_analysis/core_genome/prot_raxml/ 2>/dev/null) ]] && [ -f "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aln_trimed.faa ]; then

	mkdir -p "$outdir"/cluster_analysis/core_genome/prot_raxml

	printf "\nBuilding Phylogenomic tree from core genes (proteins) using RAxML\n"

	raxmlHPC-PTHREADS -f a -s "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aln_trimed.faa -x 127 -# 100 -m PROTGAMMAJTT -n "$prefix"_core_prot_100b -p 129 -T "$threads"

	mv RAxML* "$outdir"/cluster_analysis/core_genome/prot_raxml/

fi


################################################################################
## Build phylogenomic tree from core genes/nucleotide
################################################################################

if [[ ! $(ls -A "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/ 2>/dev/null) ]]; then

	printf "\nCleaning up core genes (nucleotides) for phylogenomic tree analysis\n"

	let "nucl_core_inc=0"

	while read core_gene; do

		let "nucl_core_inc++"

		awk -F"|" '/^>/ {gsub("ID:Prokka:",""); gsub("_[0-9]{5} "," "); print $1; next}1' "$outdir"/cluster_analysis/nucl_clst_intersection/"$core_gene" > "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"$core_gene"_tmp

		# this while loop removes multiple/duplicate genes from the same genome
		# the longest gene will be kept while others removed
		while read num acc; do

			if [ "$num" -gt "1" ]; then

				awk -v acc="$acc" 'BEGIN{RS=">"; y = "0"; z = "0"}$1 ~ acc{x = length($NF); {if(x > y) y = x}}$1 ~ acc{if(length($NF) == y && z == "0") z = "1"; gene=$1"\n"$NF} END{print ">"gene}' "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"$core_gene"_tmp >> "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"$core_gene"

			else

				awk -v acc="$acc" 'BEGIN{RS=">"}$1 ~acc{print ">"$1"\n"$NF}' "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"$core_gene"_tmp >> "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"$core_gene"

			fi

		done < <(awk -F"|" '/^>/ {gsub(">",""); print $1}' "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"$core_gene"_tmp | sort | uniq -c)

		rm "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"$core_gene"_tmp

		printf "\nAligning core gene $nucl_core_inc ($core_gene) using muscle"

		muscle -quiet -in "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"$core_gene" -out "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/"${core_gene/.fna/_aligned.fna}"

	done < "$outdir"/cluster_analysis/pancore_genome/gh_nucl/pangenome_matrix_t0__core_list.txt

	printf "\nConcatenating aligned genes (nucleotides) into single alignment\n"

	catfasta2phyml.pl -f "$outdir"/cluster_analysis/core_genome/gh_nucl_omcl/*_aligned.fna > "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aligned.fna

	# use trimal to remove poorly represented regions of the alignment

	trimal -in "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aligned.fna -out "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aln_trimed.fna -gt 0.8 -st 0.001

	readal -in "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aln_trimed.fna -info
	readal -in "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aln_trimed.fna -info > "$outdir"/cluster_analysis/core_genome/gh_nucl_core_alignment_info.txt
fi

if [[ ! $(ls -A "$outdir"/cluster_analysis/core_genome/nucl_raxml/ 2>/dev/null) ]] && [ -f "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aln_trimed.fna ]; then

	mkdir -p "$outdir"/cluster_analysis/core_genome/nucl_raxml

	printf "\nBuilding Phylogenomic tree from core genes (nucleotide) using RAxML\n"

	raxmlHPC-PTHREADS -f a -s "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aln_trimed.fna -x 127 -# 100 -m GTRGAMMA -n "$prefix"_core_nucl_100b -p 129 -T "$threads"

	mv RAxML* "$outdir"/cluster_analysis/core_genome/nucl_raxml/

fi

################################################################################
## Build phylogenomic tree from syntenic core genes (prot)
################################################################################

mkdir -p "$outdir"/cluster_analysis/core_genome/gh_synt_omcl

# get the core genome gene clusters using the core core_list.txt produced earlier
# use awk commands to clean up gene files and remove duplicate genes to the same genome in the files

if [[ ! $(ls -A "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/ 2>/dev/null) ]]; then

	printf "\nCleaning up core syntenic genes (proteins) for phylogenomic tree analysis\n"

	let "core_synt_inc=0"

	while read core_gene; do

		let "core_synt_inc++"

		awk -F"|" '/^>/ {gsub("ID:Prokka:",""); gsub("_[0-9]{5} "," "); print $1; next}1' "$outdir"/cluster_analysis/synt_clst_intersection/"$core_gene" > "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"$core_gene"_tmp

		while read num acc; do

			if [ "$num" -gt "1" ]; then

				awk -v acc="$acc" 'BEGIN{RS=">"; y = "0"; x = "0"}$1 ~ acc{x = length($NF); {if(x > y) y = x}}$1 ~ acc{if(length($NF) == y && z == "0") z = "1"; gene=$1"\n"$NF} END{print ">"gene}' "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"$core_gene"_tmp >> "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"$core_gene"

			else

				awk -v acc="$acc" 'BEGIN{RS=">"}$1 ~ acc{print ">"$1"\n"$NF}' "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"$core_gene"_tmp >> "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"$core_gene"

			fi

		done < <(awk -F"|" '/^>/ {gsub(">",""); print $1}' "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"$core_gene"_tmp | sort | uniq -c)

		rm "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"$core_gene"_tmp

		printf "\nAligning core gene $core_synt_inc ($core_gene) using muscle"

		muscle -quiet -in "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"$core_gene" -out "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/"${core_gene/.faa/_aligned.faa}"

	
	done < "$outdir"/cluster_analysis/pancore_genome/gh_synt/pangenome_matrix_t0_s__core_list.txt

	printf "\nConcatenating aligned syntenic homologs (proteins) into a single alignment\n"

	# use catfasta2phyml.pl to concatenate the aligned core genes into a single alignment
	catfasta2phyml.pl -f "$outdir"/cluster_analysis/core_genome/gh_synt_omcl/*_aligned.faa > "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aligned.faa

	# use trimal to remove poorly represented regions of the alignment

	trimal -in "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aligned.faa -out "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aln_trimed.faa -gt 0.8 -st 0.001

	readal -in "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aln_trimed.fna -info
	readal -in "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aln_trimed.fna -info > "$outdir"/cluster_analysis/core_genome/gh_synt_core_alignment_info.txt
fi

if [[ ! $(ls -A "$outdir"/cluster_analysis/core_genome/synt_raxml/ 2>/dev/null) ]] && [ -f "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aln_trimed.faa ]; then

	mkdir -p "$outdir"/cluster_analysis/core_genome/synt_raxml

	printf "\nBuilding Phylogenomic tree from core synt genes (proteins) using RAxML\n"

	raxmlHPC-PTHREADS -f a -s "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aln_trimed.faa -x 127 -# 100 -m PROTGAMMAJTT -n "$prefix"_core_synt_100b -p 129 -T "$threads"

	mv RAxML* "$outdir"/cluster_analysis/core_genome/synt_raxml/

fi
################################################################################
## Rename nodes of core genome phylogenomic trees
################################################################################

#if [[ ! $(ls -A "$outdir"/cluster_analysis/core_genome/core_genome_trees/ 2>/dev/null) ]]; then

mkdir -p "$outdir"/cluster_analysis/core_genome/core_genome_trees

if [ ! -f "$outdir"/cluster_analysis/core_genome/core_genome_trees/"$prefix"_prot_core_genome.nwk ]; then

	if [ -f "$outdir"/cluster_analysis/core_genome/prot_raxml/RAxML_bipartitions."$prefix"_core_prot_100b ]; then

		core_prot_tree="$(cat "$outdir"/cluster_analysis/core_genome/prot_raxml/RAxML_bipartitions."$prefix"_core_prot_100b)"

		while read accn newn; do

			newn="$(printf "$newn" | sed 's/\((\|)\|{\|}\|\[\|\]\|:\|;\|,\)//g' | sed 's/ /_/g')"

			core_prot_tree="${core_prot_tree/"$accn"/"$newn"_"$accn"}"


		done < "$outdir"/genome_information/"$prefix"_tax_naming.txt

		printf "%s\n" $core_prot_tree > "$outdir"/cluster_analysis/core_genome/core_genome_trees/"$prefix"_prot_core_genome.nwk

	fi

fi

if [ ! -f "$outdir"/cluster_analysis/core_genome/core_genome_trees/"$prefix"_nucl_core_genome.nwk ]; then

	if [ -f "$outdir"/cluster_analysis/core_genome/nucl_raxml/RAxML_bipartitions."$prefix"_core_nucl_100b ]; then

		core_nucl_tree="$(cat "$outdir"/cluster_analysis/core_genome/nucl_raxml/RAxML_bipartitions."$prefix"_core_nucl_100b)"

		while read accn newn; do

			newn="$(printf "$newn" | sed 's/\((\|)\|{\|}\|\[\|\]\|:\|;\|,\)//g' | sed 's/ /_/g')"

			core_nucl_tree="${core_nucl_tree/"$accn"/"$newn"_"$accn"}"

		done < "$outdir"/genome_information/"$prefix"_tax_naming.txt

		printf "%s\n" $core_nucl_tree > "$outdir"/cluster_analysis/core_genome/core_genome_trees/"$prefix"_nucl_core_genome.nwk
	
	fi

fi

if [ ! -f "$outdir"/cluster_analysis/core_genome/core_genome_trees/"$prefix"_synt_core_genome.nwk ]; then

	if [ -f "$outdir"/cluster_analysis/core_genome/synt_raxml/RAxML_bipartitions."$prefix"_core_synt_100b ]; then

		core_synt_tree="$(cat "$outdir"/cluster_analysis/core_genome/synt_raxml/RAxML_bipartitions."$prefix"_core_synt_100b)"

		while read accn newn; do

			newn="$(printf "$newn" | sed 's/\((\|)\|{\|}\|\[\|\]\|:\|;\|,\)//g' | sed 's/ /_/g')"

			core_synt_tree="${core_synt_tree/"$accn"/"$newn"_"$accn"}"


		done < "$outdir"/genome_information/"$prefix"_tax_naming.txt

		printf "%s\n" $core_synt_tree > "$outdir"/cluster_analysis/core_genome/core_genome_trees/"$prefix"_synt_core_genome.nwk

	fi

fi

################################################################################
## Finish additions to analysis log for core genome analysis
################################################################################

printf "\n  Core genome phylogenomic analysis\n  --------------------------------------------------------------------------------\n  The clusters identified as part of the core genome from the standard set of\n  homologue clusters were used to perform a phylogenomic analysis of the selected set\n  of genomes. Before the core genome homologues could be aligned duplicate genes\n  (from the same genome) had to be remove and gene headers/identifiers ammended\n  to prevent subsequent errors during alignment/analysis. The core genome homologues\n  were then aligned using MUSCLE and alignements were subsequently concatenated\n  into a single alignment using catfasta2phyml.pl. The concatenated alignment\n  was then processed with trimal to remove positions with less than 80%% representation.\n  Maximum likelihood phylogenetic analysis was then performed using the trimmed\n  alignment using RAXML and the tip labels were renamed using the\n  "$prefix"_tax_naming.txt file. The finished core genome ML tree files were\n  then stored in the directory:\n\n  "$outdir"/cluster_analysis/core_genome/core_genome_trees\n\n\n  Core genome protein homologue alignment\n  ---------------------------------------\n\n" >> "$outdir"/refine_log_"$date".txt


awk '/^>/{i++}!/>/{if(i>1) exit; print}' "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aln_trimed.faa | awk '{seqlen += length($0)}END{print "  \tLength\t= "seqlen" AA\n"}' >> "$outdir"/refine_log_"$date".txt

awk 'BEGIN{RS=">";ORS=""}{$1 =""; gsub(" ",""); print $0}' "$outdir"/cluster_analysis/core_genome/gh_prot_comb_aln_trimed.faa | sed -e 's/-//g' -e 's/./&\n/g' | sort | uniq -c | awk 'BEGIN{print "  \tComposition\n"}{AA[NR] = $2; count[NR] = $1}{sum += $1}END{for (i in AA) {printf "        %s  =  %.2f %\n", (AA[i]), (count[i]/sum*100)}; print "\n"}' >> "$outdir"/refine_log_"$date".txt

printf "  The parameters used in the phylogenetic analysis are listed below:\n\n  \tAlignment\t\t=\tAmino acid\n  \tMethod\t\t\t=\tMaximum Likelihood\n  \tEvolutionary model\t=\tGamma distritubion\n  \tSubstitution matrix\t=\tJones, Taylor & Thornton (JTT)\n  \tTest of phylogeny\t=\tBootstrap\n  \tReplicates\t\t=\t100\n\n  \tFinished renamed tree: "$prefix"_prot_core_genome.nwk\n\n\n" >> "$outdir"/refine_log_"$date".txt

printf "  Core genome nucleotide homologue alignment\n  ------------------------------------------\n\n" >> "$outdir"/refine_log_"$date".txt

awk '/^>/{i++}!/>/{if(i>1) exit; print}' "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aln_trimed.fna | awk '{seqlen += length($0)}END{print "  \tLength\t= "seqlen" bp\n"}' >> "$outdir"/refine_log_"$date".txt

awk 'BEGIN{RS=">";ORS=""}{$1 =""; gsub(" ",""); print $0}' "$outdir"/cluster_analysis/core_genome/gh_nucl_comb_aln_trimed.fna | sed -e 's/-//g' -e 's/./&\n/g' | sort | uniq -c | awk 'BEGIN{print "  \tComposition\n"}{nucl[NR] = $2; count[NR] = $1}{sum += $1}END{for (i in nucl){printf "        %s  =  %.2f %\n", (nucl[i]), (count[i]/sum*100)}; print "\n"}' >> "$outdir"/refine_log_"$date".txt

printf "  The parameters used in the phylogenetic analysis are listed below:\n\n  \tAlignment\t\t=\tNucleotide\n  \tMethod\t\t\t=\tMaximum Likelihood\n  \tEvolutionary model\t=\tGamma distritubion\n  \tSubstitution model\t=\tGTR\n  \tTest of phylogeny\t=\tBootstrap\n  \tReplicates\t\t=\t100\n\n  \tFinished renamed tree: "$prefix"_nucl_core_genome.nwk\n\n\n" >> "$outdir"/refine_log_"$date".txt

printf "  Core genome syntenic homologue (prot) alignment\n  -----------------------------------------------\n\n" >> "$outdir"/refine_log_"$date".txt

awk '/^>/{i++}!/>/{if(i>1) exit; print}' "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aln_trimed.faa | awk '{seqlen += length($0)}END{print "  \tLength\t=\t"seqlen" AA\n"}' >> "$outdir"/refine_log_"$date".txt

awk 'BEGIN{RS=">";ORS=""}{$1 =""; gsub(" ",""); print $0}' "$outdir"/cluster_analysis/core_genome/gh_synt_comb_aln_trimed.faa | sed -e 's/-//g' -e 's/./&\n/g' | sort | uniq -c | awk 'BEGIN{print "  \tComposition\n"}{AA[NR] = $2; count[NR] = $1}{sum += $1}END{for (i in AA) {printf "        %s  =  %.2f %\n", (AA[i]), (count[i]/sum*100)}; print "\n\n\n"}' >> "$outdir"/refine_log_"$date".txt

printf "  The parameters used in the phylogenetic analysis are\n  listed below:\n\n  \tAlignment\t\t=\tAmino acid\n  \tMethod\t\t\t=\tMaximum Likelihood\n  \tEvolutionary model\t=\tGamma distritubion\n  \tSubstitution matrix\t=\tJones, Taylor & Thornton (JTT)\n  \tTest of phylogeny\t=\tBootstrap\n  \tReplicates\t\t=\t100\n\n  \tFinished renamed tree: "$prefix"_synt_core_genome.nwk\n\n\n" >> "$outdir"/refine_log_"$date".txt

# copy analysis report from indir then add refine_cluster results

cp "$indir"/"$taxname"_analysis_report.txt "$outdir"/"$prefix"_analysis_report.txt

cat "$outdir"/refine_log_"$date".txt >> "$outdir"/"$prefix"_analysis_report.txt

printf "\nrefine_cluster_genes.sh Complete\n\n"

exit








