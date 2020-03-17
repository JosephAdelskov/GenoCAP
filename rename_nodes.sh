#! /bin/bash

usage() { printf "\nrename_nodes.sh\n\nUSAGE:\t[-h] [-i TREE_FILE] [-n NAME_FILE] [-o OUT_FILE]\n\n"; 
   exit;}

while getopts "hi:o:n:" opt; do
   case $opt in
	  i)
		 infile=$OPTARG
		 ;;
	  o)
		 outfile=$OPTARG
		 ;;
	  n)
		 nfile=$OPTARG
		 ;;
	  h | *)
		 usage
		 ;;
   esac
done

if [ -z $infile ] || [ -z $outfile ] || [ -z $nfile ]; then
   usage
   exit
fi

mytree="$(cat "$infile")"


while read oldn newn; do

	newn="$(printf "$newn" | sed -e 's/ /_/g' -e 's/[,|:|(|)]//g' )"

	mytree="${mytree/$oldn)/"$newn"_"$oldn")}"
	mytree="${mytree/$oldn,/"$newn"_"$oldn",}"

done < $nfile

printf "$mytree" > "$outfile"


exit
