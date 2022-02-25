#!/usr/bin/env bash

# Stop at any error
set -eo pipefail

cwd=$(pwd)

echo 
SECONDS=0
g1=false
g2=false
output_folder=ptemp
sample_size=100000
show_help=false 

# read the options
cmd=$0" "$@
TEMP=`getopt -o 1:2:s:o:h: -l genome1:,genome2:,sample-size:,output-folder:help: -n 'kinv' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -1|--genome1) g1=$2 ; shift 2 ;;
        -2|--genome2) g2=$2 ; shift 2 ;;
        -o|--output-folder) output_folder=$2 ; shift 2 ;;
        -s|--sample-size) sample_size=$2 ; shift 2 ;;
        -h|--help) show_help=true ; shift ;;
        --) shift ; break ;;
        *) echo "$2" "Internal error!" ; exit 1 ;;
    esac
done


if [[
     $g1 == false
    || $g2 == false
   ]];
then
    show_help=true
    >&2 echo "Mandatory arguments -1, -2"
fi



if [ $show_help == true ];
then
    echo "options:"
    echo "    -1, --genome1            file    First genome to compare"
    echo "    -2, --genome2            file    Second genome to compare"
    echo "    -s, --sample-size        number  # of kmers to sample"
    echo "    -o, --output-folder      dir     Output folder"
    exit
fi

mkdir -p $output_folder
mkdir -p temppp

gs1="${g1%.f*}"  # Removes .fa
gs2="${g2%.f*}"
name1="${gs1##*/}"  # Removes path until filename
name2="${gs2##*/}"

echo
echo "==> Genome1: $name1"
echo "==> Genome2: $name2"
echo "==> Sample size: $sample_size"
echo "==> Output folder: $output_folder"
echo
echo "=================================================="
echo
cd $output_folder
echo $(pwd)
if ! [[ -f "$name1"_all ]]; then
        echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Extracting kmers for $name1"
        FastK -k31 -t1 -Ptemppp $g1
        Symmex "$gs1".ktab "$gs1"_2.ktab
        Tabex "$gs1"_2.ktab LIST > "$cwd"/"$output_folder"/"$name1"_all

else
    echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Found kmers for $name1"
fi

if ! [[ -f "$name2"_all ]]; then
        echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Extracting kmers  for $name2"
        FastK -k31 -t1 -Ptemppp $g2
        Symmex "$gs2".ktab "$gs2"_2.ktab
        Tabex "$gs2"_2.ktab LIST > "$cwd"/"$output_folder"/"$name2"_all
else
    echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Found kmers for $name2"
fi



echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Getting kmers with freq = 1"
if ! [[ -f "$name1"_all_unique ]]; then
awk -F" " '{if ($4 == 1) print $2}' "$name1"_all > "$name1"_all_unique
fi
if ! [[ -f "$name2"_all_unique ]]; then
awk -F" " '{if ($4 == 1) print $2}' "$name2"_all > "$name2"_all_unique
fi

echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Getting kmers common to both samples"
if ! [[ -f unique_kmers_"$name1"_"$name2" ]]; then
    comm -12 "$name1"_all_unique "$name2"_all_unique > unique_kmers_"$name1"_"$name2"
fi

echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Sampling $sample kmers"
shuf -n $sample_size unique_kmers_"$name1"_"$name2" > unique_kmers_"$name1"_"$name2"_sample 

echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Converting to fasta"
awk '{print "> seq " ++count ORS $0}' unique_kmers_"$name1"_"$name2"_sample > unique_kmers_"$name1"_"$name2"_sample.fa

cd $cwd
if ! [[ -f "$g1".bwt ]]; then
        echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Creating index for $name1"
        bwa index $g1
fi


echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Aligning kmers to $name1"
bwa mem -a $g1 "$output_folder"/unique_kmers_"$name1"_"$name2"_sample.fa > "$output_folder"/unique_kmers_"$name1"_"$name2"_sample_aligned.sam

cd $output_folder

echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Removing header"
sed '/^@/d' < unique_kmers_"$name1"_"$name2"_sample_aligned.sam > unique_kmers_"$name1"_"$name2"_sample_aligned_headless.sam

echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Unique values of first colum"
awk '{ a[$1]++ } END { for (b in a) { print b } }' unique_kmers_"$name1"_"$name2"_sample_aligned_headless.sam



echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Remove sequences with more than one hit"
awk 'NR > 1 && !(/\t272\t/ || /\t256\t/) { print p } { p = $0 } END { print}' \
    unique_kmers_"$name1"_"$name2"_sample_aligned_headless.sam \
    > unique_kmers_"$name1"_"$name2"_sample_aligned_headless_clean_t.sam

sed '/^\t\(272\|256\|4\)/d' < unique_kmers_"$name1"_"$name2"_sample_aligned_headless_clean_t.sam \
    > unique_kmers_"$name1"_"$name2"_sample_aligned_headless_clean.sam


echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==>  Remove surplus fields"
awk '{print $1 "\t" $3 "\t" $9}' unique_kmers_"$name1"_"$name2"_sample_aligned_headless_clean.sam \
    | sort -k2 -n > unique_kmers_"$name1"_"$name2"_sample_aligned_headless_clean_table


echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==>  Get list for realignment"
awk '{print $3}' unique_kmers_"$name1"_"$name2"_sample_aligned_headless_clean_table \
    > unique_kmers_"$name1"_"$name2"_both_sample_2


echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==>  Convert to fasta"
awk '{print "> seq " ++count ORS $0}' unique_kmers_"$name1"_"$name2"_both_sample_2 \
    > unique_kmers_"$name1"_"$name2"_both_sample_2.fa

cd $cwd
if ! [[ -f "$g2".bwt ]]; then
        echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==> Creating index for $name2"
        bwa index $g2
fi


echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==>  Align kmers to $name2"
bwa mem -a $g2 "$output_folder"/unique_kmers_"$name1"_"$name2"_both_sample_2.fa \
    > "$output_folder"/unique_kmers_"$name1"_"$name2"_both_sample_aligned_2.sam

cd $output_folder

echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==>  Remove header"
sed '/^@/d' < unique_kmers_"$name1"_"$name2"_both_sample_aligned_2.sam \
    > unique_kmers_"$name1"_"$name2"_both_sample_aligned_2_headless.sam


echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==>  Unique values of first colum"
awk '{ a[$1]++ } END { for (b in a) { print b } }' unique_kmers_"$name1"_"$name2"_both_sample_aligned_2_headless.sam


echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==>  Remove sequences with more than one hit"
awk 'NR > 1 && !(/\t272\t/ || /\t256\t/) { print p } { p = $0 } END { print}' \
    unique_kmers_"$name1"_"$name2"_both_sample_aligned_2_headless.sam \
    > unique_kmers_"$name1"_"$name2"_both_sample_aligned_2_clean_t.sam
sed '/^\t\(272\|256\|4\)/d' < unique_kmers_"$name1"_"$name2"_both_sample_aligned_2_clean_t.sam \
    > unique_kmers_"$name1"_"$name2"_final_clean.sam


echo -e "$(($SECONDS / 3600)):$(($SECONDS / 60)):$(($SECONDS % 60))\t ==>  Remove surplus fields"
awk '{print $1 "\t" $3 "\t" $9}' unique_kmers_"$name1"_"$name2"_final_clean.sam \
    | sort -k2 -n > unique_kmers_"$name1"_"$name2"_final_clean_table

# Plot the results
Rscript ../inv_viz.R unique_kmers_"$name1"_"$name2"_final_clean_table \
    unique_kmers_"$name1"_"$name2"_sample_aligned_headless_clean_table \
    "$name1" \
    "$name2"

mv *.png $cwd

cd $cwd

rm -r temppp
rm -r $output_folder

