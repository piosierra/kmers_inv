# Workflow to find very large inversions with kmers.

```
# Obtain the kmers
FastK -k31 -t1 CM012070.1.fa 
Tabex CM012070.1.ktab LIST > kmers_all

FastK -k31 -t1 DA-386_06-X.fa 
Tabex DA-386_06-X.fa LIST > kmers_all_2

# Get only lines with freq = 1
awk -F" " '{if ($4 == 1) print $2}' kmers_all > unique_kmers
awk -F" " '{if ($4 == 1) print $2}' kmers_all_2 > unique_kmers2

# Get just lines common to both samples
comm -12 <( sort unique_kmers ) <( sort unique_kmers2 ) > unique_kmers_both


# Get just a sample
shuf -n 1000000 unique_kmers_both > unique_kmers_both_sample 

# Convert to fasta
awk '{print "> seq " ++count ORS $0}' unique_kmers_both_sample > unique_kmers_both_sample.fa

# Align
bwa mem -a CM012070.1.fa unique_kmers_both_sample.fa > kmers_aligned.sam

# Remove header
sed '/^@/d' < kmers_aligned.sam > ukmers.sam

# See unique values of first colum
awk '{ a[$1]++ } END { for (b in a) { print b } }' ukmers.sam

# Remove lines with more than one hit (next line is 256 or 272, then remove rest of 256/272)
awk 'NR > 1 && !(/\t272\t/ || /\t256\t/) { print p } { p = $0 } END { print}' ukmers.sam > ukmers_clean1.sam
sed '/^\t\(272\|256\|4\)/d' < ukmers_clean1.sam > ukmers_clean.sam

# Remove surplus fields
awk '{print $1 "\t" $3 "\t" $9}' ukmers_clean.sam | sort -k2 -n > ukmers_clean_table

# Get list for realignment
awk '{print $3}' ukmers_clean_table > unique_kmers_both_sample_2

# Convert to fasta
awk '{print "> seq " ++count ORS $0}' unique_kmers_both_sample_2 > unique_kmers_both_sample_2.fa

# Align
bwa mem -a DA-386_06-X.fa unique_kmers_both_sample_2.fa > kmers_aligned_2.sam

# Remove header
sed '/^@/d' < kmers_aligned_2.sam > ukmers_ali.sam

# See unique values of first colum
awk '{ a[$1]++ } END { for (b in a) { print b } }' ukmers_ali.sam

# Remove lines with more than one hit (next line is 256 or 272, then remove rest of 256/272)
awk 'NR > 1 && !(/\t272\t/ || /\t256\t/) { print p } { p = $0 } END { print}' ukmers_ali.sam > ukmers_ali_clean1.sam
sed '/^\t\(272\|256\|4\)/d' < ukmers_ali_clean1.sam > ukmers_ali_clean.sam

# Remove surplus fields
awk '{print $1 "\t" $3 "\t" $9}' ukmers_ali_clean.sam | sort -k2 -n > ukmers_ali_clean_table



Final:

./kinv.sh -1 ../assemblies/2RL/DA-416_04-2RL.fa -2 ../assemblies/2RL/DA-402_03-2RL.fa -o temp2

```
