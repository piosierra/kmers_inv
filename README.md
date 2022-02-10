# Identifying very large inversions in full genomes.

### 1- Identify large inversions on fasta files based on mapping of unique k-mers. Generates a dotplot of the alignments. 
Requires  FastK, bwa, R.
```
./kinv.sh \
        -1 [First fasta file] \
        -2 [Second fasta file] \
```

### 2- Detect inversions between two paths of a GFA file. Generates two bed files with the coordinates of the inversions for each path.  
Requires R.

```
./inv_from_paths \
        -g [path to gfa file] \
        -a [name of path 1] \
        -b [name of path 2] \
        -c [name of chr, used just to include the info on the bed file] \
        -p [:optional: min percentaje of inverted sequence on total inversion :def=0.8] \
        -s [:optional: min size of inversion :def=50:]
```