# bedpeintersect
A script for multi-intersect and annotation of HiC bedpee files

The script intersect the genome contacts represented in a bedpe file to all the bedfiles stored
in a folder. Also annotate the closest gen each genomic region.

## Requirements

- bedtools
- pybedtools


## Files preparation.

### Annotation files.

Currently this repo only ships processed annotation files for human genome GRCh38

The Ensemble annotation .gtf file was processed as follows:

1. Add the 'chr' string to the chromosome names.
2. Separate the 'gene' features in one file and the 'exon' and 'utr' features in other file.
3. Sort the both files with `bedtools sort`.

### Bedpe input file.

The .bedpe input file must contain only the 6 fields correspondint to a genomic contact, i.e. chromosome name, start position and end position for both anchors of the contact.






