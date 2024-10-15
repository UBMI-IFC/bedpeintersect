# bedpeintersect
A script for multi-intersect and annotation of HiC bedpee files

The script intersect the genome contacts represented in a bedpe file to all the bedfiles stored
in a folder. Also annotate the closest gen each genomic region.

## Requirements

### python
- pandas
- pybedtools

### External tools
- bedtools

### Conda environment.

The conda or mamba environment to run this script can be created as follows.

```bash
conda create -n bedpeintersect-env -c bioconda -c conda-forge pybedtools pandas
```

or 

```bash
mamba create -n bedpeintersect-env -c bioconda -c conda-forge pybedtools pandas
```


## Files preparation.

### Annotation files.

Currently this repo only ships processed annotation files for human genome GRCh38. This files path are 
hardcoded

The Ensemble annotation .gtf file was processed as follows:

1. Add the 'chr' string to the chromosome names.
2. Separate the 'gene' features in one file and the 'exon' and 'utr' features in other file.
3. Sort  both files with bedtools (`bedtools sort`).

### Bedpe input file.

The bedpe input file must contain only the 6 fields correspondint to a genomic contact/interaction, i.e.:

  - anchor 1 chromosome name 
  - anchor 1 start position 
  - anchor 1 end position
  - anchor 2 chromosome name 
  - anchor 2 start position 
  - anchor 2 end position


### Bed files folder.

All the bedfiles to intersect must be in one folder.

## Outputs.

The script generates 3 files:

1. "'melt_'": This is the melted version of the original bedpe. Thisi is, each anchor of each contact 
is in a separate line. Also an identifier name for each contact is created.
2. "'annotated_'": This is the oríginal bedbe with the name of each contact added and the results 
for each bed file.
3. "'table_'": Detailed results in table format, each line is an anchor. Contains the closest genomic
feature for each genomic region.



