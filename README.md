# Long-read Phasing Pipeline
Jack Humphrey

Takes aligned BAM files from long-read RNA-seq (ONT or PacBio) and runs the LongCallR phasing tool.

Computes phasing rate (the proportion of mapped reads that can be phased into haplotypes)

Performs allele-specific expression, allele-specific splicing, and allele-specific editing (experimental)


## Dependencies:

python (>=3.8)
- argparse
- os
- sys
- pysam
- concurrent.futures

Samtools - full path to the executable should be in config

LongCallR - this must be installed and the full path to the folder should be in the config

Reference FASTA and GTF files - setup to use GRCh38 and GENCODE, using "chr1" type chromosome names

Assumes that each sample corresponds to an aligned long-read RNA-seq BAM file - use Minimap2 with the suggested settings for the library type.

## Future additions:

Collation of results across samples

Integration with GWAS summary statistics 

## Metadata

expects a two column TSV file with columns "sample" and "bam_path" - the absolute path to BAM files

## Running the example dataset

example/ contains examples of the config.yaml, the metadata.tsv, and a fragment of a PacBio long-read RNA-seq file (a few genes on chr21).

to test that the pipeline works, run:

```
snakemake -s Snakefile --configfile example/example_config.yaml -npr -c1
```
