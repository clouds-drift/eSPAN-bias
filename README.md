# eSPAN-bias
Bias calculation for eSPAN experiments

eSPAN (enrichment and sequencing of protein associated nascent DNA) method is able to analyze the enrichment of proteins of interest, including histones and its modifications, at replicating chromatin in a strand-specific manner in mammalian cells. Strand specific read coverage on plus and minus strands can be calculated from the clean bam file by [deeptools](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html). Then read coverage bias can be analyzed by 'eSPAN-bias' scripts.

## Requirements
* deeptools
* R packages: "openxlsx", "GenomicRanges", "Rsamtools", "rtracklayer", "reshape2", "zoo", "ggplot2","plyr"
* R packages needed for 'pptx' output: "officer","rvg"

## Usage
The eSPAN-bias scripts provides the analysis of bias calculation based on bigwig files:
* cal_bias.R: Calculate read coverage bias for each sample
* normalize_bias.R: Normalize the eSPAN bias against the BrdU IP bias to exclude the background effects in replication
* draw_bias.R: Draw plot for bias results and group samples by user’s specification

