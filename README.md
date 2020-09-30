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

Some example data from GSE112522 are deposited in folder 'example_data' for testing. Annotation files are deposited in folder 'annotation',including:

saccer3_G1_origin.bed, yeast early origins

saccer3_G2_origin.bed, yeast late origins

saccer3_TSS.bed, yeast gene TSS

***
### 1. Rscript cal_bias.R -h
```
Options:
        -w WATSON, --watson=WATSON
                Bigwig file of plus strand. (Required)
                Multiple files can be joined by ',', like 'sample1_plus.bw, sample2_plus.bw'.
        -c CRICK, --crick=CRICK
                Bigwig file of minus strand. (Required)
                Multiple files can be joined by ',', like 'sample1_minus.bw, sample2_minus.bw'. 
        -r REFERENCE_POINT, --reference_point=REFERENCE_POINT
                A BED file of interested annotations and the middle point will be aligned. (Required)
                Default is the yeast G1 replication origins.
        -m METHOD, --method=METHOD
                Method for bias calculation, should be 'partition' or 'log'. (Optional, default='partition')
                'partition' = (w+c) / (w-c)
                'log' = log2 (w / c)
        --bw=BW
                Output folder for bigwig signals of genome-wide bias.It can be loaded by IGV to check individual locus. (Optional, default='bw_bias')
        -o MATRIX, --matrix=MATRIX
                Output folder for bias matrix. Will be used for bias normalization or bias plot. (Optional, default='bias_matrix')
        -s SUMMATFILE, --sumMatFile=SUMMATFILE
                Summary of output matrix file. (Optional, default='my_bias_mat.xlsx')
        --flank=FLANK
                Flanking regions around the reference point. It should be two numeric values joined by ','. (Optional, default='-2000,2000')
                Positive value means downstream and negative value means upstream to the reference point.
        -b BIN, --bin=BIN
                Bin size for matrix calculation. (Optional, default=100)
        --smooth=SMOOTH
                Number of bins on each side will be used for smooth. (Optional, default=0)
        -t THRESHOLD, --threshold=THRESHOLD
                Regions with coverage less than this value will be sikpped. (Optional, default=0)
        -N N_THREAD, --N_thread=N_THREAD
                Number of parallel processing. (Optional, default=2)
        -f FORCE, --force=FORCE
                If not forced, existing files will be skipped for calculation again. (Optional, default=T)
```
` Rscript cal_bias.R -w example_data/GSM3072025_y929_WT_H3K4me3_eSPAN/GSM3072025_y929_WT_H3K4me3_eSPAN_plus.bw -c example_data/GSM3072025_y929_WT_H3K4me3_eSPAN/GSM3072025_y929_WT_H3K4me3_eSPAN_minus.bw -r annotation/saccer3_G1_origin.bed`
##### Output:
bw_bias/GSM3072025_y929_WT_H3K4me3_eSPAN_partion.bw, genomwide bias signal.
bias_matrix/GSM3072025_y929_WT_H3K4me3_eSPAN.xlsx, bias matrix surrounding the G1 origins.
bias_matrix/GSM3072025_y929_WT_H3K4me3_eSPAN.txt, average bias profile surrounding the G1 origins.

##### Details:
* ‘--method’ can choose to use either partition or log2 ration method for bias calculation. Partition method is recommended, since it can scale the bias to (-1, 1).

| ![Partition](./graph/partition.jpg) | <img src="./graph/logRatio.pdf" width="300px"> |
| :--: | :--: |
| **(W+C)/(W-C)** | **log2(W/C)** |

* ‘--threshold’ can filter unreliable regions with low coverage by setting a proper coverage threshold. 

| aa  |
| :--: |
| **t** |

| First Header  | Second Header |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |

