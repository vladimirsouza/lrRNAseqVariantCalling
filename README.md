# lrRNAseqVariantCalling


This repository contains all the code used for the manuscript *Transformation of alignment files improves the performance of variant callers for long-read RNA sequencing data* (DOI: XXX). The order in which the scripts were run is indicated by numbers in their names or folders.

In the manuscript, we present a pipeline to increase the performance of the variant callers `DeepVariant` and `Clair3` on Iso-Seq data. This pipeline consists of using `minimp2` (to align reads to a reference genome), `GATK`'s `SplitNCigarReads` function (to split reads at intronic regions), `flagCorrection` (a tool developed by us to manipulate BAM files output by `SplitNCigarReads` to make them adequate for deep learning-based variant callers), and a variant caller (we suggest `DeepVariant` or `Clair3`). The generic code for this pipeline is showed below.

## Tools required to be installed

Our pipeline required the following tools to be installed:
* `minimap2` (>= 2.17-r941), [see installation](https://github.com/lh3/minimap2);
* `samtools` (>= 1.9), [see installation](https://github.com/samtools/samtools);
* `GATK` (>= v4.1.9.0), [see installation](https://github.com/broadinstitute/gatk/releases);
* `R` (>= 4.0.5), [see installation](https://www.r-project.org/);
* a variant caller, we recomend `DeepVariant` (>= 1.1.0), [see installation](https://github.com/google/deepvariant/blob/r1.3/docs/deepvariant-quick-start.md) &mdash; we run it with [Singularity](https://github.com/apptainer/singularity/blob/master/INSTALL.md), but other option may work,
* or, alternatively, `Clair3` (>= v0.1-r5) [see installation](https://github.com/HKU-BAL/Clair3) &mdash; we run it with [Bioconda](https://bioconda.github.io/user/install.html), but other options may work; and
* the `R` packages:
  * Rsamtools (>= 2.4.0),
  * foreach (>= 1.5.0), and
  * doParallel (>= 1.0.15).

To install the `R` packages, run `R` and enter the following code:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsamtools")

install.packages(c("foreach", "doParallel"))
```






## to Run the pipeline

A FASTQ file with the Iso-Seq data is required. If your data is in a PacBio BAM format, you may want to use `bamToFastq` to convert it to FASTQ.

```
bamToFastq \
  -i ${INPUT_BAM} \
  -fq ${OUTPUT_FASTQ}
```

