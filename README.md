---
title: "snvmut: an R package for SNV extraction"
author: Chiara Albertini
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{my-vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(snvmut)
```

### Package introduction

The package snvmut is a Bioconductor compliant R package that allows the user to search for SNV (Single Nucleotide Variant) inside of a .vcf file, to store the mutations in a vector and to graphically plot the mutation counts.

The package has three main functions.

#### snv_extraction

This function returns a vector of SNVs and requires three parameters:

-   a VCF object

-   the reference genome as a a BSgenome object

-   the "context_lenght" peramter, which is a numerical value that corresponds to the final window of nucleotides (mutation included) that will be stored inside the vector. The mutation is in the format [REF\>ALT] and in the window are included also the upstream and downstream nucleotides as defined by the "context_lenght" parameter.

The function is optimized to take into account and solve the redundancy caused by the fact that, for example, C[G\>A]A is the same as "T[C\>T]G" on the reverse strand. This is solved by having all mutations report C or T as the REF base.

#### snv_count

This function takes as input the vector obtained with the snv_extraction function and returns a data frame with three columns:

-   "SNV" reports the mutations found in the vector

-   "Count" reports how many time the mutation was found in the vector

-   "Percentage" report the percentage of every mutation on the total number of mutations of the vector

#### snv_graphics

This function takes as input the data frame obtained with the snv_count function and returns two ggplot2 plots:

-   a bar plot

-   a pie chart

The plots allow for a quick and easy visualization of the mutation types and proportion over the total number of mutations.

### Function usage

Loading VCF files is extremely quick when using the VariantAnnotation package, as shown here. The reference genome can also be loaded using the BSgenome packaging and by specifying the preferred genome.

```{r message=FALSE}
library(BSgenome.Hsapiens.UCSC.hg38)
library(VariantAnnotation)
library(Biostrings)
library(GenomicRanges)
library(snvmut)

genome_ref <- BSgenome.Hsapiens.UCSC.hg38
path_to_vcf <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(path_to_vcf, "hg38")[1:55]

vector_snv <- snv_extraction(vcf, 5, genome_ref)
vector_snv
```

The obtained vector can then be processed using the snv_count function to obtain a data frame with the mutation counts and percentages.

```{r}
snv_table <- snv_count(vector_snv)
snv_table
```

The obtained data frame can then be processed using the snv_graphics function that assigns a color to each mutation type and plots a bar plot (displaying the mutation counts) and a pie chart (useful for a quick proportion analysis).

```{r message=FALSE}
library(ggplot2)
library(patchwork)

snv_plots <- snv_graphics(snv_table)
snv_plots
```

#### sessionInfo

```{r}
sessionInfo()
```
