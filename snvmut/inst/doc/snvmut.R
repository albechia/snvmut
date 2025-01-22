## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(snvmut)

## ----message=FALSE------------------------------------------------------------
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

## -----------------------------------------------------------------------------
snv_table <- snv_count(vector_snv)
snv_table

## ----message=FALSE------------------------------------------------------------
library(ggplot2)
library(patchwork)

snv_plots <- snv_graphics(snv_table)
snv_plots

## -----------------------------------------------------------------------------
sessionInfo()

