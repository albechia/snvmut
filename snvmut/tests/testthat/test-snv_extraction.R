library(testthat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(VariantAnnotation)

genome <- BSgenome.Hsapiens.UCSC.hg38
path_to_vcf <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(path_to_vcf, "hg38")[1:15]

context("Tests for snv_extraction")



test_that("The context length value is a single odd integer number", {
  expect_error(snv_extraction(vcf, 2, genome), "context_length must be an odd integer number")
  expect_error(snv_extraction(vcf, -1, genome), "context_length must be an odd integer number")
  expect_error(snv_extraction(vcf, "three", genome), "The second parameter is the 'context_length' and must be a single numeric value")
})

test_that("Imput parameters are of the correct type", {
  expect_error(snv_extraction(vcf, 5, "hg38"),
  "The third parameter must be a BSgenome object")
})

test_that("The parameters provided can't be more than 3", {
  expect_error(snv_extraction(vcf, 3, genome, 6),
  "Only three parameters are required: 'VCF', 'context_length', 'reference genome'")
})
