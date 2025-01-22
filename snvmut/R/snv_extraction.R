#' The snv_extraction function in the snvmut package takes as input a VCF file, a reference genome
#' and a length number. The function returns a vector of SNV mutations that is of length equal to the
#' input length parameter.
#'
#'
#' @param vcf_file a VCF object
#' @param context_length a single number which must be an odd integer
#' @param ref_genome a BSgenome object
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import VariantAnnotation
#' @import Biostrings
#' @importFrom MatrixGenerics rowRanges
#' @import GenomicRanges
#' @import IRanges
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(VariantAnnotation)
#' genome_ref <- BSgenome.Hsapiens.UCSC.hg38
#' path_to_vcf <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(path_to_vcf, "hg38")
#' extracted_snv <- snv_extraction(vcf, 5, genome_ref)
#' @usage extracted_snv <- snv_extraction(vcf_file, context_length, ref_genome)
#'
#'
#' @export

snv_extraction <- function(vcf_file, context_length, ref_genome, ...) {

  #checks!
  #check for wrong parameter number (if they are too many > error)
  extra_para <- list(...)
  if (length(extra_para) > 0) {
    stop("Only three parameters are required: 'VCF', 'context_length', 'reference genome'.")
  }

  #check if the parameters are of the correct type and in the correct order
  if (!inherits(vcf_file, "VCF")) {
    stop("The first parameter must be a VCF object.")
  }

  if (!(is.numeric(context_length) && length(context_length) == 1)) {
    stop("The second parameter is the 'context_length' and must be a single numeric value.")
  }

  if (!(context_length == as.integer(context_length) &&
        context_length %% 2 != 0 && context_length > 0)) {
    stop("context_length must be an odd integer number.")
  }

  if (!inherits(ref_genome, "BSgenome")) {
    stop("The third parameter must be a BSgenome object.")
  }

  #using the function rowRanges to filter the vcf file for snvs
  snv_loc <- rowRanges(vcf_file[start(vcf_file) == end(vcf_file) & width(unlist(alt(vcf_file))) == 1])

  #empty vector to store mutation information
  num_snv_loc <- length(snv_loc)
  list <- vector("list", num_snv_loc)

  #each snv gets processed
  for (i in 1:num_snv_loc) {
    #standardize chromosome names to start with 'chr'
    chr <- ifelse(startsWith(as.character(seqnames(snv_loc[i])), "chr"),
                  as.character(seqnames(snv_loc[i])),
                  paste0("chr", as.character(seqnames(snv_loc[i]))))

    #extract snv position
    position <- start(snv_loc[i])
    #a range gets created surrounding the snv using the "context_lenght" parameter
    snv_range <- GRanges(chr, IRanges(start = position - context_length, end = position + context_length))
    #retrieve the sequence from the reference genome for the given range
    refSeq <- getSeq(ref_genome, snv_range)
    #extract the reference and alternate alleles
    ref_allele <- as.character(snv_loc[i]$REF)
    alt_allele <- as.character(unlist(snv_loc[i]$ALT))

    #function that takes care of situations in which the reference allele is A or G
    adjust_sequence <- function(seq, ref, alt) {
      #convert them to reverse complement
      if (ref %in% c("A", "G")) {
        reversedSeq <- as.character(reverseComplement(seq))
        reversedRef <- as.character(reverseComplement(DNAString(ref)))
        reversedAlt <- as.character(reverseComplement(DNAString(alt)))
        return(list(reversedSeq, reversedRef, reversedAlt))
      } else {
        return(list(as.character(seq), ref, alt))
      }
    }

    values <- adjust_sequence(refSeq, ref_allele, alt_allele)

    #extract also the sequences around the snv
    before <- substr(values[[1]], 1, context_length)
    after <- substr(values[[1]], length(values[[1]]), context_length)

    #put together the full string for each snv
    list[[i]] <- paste0(before, "[", values[[2]], ">", values[[3]], "]", after)
  }

  #make the list a vector for easier visualization
  vector_snv <- unlist(list, use.names = FALSE)

  #fix the lenght of each element in the vector according to the "context_lenght" parameter
  short_vector_snv <- vector("character", length(vector_snv))
  for (i in seq_along(vector_snv)) {
    snv_string <- vector_snv[i]

    #find the position of the first '['
    start_pos <- regexpr("\\[", snv_string)[[1]]
    halfcont <- (context_length - 1) %/% 2

    # Extract the substring around the SNV
    short_string <- substr(snv_string,
                              start = max(1, start_pos - halfcont),
                              stop = min(nchar(snv_string), start_pos + halfcont + 4)) # +4 because it stops just after the second ']'
    short_vector_snv[i] <- short_string
  }

  return(short_vector_snv)
}
