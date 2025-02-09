#' The snv_graphics provides a graphical rapresentation of the data table generated by the function
#' snv_counts. It plots both a bar plot and a pie chart to quickly visualize the distribution and
#' proportion of the mutations.
#'
#' @param table_snv_counts a data frame like the one returned by snv_count
#' @import ggplot2
#' @import patchwork
#' @usage snv_plots <- snv_graphics(snv_table)
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(VariantAnnotation)
#' genome_ref <- BSgenome.Hsapiens.UCSC.hg38
#' path_to_vcf <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(path_to_vcf, "hg38")
#' vector_snv <- snv_extraction(vcf, 5, genome_ref)
#' snv_table <- snv_count(vector_snv)
#' snv_plots <- snv_graphics(snv_table)
#' @usage snv_plots <- snv_graphics(snv_table)
#'
#' @export


snv_graphics <- function(table_snv_counts) {

  #checks!
  #check to see if the data frame has the correct structure
  if (!("data.frame" %in% class(table_snv_counts)) ||
       !all(c("SNV", "Count", "Percentage") %in% colnames(table_snv_counts))) {
    stop("The data frame must have columns named 'SNV', 'Count', and 'Percentage'.")
  }

  #check to see if the elements in every column are of the correct type
  if (!all(grepl("^[TCGA]>[TCGA]$", table_snv_counts$SNV))) {
    stop("Each entry in the 'SNV' column must be in the format N>N, where 'N' is a nucleotide T, C, G, or A.")
  }

  if (!all(table_snv_counts$Count > 0 & table_snv_counts$Count == as.integer(table_snv_counts$Count))) {
    stop("The 'Count' column must contain positive integer numbers.")
  }

  if (!all(table_snv_counts$Percentage > 0)) {
    stop("The 'Percentage' column must contain positive numbers.")
  }


  #barplot
  p1 <- ggplot(table_snv_counts, aes(x = SNV, y = Count, fill = SNV)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "SNV type counts", x = "SNV type", y = "Count") +
    theme(legend.position = "none")

  #pie chart
  p2 <- ggplot(table_snv_counts, aes(x = "", y = Percentage, fill = SNV, label = paste0(Percentage, "%"))) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = paste0(Percentage, "%")), position = position_stack(vjust = 0.5), color = "white", size = 2) +
    theme_void() +
    labs(title = "Proportion of SNV types") +
    theme(legend.position = "none")

  #display the two plots in a single image
  p1 + p2 + plot_layout(ncol = 2)

}
