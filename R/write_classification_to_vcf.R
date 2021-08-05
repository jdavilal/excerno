#' Write classification data frame to VCF file
#'
#' Add probabilities from classification data frame to a VCF file
#'
#' @param file An VCF file
#' @param classification.df Classification data frame with mutations as a column and posteriors for signatures. Must have 96 rows.
#'
#' @examples
#'
#' library(MutationalPatterns)
#' library(tidyverse)
#' library(vcfR)
#' library(R.utils)
#' library(Biostrings)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' vcf.file <- system.file("extdata", "SIMULATED_SAMPLE_SBS4_1.vcf", package = "excerno")
#'
#' # Load in signatures
#' cosmic.sigs <- get_known_signatures()
#'
#' # Get signatures
#' signatures <- matrix(nrow = 96, ncol = 2)
#' signatures[,1] <- cosmic.sigs[,4]
#' signatures[,2] <- get_ffpe_signature()
#' rownames(signatures) <- get_mutation_types()
#' colnames(signatures) <- c("SBS4", "FFPE")
#'
#' # Get contributions
#' contribution <- matrix(nrow = 2, ncol = 1)
#' contribution[,1] <- c(0.5, 0.5)
#' rownames(contribution) <- c("SBS4", "FFPE")
#'
#' classification.df <- get_classification(signatures, contribution)
#' write_classification_to_vcf(vcf.file, classification.df)
#'
#' @export
write_classification_to_vcf <- function(file, classifications.df) {

  # Argument validations
  if (!is.character(file)) { stop("argument file must be type character") }
  if (is.null(classifications.df$mutations)) { stop("argument classifications.df must have column mutations") }

  # Read and create VCF file
  suppressWarnings(invisible(capture.output(vcf.data <- read.vcfR(file))))
  sample <- get_mutational_vector(file)

  sig.names <- colnames(classifications.df)[2:length(colnames(classifications.df))]

  # Adding info meta data
  prob.info <- paste("##INFO=<ID=PROB,Number=", toString(length(sig.names)), ",Type=Integer,Description=\"Posteriors for each signature (", sep = "")

  for (sign in sig.names) {
    prob.info <- paste(prob.info, sign, ",", sep = "")
  }

  prob.info <- str_remove(prob.info, ",$")
  prob.info <- paste(prob.info, ")\">", sep = "")
  vcf.data@meta <- c(vcf.data@meta, prob.info)

  for (mut in 1:length(sample)) {

    # Get probabilities for a mutation type
    mut.row <- classifications.df %>%
      filter(mutations == sample[mut])

    # Generate a string with probabilites
    prob.str <- ""

    for (sign in sig.names) {
      prob.str <- paste(prob.str, formatC(mut.row[[sign]], digits = 2, format = "f"), ",", sep = "")
    }

    prob.str <- str_remove(prob.str, ",$")

    # Insert prob.str to data frame
    vcf.data@fix[,"INFO"][mut] <- paste(vcf.data@fix[,"INFO"][mut], ";PROB=", prob.str, sep = "")

    # Insert pass filter to variants with higher than 0.5 probability for "FFPE"
    if (!is.null("FFPE")) {
      if (mut.row[["FFPE"]] < 0.5) {
        vcf.data@fix[,"FILTER"][mut] <- "PASS"
      }
    }
  }

  # Writing vcf file with new info
  vcf.file <- paste(str_remove(tail(str_split(file, "/")[[1]], n = 1), ".vcf$"), "_classified.vcf.gz", sep = "")
  write.vcf(vcf.data, vcf.file)
  gunzip(vcf.file, overwrite = TRUE)

  # Create list for outputting results
  output <- list()
  output$vcf.data <- vcf.data
  output$class.df <- classifications.df

  return (output)
}
