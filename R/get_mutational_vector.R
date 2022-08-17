#' Get mutational vector
#'
#' Get the mutational vector from a vcf file
#'
#' @param vcf.data The name of the vcf file
#' @return A vector of mutation strings
#' @examples
#'
#' library(MutationalPatterns)
#' library(tidyverse)
#' library(vcfR)
#' library(Biostrings)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Load file for testing
#' file <- system.file("extdata", "SIMULATED_SAMPLE_SBS4_1.vcf", package = "excerno")
#'
#' # Load in correct output
#' cosmic.sigs <- get_known_signatures()
#' cosmic.sig4 <- as.matrix(cosmic.sigs[, 4])
#' ffpe.sig <- get_ffpe_signature()
#'
#' sample.ffpe <- create_signature_sample_vector(ffpe.sig, 500)
#' sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 500)
#'
#' test.vector <- c(sample.ffpe, sample.sig4)
#' vcf.vector <- get_mutational_vector(file)
#'
#' # Load files for inputing multiple vcfR objects
#' vcf.files <- list.files(
#'   system.file("extdata", package = "excerno"),
#'   pattern = "SIMULATED_SAMPLE_SBS4_\\d.vcf", full.names = TRUE)
#'
#' sample.vectors <- get_mutational_vectors(vcf.files)
#' @export
get_mutational_vector <- function(vcf.file, ref_genome="BSgenome.Hsapiens.UCSC.hg38") {

  # Argument validation
  if (!is.character(vcf.file)) { stop("argument vcf.file is not type character") }


  suppressWarnings(
    invisible(
      capture.output(
        grl <-read_vcfs_as_granges(c(vcf.file),c(vcf.file),ref_genome,
                                   type=c("snv"),
                                   predefined_dbs_mbs =TRUE))))

  # Get the context of the mutations
  type_context <- type_context(grl[[1]], ref_genome)
  ctxt_vec <- as.character(type_context$context)
  mut_vec <- type_context$types
  # Create the string with the right format e.g. G[C>T]T
  pre <- str_sub(ctxt_vec,1,1)
  post <- str_sub(ctxt_vec,3,3)
  mutations <- str_c(pre,"[",mut_vec,"]",post)
  return (mutations)
}


#' @rdname get_mutational_vector
#' @export
get_mutational_vectors <- function(vcf.files, ref_genome="BSgenome.Hsapiens.UCSC.hg38") {
    lapply(vcf.files, get_mutational_vector, ref_genome=ref_genome)
}
