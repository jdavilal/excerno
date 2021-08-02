#' Get mutational vector
#'
#' Get the mutational vector from a vcf file
#'
#' @param vcf.data The name of the vcf file
#' @return A vector of mutation strings
#' @examples
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
get_mutational_vector <- function(vcf.file) {

  # Argument validation
  if (!is.character(vcf.file)) { stop("argument vcf.file is not type character") }

  vcf.data <- read.vcfR(vcf.file)

  # Convert to data frames
  vcf.fix <- data.frame(vcf.data@fix)
  vcf.fix$CHROM <- as.numeric(vcf.fix$CHROM)
  vcf.fix$POS <- as.numeric(vcf.fix$POS)
  vcf.gt <- data.frame(vcf.data@gt)

  mutations <- vector(mode = "character")
  seq <- list()
  mut_num <- length(vcf.fix[[1]])

  # Iterate through the mutations
  for (i in 1:mut_num) {
    row <- vcf.fix[i, ]

    # Loads in chromosome sequence if needed
    if (row$CHROM > length(seq)) {
      seq[[paste("CHR", toString(row$CHROM), sep = "")]] <- getSeq(Hsapiens, paste("chr", "1", sep = ""))
    }

    # Formatting
    mutation <- format_mutation(row$POS, seq[[row$CHROM]], toString(row$ALT))
    mutations <- c(mutations, mutation)
  }

  return (mutations)
}

#' @rdname get_mutational_vector
#' @export
get_mutational_vectors <- function(vcf.files) {

  sample.vectors <- list()

  for (i in 1:length(vcf.files)) {
    sample.vector <- get_mutational_vector(vcf.files[[i]])
    sample.vectors[[i]] <- sample.vector
  }

  return (sample.vectors)
}
