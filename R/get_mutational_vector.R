
get_mutational_vector <- function(vcf.data) {

  # Argument validation
  if (class(vcf.data) != "vcfR") { stop("argument vcf.data is not type vcfR") }

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

get_mutational_vectors <- function(vcf.data) {

  sample.vectors <- list()

  for (i in 1:length(vcf.data)) {
    sample.vector <- get_mutational_vector(vcf.data[[i]])
    sample.vectors[[i]] <- sample.vector
  }

  return (sample.vectors)
}
