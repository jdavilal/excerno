#' Create GRange object from sample data frame
#'
#' Generate a GRange object from sample data frame with additional columns for a vcf file
#'
#' @param sample.df A data frame with columns 'mutations' and 'truth'. See output of \code{\link{classify_simulated_samples()}} for example
#' @param seq A DNAString. The sequence of a chromosome.
#' @param chromosome A string of the chromosome
#' @param info A vector for the info column
#' @param quality A vector for the quality column
#' @param filter A vector for the filter column
#' @param format A vector for the format column
#' @param samples A list of vectors for the sample columns
#' @param samples.names A vector of strings for the sample column headers
#' @return A GRange object with the information inserted
#' @examples
#'
#' # Load in signatures
#' cosmic.sigs <- get_known_signatures()
#' cosmic.sig4 <- as.matrix(cosmic.sigs[,4])
#' ffpe.sig <- get_ffpe_signature()
#'
#' # Create samples
#' sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 100)
#' sample.ffpe <- create_signature_sample_vector(ffpe.sig, 100)
#'
#' # Create classification data frame
#' samples <- list(sample.sig4, sample.ffpe)
#' signatures <- list(cosmic.sig4, ffpe.sig)
#' classify.df <- classify_simulated_samples(samples, signatures)
#'
#' seq <- getSeq(Hsapiens, "chr1")
#'
#' create_gr_from_sample(classify.df, seq, "chr1")
#' create_gr_from_sample(classify.df, seq, "chr1", info, quality, filter, format, samples, sample.names)
#' @export
create_gr_from_sample <- function(sample.df, seq, chromosome, info = c(), quality = c(), filter = c(), format = c(), samples = list(), sample.names = c()) {

  # Argument sample validations
  if (!is.data.frame(sample.df)) { stop("argument sample.df is not a dataframe") }
  if (!is.element("mutations", colnames(sample.df))) { stop("argument sample.df must have column mutations") }
  if (!is.element("truth", colnames(sample.df))) { stop("argument sample.df must have column truth") }

  # Argument seq validations
  if (class(seq) != "DNAString") { stop("argument seq must be class DNAString") }

  # Argument chromosome validations
  if (!is.character(chromosome)) { stop("argument chromosome is not type character") }
  if (!str_detect(chromosome, "chr[\\d]+")) { stop("argument chromosome is not of format: chr# e.g chr17") }

  # Argument extra columns validations
  num.rows <- dim(sample.df)[1]
  if (length(info) != 0 && !is.character(info)) { stop("argument info must be vector of class character") }
  if (length(info) != 0 && length(info) != num.rows) { stop(paste("argument info must be vector of length", toString(num.rows))) }

  if (length(quality) != 0 && !is.numeric(quality)) { stop("argument quality must be vector of class numeric") }
  if (length(quality) != 0 && length(quality) != num.rows) { stop(paste("argument quality must be vector of length", toString(num.rows))) }

  if (length(filter) != 0 && !is.character(filter)) { stop("argument filter must be vector of class character") }
  if (length(filter) != 0 && length(filter) != num.rows) { stop(paste("argument filter must be vector of length", toString(num.rows))) }

  if (length(format) != 0 && !is.character(format)) { stop("argument format must be vector of class character") }
  if (length(format) != 0 && length(format) != num.rows) { stop(paste("argument format must be vector of length", toString(num.rows))) }

  if (length(samples) != 0 && !is.list(samples)) { stop("argument samples must be class list") }
  if (length(samples) != 0 && !is.character(samples[[1]])) { stop("argument samples must be list of class character") }
  if (length(samples) != 0 && length(samples[[1]]) != num.rows) { stop(paste("argument samples must be vector of length", toString(num.rows))) }

  position <- vector(mode = "double")
  info <- vector(mode = "character")
  # Original nucleotide
  ref <- vector(mode = "character")
  # New nucleotide
  alt <- vector(mode = "character")

  print("Determining position of mutations...")

  # Iteration begins at 10000 to avoid telomeres
  i <- 10000

  # Loops through every mutation in sample
  for (m in 1:length(sample.df$mutations)) {

    # Reformat mutation string and stores substitutions into ref and alt
    mut.reformatted <- str_replace_all(sample.df$mutations[m], "[[:punct:]>]", "")
    ref <- c(ref, str_sub(mut.reformatted, 2, 2))
    alt <- c(alt, str_sub(mut.reformatted, 3, 3))
    str_sub(mut.reformatted, 3, 3) <- ""

    if (length(info) == num.rows) {
      info[m] <- c(info[m], paste(";TRUTH=", sample.df$truth[m], sep = ""))
    } else {
      info <- c(info, paste("SOMATIC", ";TRUTH=", sample.df$truth[m], sep = ""))
    }

    # Searches chromosome sequence for mutation match
    while (mut.reformatted != toString(seq[i:(i+2)])) {
      i = i + 1
    }

    # Once matched, store position
    position <- c(position, i+1)

  }

  print("...Done!")

  # Creates GRange object with position, ref, and alt
  gr <- GRanges(Rle(c(chromosome), c(1)), IRanges(start = position, end = position))

  # Additional metacolumns added
  gr$ID <- "."
  gr$REF <- ref
  gr$ALT <- alt
  gr$INFO <- info

  # Adding "." to columns with no input
  if (length(quality) == 0) {
    gr$QUAL <- "."
  } else {
    gr$QUAL <- quality
  }
  if (length(filter) == 0) {
    gr$FILTER <- "."
  } else {
    gr$FILTER <- filter
  }
  if (length(format) == 0) {
    gr$FORMAT <- "."
  } else {
    gr$FORMAT <- format
  }

  if (length(samples) != 0) {
    for (i in 1:length(samples)) {
      gr@elementMetadata[sample.names[i]] <- samples[[i]]
    }
  }

  return (gr)
}
