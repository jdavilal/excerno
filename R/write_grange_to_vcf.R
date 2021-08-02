#' Write a GRange to a vcf file
#'
#' @param grange A grange object with columns id, ref, alt, qual, filter, info
#' @param file.name A string of the new file name. Make sure to end string with ".vcf"
#' @examples
#'
#' set.seed(10)
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
#' set.seed(10)
#' info <- sample("SOMATIC", 200, replace = TRUE)
#' quality <- sample(50:100, 200, replace = TRUE)
#' filter <- sample("PASS", 200, replace = TRUE)
#' format <- sample("GT:GQ", 200, replace = TRUE)
#' samples <- list(sample(paste("0/0:", 1:100, sep = ""), 200, replace = TRUE), sample(paste("0/0:", 1:100, sep = ""), 200, replace = TRUE))
#' sample.names <- c("SAMPLE1", "SAMPLE2")
#' set.seed(10)
#' classify.gr <- create_gr_from_sample(classify.df, seq, "chr1", info, quality, filter, format, samples, sample.names)
#'
#' # Funtion usage
#' file.name <- "new_file.vcf"
#' write_grange_to_vcf(classify.gr, file.name)
#' @export
write_grange_to_vcf <- function(grange, file.name) {

  # Argument validation
  if (class(grange) != "GRanges") { stop("argument grange is not class GRange") }
  if (!is.character(file.name)) { stop("argument file.name is not type character") }

  # Create empty file
  file.data <- ""
  write(file.data, file.name)

  # Loads in template for a vcf file and writes to new file
  template.path <- system.file("utils", "vcf-template.vcf", package = "excerno")
  file.data <- read_file(template.path)
  file.data <- str_remove_all(file.data, "\r")
  file.data <- str_remove(file.data, "\n$")

  # Appends to file.data meta column names
  for (name in colnames(mcols(grange))) {
    file.data <- paste(file.data,"\t", name, sep="")
  }

  # Loops through rows in grange
  for (x in 1:length(grange)) {

    # Append CHROM number and position to file.data
    file.data <- paste(file.data, "\n","1\t", start(ranges(grange)[x]), sep="")

    # Loops through each metacolumn and appends meta data to file.data
    for (col in 1:length(colnames(mcols(grange)))) {
      meta.data <- mcols(grange)[colnames(mcols(grange))[col]][[1]][x]
      file.data <- paste(file.data,"\t", meta.data, sep="")
    }
  }

  # Loads file.data to new file
  write(file.data, file.name)
}
