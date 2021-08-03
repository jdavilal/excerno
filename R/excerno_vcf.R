#' Exercerno VCF
#'
#' excerno_vcf() produces filtered vcf files. It uses NMF or nonnegative linear combination of mutation signatures to determine contribution of signatures in samples. Then Bayes' Theroem is used to classify each variant.
#'
#' @param files VCF files
#' @param artifact The signature to consider as an artifact
#' @param method A string. The method used to determine the signatures (if not given) and calculate the contributions of each signature
#' @param num.signatures Number of signatures. Necessary arugment for when method "linear" is choosen.
#' @param target.sig Matrix of the target signatures.Necessary arugment for when method "linear" is choosen.
#'
#' @return Object containing the vcf objects and classification data frame
#' @examples
#'
#' # Load in files
#' files <- list.files( system.file("extdata", package = "excerno"), pattern = "simulated_sample_sig4..", full.names = TRUE)
#'
#' exerno_vcf(files)
#' @export
excerno_vcf <- function(files, artifact = c(), method = "nmf", num.signatures = 2, target.sigs = c()) {

  # Arguments validation
  if (!is.character(files)) { stop("argument files must be type character") }
  if (!is.character(artifact)) { stop("argument artifact must be type character") }
  if (method != "nmf" && method != "linear") { stop("argument method must be \"nmf\" or \"linear\"") }
  if (method == "linear" && is.null(target.sigs)) { stop("argument target.sigs must be non-empty if method equals linear") }
  if (method == "nmf" && length(files) < 2) { stop("argument files must be length greater than 2 if method equals nmf") }

  # Read vcf files
  vcf.data <- list()
  num.samples <- length(files)

  for (i in 1:length(files)) {
    # Hide print
    suppressWarnings(invisible(capture.output(vcf.data[i] <- read.vcfR(files[i]))))
  }

  print_info("Creating mutational vectors")
  samples <- get_mutational_vectors(vcf.files)

  # Perform method to get present signatures and contributions of each sample
  if (method == "nmf") {
    print_info("Performing NMF")

    # Argument validation for method = nmf

    samples.df <- data.frame(mutations = get_mutation_types())

    # Add signatures to df
    for (i in 1:num.samples) {
      sample.tbl <- table(samples[[i]])/sum(table(samples[[i]]))
      sample.df <- data.frame(mutations = rownames(sample.tbl))
      sample.df[paste("SAMPLE", toString(i), sep = "")] <- as.vector(sample.tbl)
      samples.df <- left_join(samples.df, sample.df, by = "mutations")
    }

    samples.df[is.na(samples.df)] <- 0

    # Convert data frame to matrix for NMF
    sample.matrix <- as.matrix(samples.df[2:(num.samples + 1)])
    rownames(sample.matrix) <- get_mutation_types()

    # NMF and renaming columns and rows
    nmf.res <- extract_signatures(sample.matrix, rank = num.signatures, nrun = 10)
    sig.names <- find_signature_names(nmf.res$signatures)
    colnames(nmf.res$signatures) <- sig.names
    rownames(nmf.res$contribution) <- sig.names

    signatures <- nmf.res$signatures
    contribution <- nmf.res$contribution
  } else if (method == "linear") {
    print_info("Performing Linear combination method")

    # Defines a function that joins two data frames and replaces NA with 0
    left_join_NA <- function(x, y, by) {
      left_join(x = x, y = y, by = by) %>%
        mutate_each(funs(replace(., which(is.na(.)), 0)))
    }

    # Creates a vector of all 96 snv mutation types
    mutations <- get_mutation_types()
    num.signatures <- length(colnames(target.sigs))
    sig.names <- colnames(target.sigs)

    contribution <- matrix(nrow = num.signatures, ncol = num.samples)
    rownames(contribution) <- sig.names

    for (i in 1:num.samples) {

      # Converts mutations_vector into a table of probabilities
      tab <- table(samples[[i]])/sum(table(samples[[i]]))

      # Converts table of probabilities to a data frame and defines column names
      mutation.probs.df <- data.frame(tab)
      colnames(mutation.probs.df) <- c("mutations", "frequencies")

      # Creates a data.frame of mutation types to join with sample probabilities
      df <- data.frame(mutations)

      # Joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
      mut.profile.df <- left_join_NA(df, mutation.probs.df, by = "mutations")

      # Transform mut.profile.df to mutational profile matrix, add pseudocount to column of probabilities
      mut.profile<- as.matrix(mut.profile.df$frequencies) + 0.0001
      rownames(mut.profile)= mutations

      fits.res <- fit_to_signatures(mut.profile, target.sigs)

      contribution[,i] <- fits.res$contribution
    }
    signatures <- target.sigs
    colnames(contribution) <- paste("SAMPLE", seq(1:num.samples), sep = "")
  }

  print_info("Generating classification data frames")
  # Determine signature probabilities from each sample
  classifications.df <- get_classifications(signatures, contribution)

  print_info("Loading in values to vcf files")
  # Insert probabilities into original vcfR object
  # samples[[1]][1]formatC(x, digits = 8, format = "f")
  for (i in 1:num.samples) {

    # Adding info meta data
    prob.info <- paste("##INFO=<ID=PROB,Number=", toString(length(sig.names)), ",Type=Integer,Description=\"Posteriors for each signature (", sep = "")

    for (sign in sig.names) {
      prob.info <- paste(prob.info, sign, ",", sep = "")
    }

    prob.info <- str_remove(prob.info, ",$")
    prob.info <- paste(prob.info, ")\">", sep = "")
    vcf.data[[i]]@meta <- c(vcf.data[[i]]@meta, prob.info)

    for (mut in 1:length(samples[[i]])) {

      # Get probabilities for a mutation type
      mut.row <- classifications.df[[i]] %>%
        filter(mutations == samples[[i]][mut])

      # Generate a string with probabilites
      prob.str <- ""

      for (sign in sig.names) {
         prob.str <- paste(prob.str, formatC(mut.row[[sign]], digits = 2, format = "f"), ",", sep = "")
      }

      prob.str <- str_remove(prob.str, ",$")

      # Insert prob.str to data frame
      vcf.data[[i]]@fix[,"INFO"][mut] <- paste(vcf.data[[i]]@fix[,"INFO"][mut], ";PROB=", prob.str, sep = "")

      # Insert pass filter to variants with higher than 0.5 probability for artifacts
      if (!is.null(artifact)) {
        if (mut.row[[artifact]] < 0.5) {
          vcf.data[[i]]@fix[,"FILTER"][mut] <- "PASS"
        }
      }
    }

    # Writing vcf file with new info
    vcf.file <- paste(str_remove(tail(str_split(files[i], "/")[[1]], n = 1), ".vcf$"), "_classified.vcf.gz", sep = "")
    write.vcf(vcf.data[[i]], vcf.file)
    gunzip(vcf.file, overwrite = TRUE)
  }

  # Create list for outputting results
  output <- list()
  output$vcf.data <- vcf.data
  output$class.df <- classifications.df

  return (output)
}
