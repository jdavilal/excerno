#' Get classification of each mutation
#'
#' Get classification of each mutation using Bayes' Theorem.
#'
#' @param samples List of mutation samples, each coming from a different mutational signature
#' @param signatures List of the mutational signature that each sample came from. Each element should correspond an element in argument samples at the same index.
#' @param signature.names Vector of the signature names. Each element should correspond an element in argument signatures at the same index.
#' @return A data frame with the classification of each mutation in the samples
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
#' samples <- list(sample.sig4, sample.ffpe)
#' signatures <- list(cosmic.sig4, ffpe.sig)
#'
#' classify_simulated_samples(samples, signatures)
#' classify_simulated_samples(samples, signatures, c("SBS4", "FFPE"))
#' @export
classify_simulated_samples <- function(samples, signatures, signature.names = c()) {

  # Argument validation
  if (length(samples) < 2) { stop("argument samples must contain more than 1 sample") }
  if (length(signatures) < 2) { stop("argument signatures must contain more than 1 signature") }
  if (length(signatures) != length(samples)) { stop(paste("argument signatures must be of length", length(samples))) }
  if (length(signature.names) != length(samples) && !is.null(signature.names)) { stop(paste("argument signature.name must be of length", length(samples))) }

  sig.names <- vector(mode = "character")

  # Finding signature names if not given
  if (is.null(signature.names)) {
    for (i in 1:length(signatures)) {
      sig.names[i] <- find_signature_name(signatures[[i]])
    }
  } else {
    sig.names <- signature.names
  }

  # Empty vectors for names of mutation sources (signatures) and randomly generated samples
  names <- vector(mode = "character")
  signature.sample <- vector(mode = "character")

  for (i in 1:length(samples)) {
    # Adds each signature name for every element/profile in the sample to empty vector
    names <- append(names, rep(sig.names[i], length(samples[[i]])))
    # Adds each element/profile in the sample to empty vector
    signature.sample <- append(signature.sample, samples[[i]])
  }

  # Puts the long generated mutations and their source into a data frame and names columns
  df.final <- data.frame(signature.sample, names)
  colnames(df.final) <- c("mutations", "truth")


  # Converts sample vector into a table of probabilities
  tab <- table(signature.sample)/sum(table(signature.sample))

  # Converts table of probabilities to a data frame and defines column names
  sample.probs.df <- data.frame(tab)
  colnames(sample.probs.df) <- c("mutations", "frequencies")

  # Defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
    left_join(x = x, y = y, by = by) %>%
      mutate_each(funs(replace(., which(is.na(.)), 0)))
  }

  # Creates a vector of mutation types
  mutations <- get_mutation_types()

  # Creates a data.frame of mutation types to join with sample probabilities
  df <- data.frame(mutations)

  # Joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  sample.df <- left_join_NA(df, sample.probs.df, by = "mutations")

  # Add pseudocount to column of probabilities
  sample.matrix<- as.matrix(sample.df$frequencies) + 0.0001
  rownames(sample.matrix)= mutations

  # Computes exposure matrix using sample matrix and true signature matrix
  linear.res <- fit_to_signatures(sample.matrix, do.call(cbind, signatures))
  # Add element with original signatures so that Bayes function accepts this object
  linear.res$signatures <-do.call(cbind, signatures)

  # Name the rows and columns that represent mutational signatures for compatibility with Bayes function
  colnames(linear.res$signatures) <- sig.names
  rownames(linear.res$contribution) <- sig.names


  # Empty list to hold probability of each signature given randomly generated mutation
  prob.list <- list()

  for (i in length(samples)){
    # Each element of empty list is vector of probabilities for every element/profile from mutational sample
    prob.list[[i]]<-vector()
  }

  for (m in mutations) {
    # Each mutation type (string) is sent to Bayes function along with "nmf" results
    probs <- extract_all_prob(m, linear.res$signatures, linear.res$contribution)
    for (n in 1:length(samples)){
      # Probability for given mutation is sent to all corresponding signatures then proceeds to next mutation
      prob.list[[n]][m] <- probs[[n]]
    }
  }

  # Convert posterior probability list into a dataframe
  prob.list.df = data.frame(prob.list)
  rownames(prob.list.df) <- NULL
  # Each column represents a signature
  colnames(prob.list.df) <- sig.names
  # Add column with mutation names
  prob.list.df$mutations <- mutations

  # Join posterior probabilities with randomized mutations+source (keep only probabilities that
  # Correspond to an existing randomized mutation)
  df.final <- left_join(df.final, prob.list.df, by = "mutations")

  # Empty vector to hold Bayesian classifier
  classify <- vector(mode = "character")

  for (row in 1:nrow(df.final)){
    # Store posterior probabilities in each row for all signatures
    probabilities <- df.final[row,3:(2+length(samples))]
    tryCatch(
      expr = {
        # Finds the index for the maximum probability and retrieves name of signature that is the max
        classifier <- colnames(probabilities)[which.max(probabilities[1,])]
        # Each element is signature name (string) that has the maximum posterior probability
        classify[row]<- classifier
      },
      error = function(e){
        # When mutation that should not have been generated appears, probabilities and classifier are NA
        classifier <- NA
      }
    )
  }

  # Add column to the data frame that represents Bayesian classifier
  df.final$classify <- classify

  # Add column to the data frame that determines if classifier is different from source
  df.final <- df.final %>%
    # Misclassification when the source string does not match the Bayesian classifier
    mutate(misclassification = ifelse(classify != truth, 1, 0))

  # Remove mutations that should not have been randomly generated
  df.final <- df.final %>%
    # Mutations that don't appear in a signature are evaluated as NA in posteriors & classification
    drop_na()

  df.final <- df.final[sample(nrow(df.final)),]

  return (df.final)

}
