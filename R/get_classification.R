#' Get classification using Bayes' Thereom
#'
#' Determine probability of where each mutation came from using Bayes' Thereom.
#'
#' @param signatures A matrix of signatures
#' @param contribution A matrix of the contributions coming from each signature
#' @return A data frame with the probabilities for each mutation.
#'
#' @examples
#'
#' # Load in signatures
#' cosmic.sigs <- get_known_signatures()
#'
#' # Get signatures
#' signatures <- matrix(nrow = 96, ncol = 2)
#' signatures[,1] <- cosmic.sigs[,4]
#' signatures[,2] <- get_ffpe_signature()
#' colnames(signatures) <- c("SBS4", "FFPE")
#' rownames(signatures) <- get_mutation_types()
#'
#' # Get contributions
#' contribution <- matrix(nrow = 2, ncol = 1)
#' contribution[,1] <- c(0.5, 0.5)
#' rownames(contribution) <- c("SBS4", "FFPE")
#'
#' # Function usage
#' get_classification(signatures, contribution)
#'
#' # For more than one sample
#' contributions <- matrix(nrow = 2, ncol = 4)
#' contributions[,] <- sample(0.5, 4, replace = TRUE)
#' get_classifications(signatures, contributions)
#'

get_classification <- function(signatures, contribution) {

  # Argument validations for signatures
  if (!is.matrix(signatures)) { stop("argument signatures must be class matrix")}
  if (dim(signatures)[1] != 96) { stop("argument signatures must have a row length of 96") }
  if (dim(signatures)[2] < 2) { stop("argument signatures must have more than one signature") }
  if (is.null(rownames(signatures))) { stop("argument signatures has no row names") }
  if (all(rownames(signatures) != get_mutation_types())) { stop("row names for argument signatures does not equal get_mutation_types()") }
  if (is.null(colnames(signatures))) { stop("argument signature does not have have column names defined") }

  # Argument validations for contribution
  if (!is.matrix(contribution)) { stop("argument contribution must be class matrix") }
  if (is.null(rownames(contribution))) { stop("argument contribution must have row names defined") }

  # Argument validations for both
  if (!all(colnames(signatures) == rownames(contribution))) { stop("column names of argument signature must equal row names of argument contribution") }

  output <- data.frame(mutations = get_mutation_types())
  prob.list <- list()
  sig.names <- colnames(signatures)

  # Each element of empty list is vector of probabilities for every element/profile from mutational sample
  for (name in sig.names) {
    prob.list[[name]] <- vector()
  }

  # Calculate bayes' theorem for all mutations
  for (m in get_mutation_types()) {
    probs <- extract_all_prob(m, signatures, contribution)

    # Assign probabilities to correct dataframe
    for (s in 1:length(sig.names)) {
      prob.list[[sig.names[s]]][m] <- probs[[sig.names[s]]]
    }
  }

  # Assign probabilities to output object
  for (name in sig.names) {
    output[[name]] <- prob.list[[name]]
  }

  return (output)
}

#' @rdname get_classification
#' @export
get_classifications <- function(signatures, contribution) {

  # Argument validations
  if (is.null(colnames(contribution))) { stop("argument contribution does not have column names defined") }
  if (length(colnames(contribution)) == 1) { warning("suggest use of get_classification() instead for only one sample") }

  sample.size <- length(colnames(contribution))
  output.list <- list()

  for (i in 1:sample.size) {
    output.list[[colnames(contribution)[i]]] <- get_classification(signatures, as.matrix(contribution[,i]))
  }

  return (output.list)
}
