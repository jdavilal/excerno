
get_classification <- function(signatures, contribution) {

  output <- data.frame(mutations = get_mutation_types())
  prob.list <- list()
  sig.names <- colnames(signatures)

  # Each element of empty list is vector of probabilities for every element/profile from mutational sample
  for (name in sig.names) {
    prob.list[[name]] <- vector()
  }

  # Calculate bayes' theorem for all mutations
  for (m in get_mutation_types()) {
    probs <- extract_all_prob(m, 1, signatures, contribution)

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

get_classifcations <- function(signatures, contribution) {

  # Argument validations
  if (!is.matrix(contribution)) { stop("argument contribution must be class matrix") }
  if (is.null(colnames(signatures))) { stop("argument signature must have column names defined") }
  if (is.null(rownames(contribution))) { stop("argument contribution must have row names defined") }
  if (!all(colnames(signatures) == rownames(contribution))) { stop("column names of argument signature must equal row names of argument contribution") }
  if (length(colnames(contribution)) == 1) { warning("Suggest use of get_classification() instead for only one sample") }

  sample.size <- length(colnames(contribution))
  output.list <- list()

  for (i in 1:sample.size) {
    output.list[[colnames(contribution)[i]]] <- get_classification(signatures, as.matrix(contribution[,i]))
  }

  return (output.list)
}
