#' Extract probabilities using bayes' thereom
#'
#' extract_all_prob() implements bayes' thereom to get the probabilites of a mutation coming from a specific signature
#'
#' @param mutation Input string of the mutation type to look at. Formatted examples: C[C>A]T, G[T>C]A
#' @param signatures Matrix of the signatures relevent in sample
#' @param contributions Matrix of how much each signature contributed to the samples
#' @return A vector of the posteriors
#' @examples
#'
#'cosmic.sigs <- get_known_signatures()
#'
#' # Get signatures
#' signatures <- matrix(nrow = 96, ncol = 2)
#' signatures[,1] <- cosmic.sigs[,4]
#' signatures[,2] <- get_ffpe_signature()# Naming columns and rows
#' colnames(signatures) <- c("SBS4", "FFPE")
#' rownames(signatures) <- get_mutation_types()
#'
#' # Get contributions
#' contribution <- matrix(nrow = 2, ncol = 1)
#' contribution[,1] <- c(0.5, 0.5)
#' rownames(contribution) <- c("SBS4", "FFPE")
#'
#' mutation <- "A[C>T]A"
#' extract_all_prob(mutation, signatures, contribution)
#' @export
extract_all_prob <- function(mutation, signatures, contribution) {

  # Argument validation for mutation
  if (!is.character(mutation)) { stop("argument mutation must be type character")}
  if (!str_detect(mutation, "[ACTG]\\[[CT]>[ACTG]\\][ACTG]")) { stop("argument mutation must be formatted: *[*>*]*")}

  # Argument validation for signatures
  if (dim(signatures)[1] != 96) { stop("argument signatures must have a row length of 96") }
  if (is.null(rownames(signatures))) ( stop("argument signatures has no row names") )
  if (!all(rownames(signatures) == get_mutation_types())) { stop("Row names for argument signatures must equal get_mutation_types()") }

  if (is.null(rownames(contribution))) ( stop("argument contribution has no row names") )
  if (!all(colnames(signatures) == rownames(contribution))) { stop("argument contribution row names and argument signatures column names does not match") }

  # Creates empty list of length equal to amount of signatures
  prob = vector(mode = "list", length=length(colnames(signatures)))

  # initailize accumulator variable and denominator of Bayes Theorem
  i <- 1
  denominator = 0

  # Fills prob with prior and likelihood for each signature
  for (signature in colnames(signatures)) {
    # Naming each element of prob
    names(prob)[i] <- signature

    # Pr(Mutation|Signature) - likelihood
    a <- signatures[mutation,signature]/sum(signatures[,signature])

    # Pr(Signature) - prior
    b <- contribution[signature,]/sum(contribution)

    # Appending likelihood and prior into prob
    prob[[i]] <- c(a, b)

    # Adding product of likelihood and prior to denominator
    denominator = denominator + a * b

    # Accumulating
    i <- i + 1
  }

  # Intialize acummulator and vector to store posteriors
  i <- 1
  posteriors <- vector(mode = "list", length = length(colnames(signatures)))

  # Fills posteriors
  for (x in colnames(signatures)) {
    names(posteriors)[i] <- x
    # Calculates numerator of Bayes formula for each signature
    numerator <- prob[[x]][1] * prob[[x]][2]
    # Finishes Bayes calculation and inserts value into the vector posteriors
    posteriors[[i]] <- numerator / denominator

    # Accumulating
    i <- i + 1
  }

  return (posteriors)
}
