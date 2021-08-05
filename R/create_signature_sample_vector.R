#'Create signature sample vector
#'
#' Simulate sample to match signature distribution (vector of all mutations)
#' @param signature A signature matrix
#' @param num.mut Number of mutations generated
#'
#' @return A vector of type character
#' @examples
#'
#' library(MutationalPatterns)
#' library(tidyverse)
#'
#' # Load signature from COSMIC
#' cosmic.sigs <- get_known_signatures()
#' cosmic.sig4 <- as.matrix(cosmic.sigs[,4])
#'
#' # Create vector
#' create_signature_sample_vector(cosmic.sig4, 100)
#' @export
create_signature_sample_vector <- function(signature, num.mut = 100){

  # Argument validation
  if (!is.matrix(signature)) { stop("argument signature must be type matrix or numeric vector") }
  if (num.mut < 50) { warning("Low amount of mutations will not reflect original signature as well")}

  #creates a vector of mutation types
  mutations <- get_mutation_types()
  #creates a vector of probabilities of each mutation type
  probabilities = as.vector(signature[,1])

  #creates sample of desired number of total mutations to match the distribution of mutations in the signature
  signature.sample <- sample(mutations, size = num.mut, replace = TRUE, prob = probabilities)

  #returns a vector with length = number.of.mutations, each element is a mutation type (string)
  return(signature.sample)
}
