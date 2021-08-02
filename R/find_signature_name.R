#' Find signature name
#'
#' Find the cosmic name of the signature that closest resembles a signature
#' @param signature A signature matrix
#' @return A string of the most similar signature
#' @examples
#'
#' # Load in a signature
#' cosmic.sigs <- get_known_signature()
#' signature.4 <- as.matrix(cosmic.sigs[,4])
#' rownames(signature.4) <- get_mutation_types()
#'
#' # Usage
#' find_signature_name(signature.4)
#'
#' # Create matrix with multiple signatures
#' signatures <- cosmic.sigs[,1:3]
#' find_signature_names(cosmic.sigs[,1:3])
#' @export
find_signature_name <- function(signature) {

  sig.matrix <- signature

  # Argument validation
  if (!is.matrix(sig.matrix)) { stop("argument signature must be class matrix") }
  if (length(sig.matrix) != 96) { stop("argument signature must be length 96") }
  if (!all(rownames(sig.matrix) == get_mutation_types())) { rownames(sig.matrix) <- get_mutation_types() }

  # Load in cosmic signatures
  cosmic.sigs <- get_known_signatures()
  rownames(cosmic.sigs) <- get_mutation_types()
  cosmic.sigs <- cbind(cosmic.sigs, get_ffpe_signature())
  colnames(cosmic.sigs)[length(colnames(cosmic.sigs))] <- "FFPE"

  most.similar <- 0
  best.sig <- ""

  for (i in 1:length(colnames(cosmic.sigs))) {

    # Compare previos biggest cosine similarity value with next signature cosine similarity value
    if (most.similar < cos_sim(as.vector(sig.matrix), cosmic.sigs[,i])) {
      most.similar <- cos_sim(as.vector(sig.matrix), cosmic.sigs[,i])
      best.sig <- colnames(cosmic.sigs)[i]
    }
  }

  return (best.sig)
}

#' @rdname find_signature_name
#' @export
find_signature_names <- function(signatures) {

  sig.names <- vector(mode = "character")

  for (i in 1:dim(signatures)[2]) {
    sig.names <- c(sig.names, find_signature_name(as.matrix(signatures[,i])))
  }

  return (sig.names)
}
