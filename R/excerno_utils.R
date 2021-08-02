#' Get mutation types
#'
#' @return a vector of the 96 mutation types
#' @export
get_mutation_types <- function() {
  file <- system.file("utils", "mut_sig.order.txt", package = "excerno")
  data <- read.table(file, sep = "\t", header = T)

  return (data$ext.context)
}

#' Get our FFPE signature
#'
#' @return a matrix of our ffpe signature
#' @export
get_ffpe_signature <- function() {

  # Read data
  file <- system.file("utils", "ffpe.signature.txt", package = "excerno")
  data <- read.table(file, sep = "\t", header = T)

  # Format data into matrix
  ffpe.matrix <- matrix(nrow = 96)
  rownames(ffpe.matrix) <- get_mutation_types()
  ffpe.matrix[,1] <- data$prob

  return (ffpe.matrix)
}

# Print log message
print_info <- function(message) {
  output <- paste("INFO [", str_split(date(), " ")[[1]][4], "] " , message, sep = "")
  print(output)
}
