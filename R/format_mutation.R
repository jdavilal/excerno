
format_mutation <- function(position, sequence, alternate = "") {

  # Argument validation
  if (!is.double(position)) { stop("argument position is not of type double") }
  if (class(sequence) != "DNAString") { stop("argument sequence is not of type DNAString") }
  if (!is.character(alternate)) { stop("argument alternate is not of type character") }
  if (str_length(alternate) > 1) { stop("argument alternate must have only one character or none at all") }

  alternate <- toupper(alternate)
  if (!(alternate == "A" || alternate == "T" || alternate == "C" || alternate == "G")) { stop("argument alternate must be A, T, C, or G only") }

  # Returns value base on argument alternate
  if (str_length(alternate) == 1) {
    mutation <- paste(toString(sequence[position - 1]), "[", toString(sequence[position]), ">", alternate, "]", toString(sequence[position + 1]), sep = "")
  } else {
    mutation <- paste(toString(sequence[position - 1]), toString(sequence[position]), toString(sequence[position + 1]), sep = "")
  }

  return (mutation)
}
