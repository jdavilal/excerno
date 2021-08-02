library(excerno)

cosmic.sigs <- get_known_signatures()

# Get signatures
signatures <- matrix(nrow = 96, ncol = 2)
signatures[,1] <- cosmic.sigs[,4]
signatures[,2] <- get_ffpe_signature()

# Get contributions
contribution <- matrix(nrow = 2, ncol = 1)
contribution[,1] <- c(0.5, 0.5)

test_that("Arguments are valid", {

  # Argument validations for mutation
  expect_error(extract_all_prob(123, signatures, contribution), "argument mutation must be type character")
  expect_error(extract_all_prob("1234", signatures, contribution), "argument mutation must be formatted: *[*>*]*")
  expect_error(extract_all_prob("1[1>1]1", signatures, contribution), "argument mutation must be formatted: *[*>*]*")
  expect_error(extract_all_prob("A[B>C]D", signatures, contribution), "argument mutation must be formatted: *[*>*]*")

  # Argument validations for signatures
  expect_error(extract_all_prob("A[C>T]A", matrix(c(1,2,3,4), nrow = 4, ncol = 1), contribution), "argument signatures must have a row length of 96")
  expect_error(extract_all_prob("A[C>T]A", signatures, contribution), "argument signatures has no row names")
  rownames(signatures) <- sample("test", 96, replace = TRUE)
  expect_error(extract_all_prob("A[C>T]A", signatures, contribution), "Row names for argument signatures must equal get_mutation_types()")

  # Argument validations for contribution
  rownames(signatures) <- get_mutation_types()
  colnames(signatures) <- c("SBS4", "FFPE")
  expect_error(extract_all_prob("A[C>T]A", signatures, contribution), "argument contribution has no row names")
  rownames(contribution) <- sample("test", 2, replace = TRUE)
  expect_error(extract_all_prob("A[C>T]A", signatures, contribution), "argument contribution row names and argument signatures column names does not match")
})

test_that("Output is correct", {

  # Naming columns and rows
  colnames(signatures) <- c("SBS4", "FFPE")
  rownames(signatures) <- get_mutation_types()

  rownames(contribution) <- c("SBS4", "FFPE")

  mutation <- "A[C>T]A"
  probs <- extract_all_prob(mutation, signatures, contribution)

  expect_true(round(probs$SBS4, 2) == 0.09)
  expect_true(round(probs$FFPE, 2) == 0.91)
})
