library(excerno)

# Load in signatures
cosmic.sigs <- get_known_signatures()

# Get signatures
signatures <- matrix(nrow = 96, ncol = 2)
signatures[,1] <- cosmic.sigs[,4]
signatures[,2] <- get_ffpe_signature()

# Get contributions
contribution <- matrix(nrow = 2, ncol = 1)
contribution[,1] <- c(0.5, 0.5)

# Argument validations for more than one sample
contributions <- matrix(nrow = 2, ncol = 4)
contributions[,] <- sample(0.5, 4, replace = TRUE)

test_that("Arguments are valid", {

  # Argument validations for signatures
  expect_error(get_classification(c(1, 2, 3, 4), contribution), "argument signatures must be class matrix")
  expect_error(get_classification(matrix(c(1,2,3,4)), contribution), "argument signatures must have a row length of 96")
  expect_error(get_classification(as.matrix(cosmic.sigs[,4]), contribution), "argument signatures must have more than one signature")
  expect_error(get_classification(signatures, contribution), "argument signatures has no row names")
  rownames(signatures) <- sample("test", 96, replace = TRUE)
  expect_error(get_classification(signatures, contribution), "row names for argument signatures does not equal get_mutation_types()")
  rownames(signatures) <- get_mutation_types()
  expect_error(get_classification(signatures, contribution), "argument signature does not have have column names defined")
  colnames(signatures) <- c("SBS4", "FFPE")

  # Argument validations for contributions
  expect_error(get_classification(signatures, c(1, 2, 3, 4)), "argument contribution must be class matrix")
  expect_error(get_classification(signatures, contribution), "argument contribution must have row names defined")
  rownames(contribution) <- c("TEST", "FFPE")
  expect_error(get_classification(signatures, contribution), "column names of argument signature must equal row names of argument contribution")
  rownames(contribution) <- c("SBS4", "FFPE")

  expect_error(get_classifications(signatures, contributions), "argument contribution does not have column names defined")
})

test_that("Output is correct", {

  colnames(signatures) <- c("SBS4", "FFPE")
  rownames(signatures) <- get_mutation_types()
  colnames(contribution) <- "SAMPLE1"
  rownames(contribution) <- c("SBS4", "FFPE")
  rownames(contributions) <- c("SBS4", "FFPE")
  colnames(contributions) <- c("SAMPLE1", "SAMPLE2", "SAMPLE3", "SAMPLE4")

  # Test for one sample
  classification.df <- get_classification(signatures, contribution)

  row.mut <- classification.df %>%
    filter(mutations == "A[C>T]A") %>%
    head(1)

  expect_true(round(row.mut$SBS4, 2) == 0.09)
  expect_true(round(row.mut$FFPE, 2) == 0.91)

  # Test for multiple samples
  classification.df <- get_classifications(signatures, contributions)

  row.mut <- classification.df[[1]] %>%
    filter(mutations == "A[C>T]A") %>%
    head(1)
  expect_true(round(row.mut$SBS4, 2) == 0.09)
  expect_true(round(row.mut$FFPE, 2) == 0.91)

  row.mut <- classification.df[[2]] %>%
    filter(mutations == "A[C>T]A") %>%
    head(1)
  expect_true(round(row.mut$SBS4, 2) == 0.09)
  expect_true(round(row.mut$FFPE, 2) == 0.91)

  row.mut <- classification.df[[3]] %>%
    filter(mutations == "A[C>T]A") %>%
    head(1)
  expect_true(round(row.mut$SBS4, 2) == 0.09)
  expect_true(round(row.mut$FFPE, 2) == 0.91)

  row.mut <- classification.df[[4]] %>%
    filter(mutations == "A[C>T]A") %>%
    head(1)
  expect_true(round(row.mut$SBS4, 2) == 0.09)
  expect_true(round(row.mut$FFPE, 2) == 0.91)
})
