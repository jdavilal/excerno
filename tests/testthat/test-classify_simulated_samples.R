library(excerno)

set.seed(10)

# Load in signatures
cosmic.sigs <- get_known_signatures()
cosmic.sig4 <- as.matrix(cosmic.sigs[,4])
cosmic.sig6 <- as.matrix(cosmic.sigs[,6])
ffpe.sig <- get_ffpe_signature()

# Create samples
sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 100)
sample.sig6 <- create_signature_sample_vector(cosmic.sig6, 100)
sample.ffpe <- create_signature_sample_vector(ffpe.sig, 100)

samples <- list(sample.sig4, sample.ffpe)
signatures <- list(cosmic.sig4, ffpe.sig)

test_that("Arguments are valid", {
  expect_error(classify_simulated_samples(list(), signatures), "argument samples must contain more than 1 sample")
  expect_error(classify_simulated_samples(samples, list()), "argument signatures must contain more than 1 signature")
  expect_error(classify_simulated_samples(list(sample.sig4, sample.sig6, ffpe.sig), signatures), "argument signatures must be of length")
  expect_error(classify_simulated_samples(samples, signatures, c("Hello")), "argument signature.name must be of length 2")
})

test_that("Output is correct", {
  classify.df <- classify_simulated_samples(samples, signatures)
  mut.row <- classify.df %>%
    filter(mutations == "A[C>T]A") %>%
    head(1)

  expect_equal(mut.row$truth, "FFPE")
  expect_equal(round(mut.row$SBS4, 2), 0.11)
  expect_equal(round(mut.row$FFPE, 2), 0.89)
  expect_equal(mut.row$classify, "FFPE")
  expect_equal(mut.row$misclassification, 0)
})
