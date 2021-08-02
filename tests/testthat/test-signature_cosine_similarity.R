library(excerno)

set.seed(10)
cosmic.sig4 <- as.matrix(get_known_signatures()[,4])
sample.sig4 <- create_signature_sample_vector(cosmic.sig4)

test_that("Arguments are valid", {
  expect_error(signature_cosine_similarity(list(c(1, 2, 3)), cosmic.sig4), "argument mutations.vector must be type character")
  expect_error(signature_cosine_similarity(c("2", "2", "3"), c("2", "2", "3")), "argument signature.matrix must be class matrix")
  expect_error(signature_cosine_similarity(c("2", "2", "3"), matrix(ncol = 97)), "argument signature.matrix must be length 96")
})

test_that("Output is correct", {
  expect_equal(round(signature_cosine_similarity(sample.sig4, cosmic.sig4), 2), 0.91)
})
