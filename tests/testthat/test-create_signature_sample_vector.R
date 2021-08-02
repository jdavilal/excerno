library(excerno)

cosmic.sigs <- get_known_signatures()
cosmic.sig4 <- as.matrix(cosmic.sigs[,4])

test_that("Arguments are valid", {
  expect_error(create_signature_sample_vector("cosmic.sig4", 100), "argument signature must be type matrix or numeric vector")
  expect_warning(create_signature_sample_vector(cosmic.sig4, 10), "Low amount of mutations will not reflect original signature as well")
})

test_that("Output value is correct", {
  set.seed(10)
  mut.vector <- create_signature_sample_vector(cosmic.sig4, 200)
  expect_equal(round(signature_cosine_similarity(mut.vector, cosmic.sig4), 2), 0.95)
})
