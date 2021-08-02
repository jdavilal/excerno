library(excerno)

test.file <- system.file("extdata", "simulated_sample_sig4.vcf", package = "excerno")

test_that("Arguments are valid", {
  expect_error(get_mutational_vector(c(1, 2, 3, 4, 5)), "argument vcf.data is not type vcfR")
})

test_that("Output is correct", {

  # Load file for testing
  file <- system.file("extdata", "SIMULATED_SAMPLE_SBS4_1.vcf", package = "excerno")
  vcf.data <- read.vcfR(file)

  # Load in correct output
  cosmic.sigs <- get_known_signatures()
  cosmic.sig4 <- as.matrix(cosmic.sigs[, 4])
  ffpe.sig <- get_ffpe_signature()
  set.seed(10)
  sample.ffpe <- create_signature_sample_vector(ffpe.sig, 500)
  sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 500)

  test.vector <- c(sample.ffpe, sample.sig4)
  vcf.vector <- get_mutational_vector(vcf.data)

  expect_true(all(table(test.vector) == table(vcf.vector)))
})
