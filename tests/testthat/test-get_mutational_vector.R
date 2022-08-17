library(excerno)

test_that("Arguments are valid", {
  expect_error(get_mutational_vector(c(1, 2, 3, 4, 5)), "argument vcf.file is not type character")
})

test_that("Output is correct", {

  # Load file for testing
  file <- system.file("extdata", "SIMULATED_SAMPLE_SBS4_1.vcf", package = "excerno")

  # Load in correct output
  cosmic.sigs <- get_known_signatures()
  cosmic.sig4 <- as.matrix(cosmic.sigs[, 4])
  ffpe.sig <- get_ffpe_signature()
  set.seed(10)
  sample.ffpe <- create_signature_sample_vector(ffpe.sig, 500)
  sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 500)

  test.vector <- c(sample.ffpe, sample.sig4)
  vcf.vector <- get_mutational_vector(file)

  # The following tests are approximate since the simulated vcf
  # contains repeated positions which are collapsed when using
  # read_vcfs_as_granges
  expect_true(sum(table(test.vector) - table(vcf.vector))<100)

  # Load files for inputing multiple vcfR objects
  vcf.files <- list.files(
    system.file("extdata", package = "excerno"),
    pattern = "SIMULATED_SAMPLE_SBS4_\\d.vcf", full.names = TRUE)

  sample.vectors <- get_mutational_vectors(vcf.files)

  set.seed(10)
  sample.ffpe <- create_signature_sample_vector(ffpe.sig, 500)
  sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 500)
  test.vector <- c(sample.ffpe, sample.sig4)
  expect_true(sum(table(test.vector) - table(sample.vectors[[1]]))<100)

  set.seed(20)
  sample.ffpe <- create_signature_sample_vector(ffpe.sig, 800)
  sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 200)
  test.vector <- c(sample.ffpe, sample.sig4)
  expect_true(sum(table(test.vector) - table(sample.vectors[[2]]))<100)

  set.seed(30)
  sample.ffpe <- create_signature_sample_vector(ffpe.sig, 400)
  sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 600)
  test.vector <- c(sample.ffpe, sample.sig4)
  expect_true(sum(table(test.vector) - table(sample.vectors[[3]]))<100)

  set.seed(40)
  sample.ffpe <- create_signature_sample_vector(ffpe.sig, 100)
  sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 900)
  test.vector <- c(sample.ffpe, sample.sig4)
  expect_true(sum(table(test.vector) - table(sample.vectors[[4]]))<100)
})
