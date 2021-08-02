library(excerno)

test.vector <- c("T[C>A]C", "C[C>A]T", "G[C>A]C", "A[C>A]T", "C[T>A]C", "T[C>A]C", "G[C>T]C", "C[T>A]G", "T[C>A]C", "T[C>A]T", "C[T>A]A",
  "C[T>A]T", "C[C>A]G", "G[C>G]A", "G[C>A]C", "C[C>A]T", "C[T>A]G", "T[C>A]A", "T[C>A]C", "C[T>G]G", "A[C>A]C", "C[C>A]A",
  "A[T>A]C", "T[C>A]A", "C[C>A]A", "T[T>G]G", "T[C>A]A", "A[T>G]G", "C[C>A]C", "C[C>A]G", "A[C>T]A", "C[C>A]T", "C[T>A]C",
  "C[C>A]G", "C[C>A]A", "T[C>A]C", "C[T>A]G", "T[C>T]A", "T[T>C]T", "C[C>A]T", "A[C>A]T", "G[C>A]G", "T[T>A]C", "C[C>T]C",
  "C[C>T]T", "C[C>A]T", "G[C>T]A", "A[T>C]A", "A[C>A]A", "C[C>A]T", "C[T>A]A", "C[C>T]C", "A[C>A]C", "T[T>A]T", "G[T>A]G",
  "C[C>A]T", "C[C>A]A", "T[T>A]T", "C[C>A]C", "A[C>A]A", "A[C>A]C", "C[C>T]T", "T[C>A]T", "C[C>A]G", "A[C>G]A", "T[C>A]T",
  "T[C>A]A", "G[C>A]T", "C[T>A]G", "C[C>A]C", "C[C>G]T", "A[C>A]C", "T[C>A]G", "A[C>T]A", "T[C>A]A", "C[C>A]A", "T[C>A]G",
  "C[C>A]A", "A[C>A]T", "C[C>A]C", "G[T>A]G", "C[C>A]A", "A[C>T]A", "C[T>A]A", "T[C>A]C", "C[C>A]C", "C[C>A]A", "C[C>A]A",
  "A[C>G]A", "T[C>A]C", "C[C>A]A", "T[C>A]T", "C[C>A]G", "A[T>C]A", "T[C>A]T", "C[C>A]C", "G[T>A]A", "C[T>C]G", "A[C>G]G", "C[T>A]A",
  "A[C>T]A", "G[C>T]A", "C[C>T]T", "A[C>T]A", "T[C>T]C", "T[C>T]G", "A[C>T]A", "C[C>T]T", "C[C>T]A", "T[C>T]T", "T[C>T]T",
  "C[C>T]T", "T[C>T]T", "T[C>T]T", "T[C>T]A", "C[C>T]T", "T[C>T]C", "C[C>T]A", "A[C>T]T", "A[C>T]C", "T[C>T]A", "C[C>T]C",
  "C[C>T]A", "A[C>T]A", "C[C>T]A", "G[C>T]C", "A[C>T]A", "C[C>T]G", "T[C>T]A", "A[C>T]A", "G[C>T]A", "A[C>T]A", "A[C>T]C",
  "A[C>T]C", "A[C>T]A", "T[C>T]C", "C[C>T]G", "A[C>T]C", "C[C>T]C", "C[C>T]G", "C[C>T]A", "G[C>T]A", "T[C>T]A", "C[C>T]C",
  "G[C>T]A", "A[C>T]C", "C[C>T]A", "G[C>T]A", "A[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]T", "A[C>T]T", "T[C>T]A", "A[C>T]A",
  "T[C>T]A", "A[C>T]A", "A[C>T]C", "G[C>T]A", "C[C>T]T", "G[C>T]A", "A[C>T]T", "C[C>T]A", "G[C>T]C", "G[C>T]A", "G[C>T]C",
  "T[C>T]A", "C[C>T]T", "T[C>T]G", "C[C>T]T", "T[C>T]T", "A[C>T]T", "A[C>T]G", "C[C>T]C", "T[C>T]C", "G[C>T]A", "A[C>T]T",
  "C[C>T]G", "T[C>T]T", "A[C>T]T", "G[C>T]C", "G[C>T]T", "G[C>T]T", "C[C>T]C", "T[C>T]A", "A[C>T]C", "G[C>T]C", "G[C>T]A",
  "G[C>T]T", "C[C>T]T", "A[C>T]G", "A[C>T]T", "A[C>T]A", "T[C>T]C", "T[C>T]A", "A[C>T]T", "C[C>T]T", "C[C>T]C", "T[C>T]C", "A[C>T]T")

test.file <- system.file("extdata", "simulated_sample_sig4.vcf", package = "excerno")

# test_that("output vector is correct", {
#   expect_equal(min(get_mutational_vector(test.file) == test.vector), 1)
# })
