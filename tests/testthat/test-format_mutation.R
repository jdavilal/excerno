library(excerno)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

seq <- getSeq(Hsapiens, "chr1")

test_that("argument validation works", {
  expect_error(format_mutation("string", seq, "A"))
  expect_error(format_mutation(1000, "sequence", "A"))
  expect_error(format_mutation("string", seq, "AD"))
})

test_that("format_mutation returns the correct string", {
  expect_equal(format_mutation(10560, seq, "T"), "A[C>T]A")
  expect_equal(format_mutation(19213, seq, "C"), "C[A>C]A")
  expect_equal(format_mutation(20124, seq, "A"), "A[G>A]G")
})
