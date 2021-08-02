library(excerno)

set.seed(10)

# Load in signatures
cosmic.sigs <- get_known_signatures()
cosmic.sig4 <- as.matrix(cosmic.sigs[,4])
ffpe.sig <- get_ffpe_signature()

# Create samples
sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 100)
sample.ffpe <- create_signature_sample_vector(ffpe.sig, 100)

# Create classification data frame
samples <- list(sample.sig4, sample.ffpe)
signatures <- list(cosmic.sig4, ffpe.sig)
classify.df <- classify_simulated_samples(samples, signatures)

seq <- getSeq(Hsapiens, "chr1")

test_that("Arguments are valid", {

  # Validation for argument sample
  expect_error(create_gr_from_sample(list(), seq, "chr1"), "argument sample.df is not a dataframe")
  expect_error(create_gr_from_sample(classify.df[2:4], seq, "chr1"), "argument sample.df must have column mutations")

  # Validation for argument seq
  expect_error(create_gr_from_sample(classify.df, sample(1:10)), "argument seq must be class DNAString")

  # Validation for argument chromosome
  expect_error(create_gr_from_sample(classify.df, seq, 123), "argument chromosome is not type character")
  expect_error(create_gr_from_sample(classify.df, seq, "chrone"), "argument chromosome is not of format: chr# e.g chr17")
  expect_error(create_gr_from_sample(classify.df, seq, "crh1"), "argument chromosome is not of format: chr# e.g chr17")

  # Validation for extra columns
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", info = c(1)), "argument info must be vector of class character")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", info = c("SOMATIC")), "argument info must be vector of length 200")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", quality = c("Hello")), "argument quality must be vector of class numeric")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", quality = c(1, 2)), "argument quality must be vector of length 200")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", filter = c(1)), "argument filter must be vector of class character")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", filter = c("PASS")), "argument filter must be vector of length 200")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", format = c(1)), "argument format must be vector of class character")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", format = c("GT:GQ")), "argument format must be vector of length 200")

  samples.list <- list(c("10:30", "10:20"), c("20:40", "20:30"))
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", samples = c(1)), "argument samples must be class list")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", samples = list(c(12, 12), c(12, 12))), "argument samples must be list of class character")
  expect_error(create_gr_from_sample(classify.df, seq, "chr1", samples = list(c("10:30", "10:20"), c("20:40", "20:30"))), "argument samples must be vector of length 200")
})

test_that("Output is correct", {
  set.seed(10)
  classify.gr <- data.frame(create_gr_from_sample(classify.df, seq, "chr1"))

  expect_true(classify.gr$seqnames[1] == "chr1")
  expect_true(classify.gr$start[1] == 10553)
  expect_true(classify.gr$REF[1] == "C")
  expect_true(classify.gr$ALT[1] == "A")
  expect_true(classify.gr$INFO[1] == "SOMATIC;TRUTH=SBS4")
  expect_true(classify.gr$QUAL[1] == ".")
  expect_true(classify.gr$FILTER[1] == ".")
  expect_true(classify.gr$FORMAT[1] == ".")

  # Adding values to other columns
  info <- sample("SOMATIC", 200, replace = TRUE)
  quality <- sample(50:100, 200, replace = TRUE)
  filter <- sample("PASS", 200, replace = TRUE)
  format <- sample("GT:GQ", 200, replace = TRUE)
  samples <- list(sample(paste("0/0:", 1:100, sep = ""), 200, replace = TRUE), sample(paste("0/0:", 1:100, sep = ""), 200, replace = TRUE))
  sample.names <- c("SAMPLE1", "SAMPLE2")
  set.seed(10)
  classify.gr <- data.frame(create_gr_from_sample(classify.df, seq, "chr1", info, quality, filter, format, samples, sample.names))

  expect_true(classify.gr$seqnames[1] == "chr1")
  expect_true(classify.gr$start[1] == 10553)
  expect_true(classify.gr$REF[1] == "C")
  expect_true(classify.gr$ALT[1] == "A")
  expect_true(classify.gr$INFO[1] == "SOMATIC;TRUTH=SBS4")
  expect_true(classify.gr$QUAL[1] == 78)
  expect_true(classify.gr$FILTER[1] == "PASS")
  expect_true(classify.gr$FORMAT[1] == "GT:GQ")
  expect_true(classify.gr$SAMPLE1[1] == "0/0:90")
  expect_true(classify.gr$SAMPLE2[1] == "0/0:13")
})
