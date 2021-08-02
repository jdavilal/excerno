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

set.seed(10)
info <- sample("SOMATIC", 200, replace = TRUE)
quality <- sample(50:100, 200, replace = TRUE)
filter <- sample("PASS", 200, replace = TRUE)
format <- sample("GT:GQ", 200, replace = TRUE)
samples <- list(sample(paste("0/0:", 1:100, sep = ""), 200, replace = TRUE), sample(paste("0/0:", 1:100, sep = ""), 200, replace = TRUE))
sample.names <- c("SAMPLE1", "SAMPLE2")
set.seed(10)
classify.gr <- create_gr_from_sample(classify.df, seq, "chr1", info, quality, filter, format, samples, sample.names)

test_that("Arguments are valid", {
  expect_error(write_grange_to_vcf("classify.gr", "new_file.vcf"), "argument grange is not class GRange")
  expect_error(write_grange_to_vcf(classify.gr, 123), "argument file.name is not type character")
})

test_that("Written files are correct", {
  file.name <- "new_file.vcf"
  write_grange_to_vcf(classify.gr, file.name)

  # Load in vcf file
  vcf.data <- read.vcfR(file.name)

  # Convert to data frames
  vcf.fix <- data.frame(vcf.data@fix)
  vcf.fix$CHROM <- as.numeric(vcf.fix$CHROM)
  vcf.fix$POS <- as.numeric(vcf.fix$POS)
  vcf.gt <- data.frame(vcf.data@gt)

  # Check values in vcf file
  expect_true(vcf.fix$CHROM[1] == 1)
  expect_true(vcf.fix$POS[1] == 10553)
  expect_true(vcf.fix$REF[1] == "C")
  expect_true(vcf.fix$ALT[1] == "A")
  expect_true(vcf.fix$INFO[1] == "SOMATIC;TRUTH=SBS4")
  expect_true(vcf.fix$QUAL[1] == 78)
  expect_true(vcf.fix$FILTER[1] == "PASS")
  expect_true(vcf.gt$FORMAT[1] == "GT:GQ")
  expect_true(vcf.gt$SAMPLE1[1] == "0/0:90")
  expect_true(vcf.gt$SAMPLE2[1] == "0/0:13")

  # Delete recently created file
  file.remove(file.name)
})
