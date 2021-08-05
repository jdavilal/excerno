library(excerno)

vcf.file <- system.file("extdata", "SIMULATED_SAMPLE_SBS4_1.vcf", package = "excerno")

# Load in signatures
cosmic.sigs <- get_known_signatures()

# Get signatures
signatures <- matrix(nrow = 96, ncol = 2)
signatures[,1] <- cosmic.sigs[,4]
signatures[,2] <- get_ffpe_signature()
rownames(signatures) <- get_mutation_types()
colnames(signatures) <- c("SBS4", "FFPE")

# Get contributions
contribution <- matrix(nrow = 2, ncol = 1)
contribution[,1] <- c(0.5, 0.5)
rownames(contribution) <- c("SBS4", "FFPE")

classification.df <- get_classification(signatures, contribution)

test_that("Arugments are valid", {
  expect_error(write_classification_to_vcf(123, classification.df), "argument file must be type character")
  expect_error(write_classification_to_vcf(vcf.file, classification.df[2:3]), "argument classifications.df must have column mutations")
})

test_that("Output is right", {
  write_classification_to_vcf(vcf.file, classification.df)

  # Load in vcf file
  vcf.data <- read.vcfR("SIMULATED_SAMPLE_SBS4_1_classified.vcf")

  # Convert to data frames
  vcf.fix <- data.frame(vcf.data@fix)
  vcf.fix$CHROM <- as.numeric(vcf.fix$CHROM)
  vcf.fix$POS <- as.numeric(vcf.fix$POS)
  vcf.gt <- data.frame(vcf.data@gt)

  # Check values in vcf file
  expect_true(vcf.fix$CHROM[1] == 1)
  expect_true(vcf.fix$POS[1] == 10522)
  expect_true(vcf.fix$REF[1] == "C")
  expect_true(vcf.fix$ALT[1] == "T")
  expect_true(vcf.fix$INFO[1] == "SOMATIC;TRUTH=FFPE;PROB=0.07,0.93")
  expect_true(vcf.fix$QUAL[1] == 53)
  expect_true(vcf.fix$FILTER[1] == "PASS")
  expect_true(vcf.gt$FORMAT[1] == "GT:GQ")
  expect_true(vcf.gt$SAMPLE1[1] == "0/0:8")
  expect_true(vcf.gt$SAMPLE2[1] == "0/0:57")

  # Delete recently created file
  file.remove("SIMULATED_SAMPLE_SBS4_1_classified.vcf")
})
