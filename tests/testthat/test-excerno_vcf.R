library(excerno)

# Load in signatures
cosmic.sigs <- get_known_signatures()
cosmic.sig4 <- as.matrix(cosmic.sigs[,4])
ffpe.sig <- get_ffpe_signature()

# Load in vcf files
vcf.files <- list.files(system.file("extdata", package = "excerno"), pattern = "SIMULATED_SAMPLE_SBS4_\\d.vcf", full.names = TRUE)
vcf.file <- "SIMULATED_SAMPLE_SBS4_1_classified.vcf"

artifact <- "FFPE"
method <- "nmf"
num.signatures <- 2
target.sigs <- matrix(nrow = 96, ncol = 2)
target.sigs[,1] <- cosmic.sig4
target.sigs[,2] <- ffpe.sig
rownames(target.sigs) <- get_mutation_types()
colnames(target.sigs) <- c("SBS4", "FFPE")

test_that("Arguments are valid", {
  expect_error(excerno_vcf(c(1, 2, 3, 4), artifact, method), "argument files must be type character")
  expect_error(excerno_vcf(vcf.files, c(1, 2), method), "argument artifact must be type character")
  expect_error(excerno_vcf(vcf.files, artifact, "hello"), "argument method must be \"nmf\" or \"linear\"")
  expect_error(excerno_vcf(vcf.files, artifact, "linear"), "argument target.sigs must be non-empty if method equals linear")
  expect_error(excerno_vcf(vcf.file, artifact, "nmf"), "argument files must be length greater than 2 if method equals nm")

})


test_that("VCF files are correct", {

  # Testing nmf method
  excerno_vcf(vcf.files, artifact, "nmf")

  vcf.data <- read.vcfR(vcf.file)

  expect_equal(vcf.data@fix[, "INFO"][[1]], "SOMATIC;TRUTH=FFPE;PROB=0.10,0.90")
  expect_equal(vcf.data@fix[, "POS"][[1]], "10522")
  expect_equal(vcf.data@fix[, "REF"][[1]], "C")
  expect_equal(vcf.data@fix[, "QUAL"][[1]], "53")
  expect_equal(vcf.data@fix[, "FILTER"][[1]], "PASS")
  expect_equal(vcf.data@fix[, "ALT"][[1]], "T")
  expect_equal(vcf.data@fix[, "CHROM"][[1]], "1")
  expect_equal(vcf.data@gt[, "FORMAT"][[1]], "GT:GQ")
  expect_equal(vcf.data@gt[, "SAMPLE1"][[1]], "0/0:8")
  expect_equal(vcf.data@gt[, "SAMPLE2"][[1]], "0/0:57")

  # Testing nmf method
  excerno_vcf(vcf.files, artifact, "linear", target.sigs = target.sigs)

  vcf.data <- read.vcfR(vcf.file)

  expect_equal(vcf.data@fix[, "INFO"][[1]], "SOMATIC;TRUTH=FFPE;PROB=0.07,0.93")
  expect_equal(vcf.data@fix[, "POS"][[1]], "10522")
  expect_equal(vcf.data@fix[, "REF"][[1]], "C")
  expect_equal(vcf.data@fix[, "QUAL"][[1]], "53")
  expect_equal(vcf.data@fix[, "FILTER"][[1]], "PASS")
  expect_equal(vcf.data@fix[, "ALT"][[1]], "T")
  expect_equal(vcf.data@fix[, "CHROM"][[1]], "1")
  expect_equal(vcf.data@gt[, "FORMAT"][[1]], "GT:GQ")
  expect_equal(vcf.data@gt[, "SAMPLE1"][[1]], "0/0:8")
  expect_equal(vcf.data@gt[, "SAMPLE2"][[1]], "0/0:57")

  file.remove("SIMULATED_SAMPLE_SBS4_1_classified.vcf")
  file.remove("SIMULATED_SAMPLE_SBS4_2_classified.vcf")
  file.remove("SIMULATED_SAMPLE_SBS4_3_classified.vcf")
  file.remove("SIMULATED_SAMPLE_SBS4_4_classified.vcf")
})
