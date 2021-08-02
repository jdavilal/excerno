library(excerno)

cosmic.sigs <- get_known_signatures()
cosmic.sig4 <- as.matrix(cosmic.sigs[,4])

test_that("Arguments are valid", {
  expect_error(find_signature_name("hello"), "argument signature must be class matrix")
  expect_error(find_signature_name(123), "argument signature must be class matrix")


  test.matirx <- matrix(data = c(1, 2, 3))
  expect_error(find_signature_name(test.matrix))
})

test_that("Output is correct", {
  expect_equal(find_signature_name(cosmic.sig4), "SBS4")
  expect_equal(find_signature_names(cosmic.sigs[,1:4]), c("SBS1", "SBS2", "SBS3", "SBS4"))
})
