fpath <- system.file("testdata", "runLASSORegression_expected.rda", package="DepInfeR")
load(fpath)

data("targetInput")
data("responseInput")

# test serial version
result_test <- runLASSORegression(TargetMatrix = targetInput, 
                                  ResponseMatrix  = responseInput, 
                                  repeats = 3, RNGseed = 333)
test_that("Calculate target importance values", {
  expect_equal(result_test, result)
})
