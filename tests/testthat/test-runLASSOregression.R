fpath <- system.file("testdata", "runLASSORegression_expected.rda", package="DepInfeR")
load(fpath)

data("targetInput")
data("responseInput")

# test serial version
param <- BiocParallel::MulticoreParam(workers = 1, RNGseed = 333)
result_test <- runLASSORegression(TargetMatrix = targetInput, 
                                  ResponseMatrix  = responseInput, 
                                  repeats = 3, BPPARAM = param)
test_that("Calculate target importance values", {
  expect_equal(result_test, result)
})
