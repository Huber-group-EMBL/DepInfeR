load("expected/runLASSORegression_expected.rda")
data("targetInput")
data("responseInput")

set.seed(333)
# test serial version
result_test <- runLASSORegression(TargetMatrix = targetInput, 
                                  ResponseMatrix  = responseInput, repeats = 3)
test_that("Calculate target importance values", {
  expect_equal(result_test, result)
})

set.seed(123)
#test multi-core version
result_test_para <- runLASSORegression(TargetMatrix = targetInput, 
                                  ResponseMatrix  = responseInput, 
                                  repeats = 3, cores = 4)
test_that("Calculate target importance values with parallel version", {
    expect_equal(result_test_para, result_para)
})
