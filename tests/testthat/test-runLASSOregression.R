load("expected/runLASSOregression_expected.rda")
data("targetInput")
data("responseInput")

set.seed(333)
result_test <- runLASSOregression(TargetMatrix = targetInput, ResponseMatrix  = responseInput, repeats = 3)

test_that("Process targets", {
  expect_equal(result_test, result)
})

