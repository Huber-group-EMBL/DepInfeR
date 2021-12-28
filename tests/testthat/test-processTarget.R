load("expected/ProcessTargetResults_expected.rda")
data("targetMatrix")

process_test <- processTarget(targetMatrix)

test_that("Process targets", {
  expect_equal(process_test, ProcessTargetResults)
})

