fpath <- system.file("testdata", "ProcessTargetResults_expected.rda", package="DepInfeR")
load(fpath)
data("targetMatrix")

#without force keeping targets
process_test <- processTarget(targetMatrix)

test_that("Process targets", {
  expect_equal(process_test, ProcessTargetResults)
})

#force keeping targets
process_test_keepTar <- processTarget(targetMatrix,
                                      keepTargets = c("BTK","CCNT1"))
test_that("Process targets with preserving specified targets", {
    expect_equal(process_test_keepTar, ProcessTargetResults_keepTar)
})

#without removing correlated features
process_test_noRemove <- processTarget(targetMatrix,
                                      removeCorrelated = FALSE)
test_that("Process targets with preserving specified targets", {
    expect_equal(process_test_noRemove, ProcessTargetResults_noRemove)
})
