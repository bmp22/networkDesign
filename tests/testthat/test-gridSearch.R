test_that("Simple Design Found Example 1", {
  A<-JRSSExamples$ex1
  expect_equal(gridSearch(A,p=2)$AOptVal, 0.4186047, tolerance=1e-7)
  expect_equal(gridSearch(A,p=2, isoSearch = TRUE)$AOptVal, 0.4186047, tolerance=1e-7)
})
test_that("Network Design Found Example 1", {
  expect_equal(gridSearch(JRSSExamples$ex1,2, networkEffects = TRUE)$AOptVal, 0.2369146, tolerance=1e-7)

})

test_that("MRSAR Design Found Example 1", {
expect_equal(gridSearch(JRSSExamples$ex1,2,weightPrior = c(20,10,0,0,0,0.5))$AOptVal,1.10239e-05,tolerance=1e-10)
  expect_equal(gridSearch(JRSSExamples$ex1,2,weightPrior = c(20,10,0,0,0,0.5), viralOpt=FALSE)$AOptVal,0.4000108,tolerance=1e-7)
})






