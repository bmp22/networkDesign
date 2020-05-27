test_that("Simple Design Found Example 1", {
  expect_equal(gridSearch(JRSSExamples$ex1,2)$AOptVal, 0.4186047, tolerance=1e-7)
})
test_that("Network Design Found Example 1", {
  expect_equal(gridSearch(JRSSExamples$ex1,2, networkEffects = TRUE)$AOptVal, 0.2369146, tolerance=1e-7)

})


