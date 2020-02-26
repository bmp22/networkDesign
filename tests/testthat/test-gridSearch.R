test_that("Simple Design Found", {
  expect_equal(gridSearch(JRSSExamples$ex1,2)$AOptVal, 0.4186047, tolerance=1e-7)
})
