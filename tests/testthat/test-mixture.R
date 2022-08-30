test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("entropy function", {
  expect_equal(compute_shannon_entropy(c(0.5, 0.5)), 1)
})
