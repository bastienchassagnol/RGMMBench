test_that("entropy function", {
  expect_equal(compute_shannon_entropy(c(0.5, 0.5)), 1)
})


test_that("positive definitess", {
  cov_matrix <- matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  expect_true(is_positive_definite(cov_matrix))
})

test_that("MLE estimation of the covariance matrix", {
  n <- nrow(datasets::longley)
  expect_equal(cov_MLE(datasets::longley), cov(datasets::longley)*(n-1)/n)
})
