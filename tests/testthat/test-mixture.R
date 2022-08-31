test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("entropy function", {
  expect_equal(compute_shannon_entropy(c(0.5, 0.5)), 1)
})

# simulated_data <- rnmix_skewed_with_outliers(200, theta = list(
#   p = c(0.5, 0.5),
#   mu = c(0, 4), sigma = c(0.3, 0.3),
#   skew = c(0, 0)
# ))
