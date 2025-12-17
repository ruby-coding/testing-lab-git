#install.packages(c("devtools", "testthat"))
library(testthat)
source("Final_Project.R")

test_that("Final_Project.R runs without error", {
  expect_silent(source("Final_Project.R"))
})

source("Final_Project.R")

test_that("wjs_matrix columns sum to 1 and values are valid probabilities", {
  W <- wjs_matrix(qh = 0.7, qc = 0.8, max_s = 5)
  expect_equal(ncol(W), 5)
  expect_true(all(W >= 0))
  expect_true(all(W <= 1))
  expect_equal(colSums(W), rep(1, 5), tolerance = 1e-12)
})

test_that("simulate_counts preserves household totals per size", {
  set.seed(1)
  D <- simulate_counts(qh = 0.7, qc = 0.8, obs_counts = tec_77)
  expect_equal(dim(D), dim(tec_77))
  expect_equal(colSums(D), colSums(tec_77))
})

test_that("distance_two_outbreaks returns a finite nonnegative number", {
  set.seed(2)
  theta <- c(0.7, 0.8, 0.6, 0.85)
  d <- distance_two_outbreaks(theta, tec_77, tec_80)
  expect_true(is.numeric(d))
  expect_length(d, 1)
  expect_true(is.finite(d))
  expect_true(d >= 0)
})

test_that("abc_rejection_4param returns expected structure and parameter bounds", {
  set.seed(3)
  res <- abc_rejection_4param(tec_77, tec_80, N = 2000, accept_frac = 0.05)
  expect_true(is.list(res))
  expect_true(all(c("theta", "dist", "eps") %in% names(res)))
  expect_true(is.matrix(res$theta))
  expect_equal(ncol(res$theta), 4)
  expect_true(all(res$theta >= 0))
  expect_true(all(res$theta <= 1))
  expect_true(all(res$dist <= res$eps))
})
