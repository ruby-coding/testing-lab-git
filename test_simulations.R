context("Check simulations")
source("simulate_reads.R")
test_that("true seq creation function works", {
  expect_true(all(colSums(make_true_seqs(10, 5)) > 0))
  expect_true(all(colSums(make_true_seqs(20, 3)) < 3))
  expect_true(all(make_true_seqs(15, 4) %in% c(0,1)))
})