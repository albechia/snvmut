library(testthat)

context("Tests for snv_count")


test_that("The input format is correct", {
  wrong_vector <- c("ACT]T", "A[C>T]C", "T[T>C]C", "T[C>TG", "C[>T]T")
  expect_error(snv_count(wrong_vector),
               "The vector must have the same structure of the vectors returned by snv_extraction")
})

test_that("The output is correct when providing an imput that is correct", {
  good_vector <- c("A[C>T]A", "G[C>T]T", "A[C>T]A")
  expected_output <- data.frame(SNV = c("C>T"), Count = c(3), Percentage = c(100))
  expect_equal(snv_count(good_vector), expected_output)
})
