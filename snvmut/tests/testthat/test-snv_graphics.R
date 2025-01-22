library(testthat)

library(ggplot2)
library(patchwork)

context("Tests for snv_graphics")

test_that("snv_graphics stops with incorrect inputs", {
  incorrect_df <- data.frame(ABC = c("A>T", "G>C"), X = c(10, 20), "7" = c(33.3, 66.7))
  expect_error(snv_graphics(incorrect_df))
})

test_that("snv_graphics returns a plot object for correct input", {
  correct_df <- data.frame(SNV = c("A>T", "G>C"), Count = c(10, 20),Percentage = c(33.3, 66.7))
  plot <- snv_graphics(correct_df)
  expect_true(is.ggplot(plot))
})
