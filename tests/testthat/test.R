# One Series

test_that("testfile lists", {
  
  d = readRDS("inst/extdata/no_ga.rds")
  
  x = fhg_estimate(d$Q, Y = d$Ymean)
  
  expect_true(length(x) == 9)
  expect_true(nrow(x) == 2)
  expect_true(all(x$type == "Y"))
  expect_true(which.min(x$nrmse) == 1)
})

# Two Series

test_that("testfile lists", {
  
  d = readRDS("inst/extdata/no_ga.rds")
  
  x = fhg_estimate(Q = d$Q, Y = d$Ymean, V = d$V)
  
  expect_true(length(x) == 9)
  expect_true(nrow(x) == 2)
  expect_true(all(x$type == "Y"))
  expect_true(which.min(x$nrmse) == 1)
})

# Three Series no EA

test_that("testfile lists", {
  
  d = readRDS("inst/extdata/with_ga.rds")
  
  x = fhg_estimate(Q = d$Q, Y = d$Ymean, V = d$V, TW= d$TW)
  
  expect_true(length(x) == 9)
  expect_true(nrow(x) == 2)
  expect_true(all(x$type == "Y"))
  expect_true(which.min(x$nrmse) == 1)
})

# Three series w/EA