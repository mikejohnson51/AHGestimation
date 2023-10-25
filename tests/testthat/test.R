data = AHGestimation::nwis

library(dplyr)

test_that("One Series", {
  
  x =  select(data, Q, Y) %>% 
     ahg_estimate()
  
  expect_true(length(x) == 6)
  expect_true(nrow(x) == 2)
  expect_true(all(x$type == "Y"))
  expect_true(which.min(x$nrmse) == 1)
})

test_that("Two Series", {
  
  x = select(data, Q, Y, V = V) %>% 
    ahg_estimate()
  
  expect_true(length(x) == 9)
  expect_true(nrow(x) == 1)
  
})


test_that("Three Series", {
  
  x =  df = select(data, Q, Y, V, TW) %>% 
    ahg_estimate()
  
  expect_true(length(x) == 15)
  expect_true(nrow(x) == 4)
  expect_true(which.min(filter(x, viable)$tot_error) == 1)
  
  y = select(data, Q, Y, V, TW) %>% 
    ahg_estimate(verbose = TRUE)
  
  expect_true(length(y) == 15)
  expect_true(nrow(y) == 4)
})

test_that("Date Filter", {
  
  expect_warning( date_filter(select(data, Q), year = 10) )
  
  expect_true(nrow(date_filter(select(data, date, Q), year = 10)) < nrow(data))
  
  expect_true(nrow(date_filter(select(data, date, Q), year = 2)) < nrow(date_filter(select(data, date, Q), year = 10)))
  expect_true(nrow(date_filter(select(data, date, Q), year = 1)) < nrow(date_filter(select(data, date, Q), year = 1, keep_max = TRUE)))
  
  expect_equal(ncol(date_filter(select(data, date, Q), year = 10)), 2)
  
  expect_true(nrow(date_filter(data, year = 10)) < nrow(data))
  
  expect_warning(date_filter(data, year = 0))
  
})

test_that("MAD", {
  expect_true(nrow(mad_filter(data)) < nrow(data))
  expect_true(nrow(mad_filter(data,3)) < nrow(mad_filter(data,4)))
  expect_warning(mad_filter(data,1))
  
  expect_no_error( mad_filter(df = select(data, Q, Y)))
})

test_that("NLS", {
  expect_true(nrow(nls_filter(data)) < nrow(data))
  expect_true(nrow(nls_filter(data,.1)) < nrow(nls_filter(data,.7)))
  expect_warning(nls_filter(data,.01))
  
   expect_no_error(nls_filter(df = select(data, Q, Y)))
})

test_that("qva", {
  expect_true(nrow(qva_filter(data)) < nrow(data))
  expect_true(nrow(qva_filter(data,.01)) < nrow(qva_filter(data,.05)))
  expect_warning(qva_filter(data,.00001))
  
  expect_warning(qva_filter(df = select(data, Q, Y)))
})

test_that("sig", {
  expect_equal(significance_check(data), data)
  expect_error(significance_check(data, pvalue = 1e-1000))
  
  # One relation
  expect_equal(significance_check(df = select(data, Q,Y)), select(data, Q,Y))
})


test_that("hydraulics", {
 
  x =  select(data, Q, Y, V, TW) %>% 
    ahg_estimate()
  
  h = compute_hydraulic_params(x)
  
  expect_equal(nrow(h), 4)
  expect_equal(ncol(h), 7)
  expect_equal(round(h$r[1], 2), 5.01)
  
  expect_equal(round(compute_n(data),2), .14)
  
  cs = cross_section(h$r[1])
  
  expect_equal(nrow(cs), 30)
  expect_equal(ncol(cs), 4)

})

