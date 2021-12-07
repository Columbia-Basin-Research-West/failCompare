bad_time_dat=c(rep(5,3),rep(10,3))

test_that("'SEs' argument switches to FALSE when model cannot be fit",{
  my_e=quote(taglife.fn_weib3(tags.in = bad_time_dat,model.in = "weibull",tag.se = T))
  my_e2=quote(taglife.fn_weib3(tags.in = bad_time_dat,model.in = "weibull",tag.se = F))
  expect_identical(try_fit(my_call = my_e),try_fit(my_call = my_e2))})

