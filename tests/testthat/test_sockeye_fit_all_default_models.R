# library(testthat)
data(sockeye)
taglife=sockeye[,1]
f_all=fc_fit(time=taglife,model="all")

modnms=c('weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma','vitality.ku','vitality.4p')
f_indiv=sapply(modnms,fc_fit,time=taglife,simplify = F)

test_that("default models fit together and separately for sockeye",{
  expect_s3_class(f_all,class="fc_list")})


test_that("fc_combine converts full set of models to a fc_list",{
  tmp=fc_combine(f_indiv)
  expect_s3_class(tmp,class="fc_list")
  # expect_identical(object = tmp,expected = f_all)
  })

test_that("get_param_nm returns same results for single model as for vector",{
  expect_identical(
    object = rbind(get_param_nm(mod_nm=c("weibull3","gompertz"))),
    expected = get_param_nm(mod_nm=c("weibull3","gompertz")))
})

test_that("get_param_nm errors when duplicate models or improper names are entered",{
  expect_error(get_param_nm(mod_nm=c("weibull3","weibull3","gompertz")))
  expect_error(get_param_nm(mod_nm=letters[1:5]))
})
