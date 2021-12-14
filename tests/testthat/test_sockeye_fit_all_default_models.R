
library(failCompare)

data(sockeye)
taglife=sockeye[,1]

# model name vectors
flex_mods=names(failCompare:::fc_mod_ls)[names(failCompare:::fc_mod_ls) %in% names(flexsurv::flexsurv.dists)]
vit_mods=names(failCompare:::fc_mod_ls)[grep(names(failCompare:::fc_mod_ls),pattern="vitality")]
all_mods=c(flex_mods,vit_mods,"weibull3")

my_km=fc_fit(taglife,model = "kaplan-meier")
# my_km$KM_DF
# my_km$KM_mod


fmod=fc_fit(time = taglife,model=c("weibull","vitality.ku"),non_cen = rep(TRUE,length(taglife)))

nn=NULL
# fc_combine(list(aa,bb,nn))


############################################# #
### Testing fc_combine() warning
############################################# #

aa=fc_fit(time = taglife,
          model="weibull",
          SEs=T,
          non_cen = rep(TRUE,length(taglife)))


bb=fc_fit(time = taglife,
          model="weibull3",
          Hess=T,
          non_cen = rep(TRUE,length(taglife)))


cc=fc_combine(list(aa,bb))
cc$par_tab


fc_fit_single(y = taglife,
              y_sfrac=fc_surv(taglife),
              model="weibull",
              Hess=T,
              non_cen = rep(TRUE,length(taglife)),
              KM_DF=my_km$KM_DF,
              KM_mod=my_km$KM_mod)



test_that("<2 model obj in fc_combine cause an error",
          {expect_error(fc_combine(mod_ls = list(aa)))})

test_that("duplicate model names cause an error",
          {expect_error(fc_combine(mod_ls = list(aa,aa)))})

############################################# #
### Testing individual fc_fit() expresssions
############################################# #

# checking flexsurv mods
lapply(flex_mods,function(modnms){
  fc_fit(time = taglife,
                # y_sfrac=fc_surv(taglife),
                model=modnms,
                Hess=T,
                non_cen = rep(TRUE,length(taglife)))})

lapply(flex_mods,function(modnms){
  fc_fit(time = taglife,
         # y_sfrac=fc_surv(taglife),
         model=modnms,
         Hess=F,
         non_cen = rep(TRUE,length(taglife)))})

# checking weibull3 model
fc_fit(time = taglife,
              # y_sfrac=fc_surv(taglife),
              model="weibull3",
              Hess=T,
              non_cen = rep(TRUE,length(taglife)))

fc_fit(time = taglife,
              # y_sfrac=fc_surv(taglife),
              model="weibull3",
              Hess=F,
              non_cen = rep(TRUE,length(taglife)))

### vitality mods
lapply(vit_mods,function(modnms){
  fc_fit_single(y = taglife,
                y_sfrac=fc_surv(taglife),
                model=modnms,
                Hess=T,
                non_cen = rep(TRUE,length(taglife)),
                KM_DF=my_km$KM_DF,
                KM_mod=my_km$KM_mod)})


########################################## #
### Testing fc_fit_single() expresssions
########################################## #

# checking flexsurv mods
lapply(flex_mods,function(modnms){
  fc_fit_single(y = taglife,
                y_sfrac=fc_surv(taglife),
                model=modnms,
                Hess=T,
                non_cen = rep(TRUE,length(taglife)),
                KM_DF=my_km$KM_DF,
                KM_mod=my_km$KM_mod)})

lapply(flex_mods,function(modnms){
  fc_fit_single(y = taglife,
                y_sfrac=fc_surv(taglife),
                model=modnms,
                Hess=F,
                non_cen = rep(TRUE,length(taglife)),
                KM_DF=my_km$KM_DF,
                KM_mod=my_km$KM_mod)})


# checking weibull3 model
fc_fit_single(y = taglife,
              y_sfrac=fc_surv(taglife),
              model="weibull3",
              Hess=T,
              non_cen = rep(TRUE,length(taglife)),
              KM_DF=my_km$KM_DF,
              KM_mod=my_km$KM_mod)

fc_fit_single(y = taglife,
              y_sfrac=fc_surv(taglife),
              model="weibull3",
              Hess=F,
              non_cen = rep(TRUE,length(taglife)),
              KM_DF=my_km$KM_DF,
              KM_mod=my_km$KM_mod)


### vitality mods
lapply(vit_mods,function(modnms){
  fc_fit_single(y = taglife,
                y_sfrac=fc_surv(taglife),
                model=modnms,
                Hess=T,
                non_cen = rep(TRUE,length(taglife)),
                KM_DF=my_km$KM_DF,
                KM_mod=my_km$KM_mod)})


##################################### #
### Testing fc_tryfit() expresssions
##################################### #

## Directly use these

flx_mds=list()
for(i in 1:length(flex_mods)){
  flx_mds[[i]]=quote(flexsurv::flexsurvreg(survival::Surv(time=taglife,event=rep(TRUE,length(taglife))) ~ 1,
                                           dist = eval(flex_mods[i]),hessian = T))}

# funcion calls evaluated by fc_tryfit()
test_e=c(
  flx_mds,
  quote(vitality::vitality.ku(taglife,sdata = fc_surv(taglife),rc.data = F,pplot =F,silent=T,se=T)),
  quote(vitality::vitality.4p(taglife,sdata = fc_surv(taglife),rc.data = F,pplot =F,silent=T,se=T)),
  quote(taglife.fn_weib3(taglife,model.in = "weibull",tag.se=T))
)

# evaluating function calls
lapply(1:9,function(x){fc_tryfit(fit_call = eval(test_e[[x]]),model = all_mods[x],Hess=T)})



##################################### #
### Testing fc_fit() expresssions
##################################### #

# Testing that optimizer setting pass through
aa=fc_fit(time = taglife,
          model="weibull",
          non_cen = rep(TRUE,length(taglife)),control=list(maxit=1))

aa=fc_fit(time = taglife,
          model="weibull",
          non_cen = rep(TRUE,length(taglife)),control=list(trace=1))


bb=fc_fit(time = taglife,
          model="weibull3",
          Hess=T,
          non_cen = rep(TRUE,length(taglife)),control=list(trace=1))


# library(testthat)
data(sockeye)
taglife=sockeye[,1]
# 
# f_all=suppressWarnings(fc_fit(time=taglife,model="all"))
# modnms=c('weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma','vitality.ku','vitality.4p')
# f_indiv=suppressWarnings(sapply(modnms,fc_fit,time=taglife,simplify = F))
# tmp=fc_combine(f_indiv)
# 
# test_that("default models fit together and separately for sockeye",{
#   expect_s3_class(f_all,class="fc_list")
#   expect_s3_class(tmp,class="fc_list")
#   })

test_that("get_param_nm returns same results for single model as for vector",{
  expect_identical(
    object = rbind(get_param_nm(mod_nm=c("weibull3","gompertz"))),
    expected = get_param_nm(mod_nm=c("weibull3","gompertz")))
})

test_that("get_param_nm errors when duplicate models or improper names are entered",{
  expect_error(get_param_nm(mod_nm=c("weibull3","weibull3","gompertz")))
  expect_error(get_param_nm(mod_nm=letters[1:5]))
})

test_that("fc_boot will not provide predictions without 'times'",
          {expect_error(fc_boot(weib_mod,nrep = 50))})



bad_time_dat=c(rep(5,3),rep(10,3))
my_e=quote(taglife.fn_weib3(tags.in = bad_time_dat,model.in = "weibull",tag.se = T))
my_e2=quote(taglife.fn_weib3(tags.in = bad_time_dat,model.in = "weibull",tag.se = F))
