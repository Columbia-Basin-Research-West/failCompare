#' @title Nonparametric bootstrap of failure time model object
#'
#' @details Random sampling of the failure time data with replacement as a means for propagating uncertainty in predictions of survival probability and estimates of parameter sampling distributions (Tibshirani and Efron 1993).
#'
#' @param mod_obj failure time model object of class "fc_obj"
#' @param nrep number of resampling replicates
#' @param type  number of resamples
#' @param times times at which survival fraction will be estimated
#' @param tol optional tolerance setting for the estimated proportion of bootstrap data sets that cannot be fit default = 0.9
#' @param ... aruments passed to the optimizer
#'
#' @return if \code{type="pred"} survival fraction or proportion of failed subjects (nrep x times) is returned, and
#' if \code{type="par"} a matrix of bootstrap parameter estimates  dimensions (nrep x number of parameters)
#'
#' @examples 
#' data(sockeye)
#' taglife=sockeye[,"days"]
#' weib_mod=fc_fit(taglife,model="weibull")
#'
#' @seealso
#' \code{\link{fc_fit}}
#'
#' @references 
#' 
#' Efron, B. and Tibshirani, R. 1993 An Introduction to the
#' Bootstrap. Chapman and Hall, New York.
#' 
#' Townsend R., L., J. R. Skalski, P. Dillingham, T. W. Steig. 2006 Correcting bias in survival estimation resulting from tag failure in acoustic and radiotelemetry studies.
#'  Journal of Agricultural Biological and Environmental Statistics.11:183–196. 
#'
#' @export
#' 
#' 
#'
fc_boot=function(mod_obj,nrep,type="pred",times=NULL,tol=0.9,...){
  if(!(type %in% c("pred","par")) & length(type)!=1){stop("expect either 'pred' or 'par' as type arguments")}
  stopifnot(class(mod_obj)=="fc_obj")
  if(type=="pred" & is.null(times)){stop("numeric vector of 'times'required for survival predictions")}
  if(all(is.numeric(c(nrep,times,tol))) & length(nrep)!=1 & length(tol)!=1){stop("")}

  obs_dat=mod_obj[["times"]]
  counter=0
  par_mat=matrix(ncol=2,nrow=0)
  pred_mat=NULL
  message("bootstrap initialized")
  
  # matching arguments from the original model fitting function
  call=match.call(expand.dots = F)
  args=formalArgs(eval(mod_obj$model_choice))
  str(call)
  str(args)
  
  print(match(names(call),args))
  
  
  if(mod_obj$mod_choice %in% c("vitality.ku","vitality.4p")){noHess="SE=F"}
  else{noHess="hessian=F"}
        # Initial evaluation
  bt=proc.time()
  while(max(1,nrow(par_mat)) <= 50){
    counter=counter+1
    ind=sample(1:nrow(obs_dat),replace = T)
    t_raw=obs_dat[ind,1]

    t_ord=order(obs_dat[ind,1]) # sorting times
    t_new=fc_surv(t_raw[t_ord])
    boot_dat=data.frame(ind=ind[t_ord],
                        time=t_raw[t_ord],
                        surv_frac=t_new,
                        non_cen=(obs_dat$non_cen[ind])[t_ord])
    # ensure that error-catching works here
    try({fit=fc_fit(time = boot_dat$time,model = mod_obj$mod_choice,censorID=NULL,rc.value=NULL,,eval(noHess))
          par_mat=rbind(par_mat,fit$par_tab[,1])})
  }

  # Evaluate initial sim
  tol_val=nrow(par_mat)/counter
  if(tol_val<tol){stop(paste="proportion of fitted models was :",tol_val,"; below the tolerance level: ",tol)}
  time_calc=(((proc.time()-bt)/50)*(nrep-50))[3]
  t_vals=c(time_calc,time_calc/60,time_calc/3600)
  print(t_vals)
  ind=length(which(c(t_vals>60,t_vals>3600)))
  t_est=paste(round(t_vals[ind+1]),c("seconds","minutes","hours")[ind+1])
  print(t_est)

  if(time_calc>40){message("based on the initial 50 samples the expected time before completion is: ",t_est)}

  # Remaining replications
  rem_reps=nrep-50; if(rem_reps<=950){message("Low sample size, consider increasing to at least 1000 for more reliable estimates")}
  while(nrow(par_mat) <= nrep){
    ind=sample(1:nrow(obs_dat),replace = T)
    t_raw=obs_dat[ind,1]

    t_ord=order(obs_dat[ind,1]) # sorting times
    t_new=fc_surv(t_raw[t_ord])
    boot_dat=data.frame(ind=ind[t_ord],
                        time=t_raw[t_ord],
                        surv_frac=t_new,
                        non_cen=(obs_dat$non_cen[ind])[t_ord])
    try({fit=fc_fit(time = boot_dat$time,model = mod_obj$mod_choice)
      par_mat=rbind(par_mat,fit$par_tab[,1])})
  }

  out=par_mat
  if(type=="pred"){
    out=apply(par_mat,1,function(x) fc_pred(pars = x,times = times,model=mod_obj$mod_choice))}

  return(out)
  }


library(devtools)

# common functions
document()

bad_time_dat=c(rep(5,10),rep(10,10))
fc_fit(bad_time_dat,model = "weibull")

#throws serious error
fc_fit(bad_time_dat,model = "all")

names(failCompare:::fc_mod_ls)


sapply(names(failCompare:::fc_mod_ls),function(x){fc_fit(bad_time_dat,model=x)})

modnms=names(failCompare:::fc_mod_ls)
for(i in 1:9){
  fc_fit(bad_time_dat,model = modnms[[i]])}

sapply( ?function(x){fc_fit(bad_time_dat,model=x)})






deps=c("vitality","flexsurv","survival","failcompare")
install.packages(deps)


Be sure that you have installed the 3 packages on which failCompare depends.
Paste the following code into the Console before installing failCompare 
# <return> [do not include this line!!]
###(Begin code block in Lucida Console if possible())[do not include this line!!]
deps=c("vitality","flexsurv","survival","failCompare")
install.packages(deps)



# testing functions

usethis::use_testthat()



sapply(deps,function(x){rownames(installed.packages(x))})
sapply(deps,function(x){x %in% rownames(installed.packages())})

search()

# sapply(c("vitality","flexsurv","survival"),remove.packages)


build()
sapply(c("vitality","flexsurv","survival"),use_package)

setwd("C:/Users/17072/OneDrive - UW/Desktop/Rdev")
utils::install.packages(pkgs = "failCompare_0.9.0.zip",repos=NULL,dependencies=T,
                        type="binary",destdir = getwd())


deps=c("vitality","flexsurv","survival")
sapply(deps,function(x){install.packages(x)})

library(failCompare)





data(sockeye)
taglife=sockeye[,"days"]
f_all=fc_fit(time=taglife,model="all")
modnms=c('weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma','vitality.ku','vitality.4p')


### Load example dataframe
data(sockeye)
taglife=sockeye[,"days"] #define vector of times

### Fit a 2-parameter Weibull model
weib_mod=fc_fit(time=taglife,model="weibull")

### Fit a 4-parameter Vitality 2013 model
vit_mod=fc_fit(time=taglife,model="vitality.4p")

# Combine two "fc_obj" objects into a model list of class "fc_list" 
fc_combine(mod_ls = list(weib_mod,vit_mod))

list(weib_mod,vit_mod)[[1]][["mod_choice"]]
list(weib_mod,vit_mod)[[2]][["mod_choice"]]

aa=as.character(sapply(list(weib_mod,vit_mod),function(x){x$mod_choice}))
bb=match(aa,names(fc_mod_ls))
order(bb)

unlist(sapply(c("weibull"),function(x) as.character(fc_mod_ls[[x]]),simplify = F,USE.NAMES = F))


unlist(sapply(c("weibull","vitality.ku"),function(x) as.character(fc_mod_ls[[x]]),simplify = F,USE.NAMES = F))

modnms[mod_ord]

unlist(sapply(out_ls[["mod_choice"]],function(x) as.character(fc_mod_ls[[x]]),simplify = F,USE.NAMES = F))


f_indiv=sapply(modnms,fc_fit,time=taglife,simplify = F)
f_indiv=sapply(rev(modnms),fc_fit,time=taglife,simplify = F)
tmp=fc_combine(f_indiv)


names(f_all$mod_choice)
names(tmp$mod_choice)

expect_identical(f_all,tmp)


expect_identical(f_all$fit_vals,tmp$fit_vals)

f_all_alt=tmp
rownames(f_all_alt$fit_vals)=1:nrow(f_all_alt[["fit_vals"]])
expect_identical(f_all$fit_vals,f_all_alt$fit_vals)

expect_identical(head(f_all_alt$fit_vals),head(f_all$fit_vals))


head(weib_mod$fit_vals)


expect_identical(f_all$mod_choice,tmp$mod_choice)


f_all$par_tab

rownames(f_all$fit_vals)
rownames(tmp$fit_vals$model)


expect_identical(f_all$mod_choice,as.character(sapply(f_indiv,function(x){x$mod_choice})))

ls_nms=as.character(sapply(f_indiv,function(x){x$mod_choice}))
mod_ord=match(modnms,ls_nms)
modnms[mod_ord]

names(f_indiv)
f_indiv[mod_ord]

expect_identical(mod_ord,order(mod_ord))


modnms

list(
  "weibull"=c("shape","scale"),
  "weibull3"=c("shape","thrsh","scale"),
  "gompertz"=c("shape","rate"),
  "gamma"=c("shape","rate"),
  "lognormal"=c("meanlog","sdlog"),
  "llogis"=c("shape","scale"),
  "gengamma"=c("mu","sigma","Q"),
  "vitality.ku"=c("r","s","k","u"),
  "vitality.4p"=c("r","s","lambda","beta"))

expect_


#'#' lognormal parameter names 
#' get_param_nm(mod_nm=c("lognormal"))
#' 
#' # 3-parameter Weibull and Gompertz parameter names 
#' get_param_nm(mod_nm=c("weibull3","gompertz"))
#' 
#' 
#' fc_mod_ls <- list(
"weibull"=c("shape","scale"),
"weibull3"=c("shape","thrsh","scale"),
"gompertz"=c("shape","rate"),
"gamma"=c("shape","rate"),
"lognormal"=c("meanlog","sdlog"),
"llogis"=c("shape","scale"),
"gengamma"=c("mu","sigma","Q"),
"vitality.ku"=c("r","s","k","u"),
"vitality.4p"=c("r","s","lambda","beta"))





#' @examples
#' # lognormal parameter names 
#' get_param_nm(mod_nm=c("lognormal"))
#' 
#' # 3-parameter Weibull and Gompertz parameter names 
#' get_param_nm(mod_nm=c("weibull3","gompertz")) 






eval(quote(failCompare::fc_fit),parent.frame())






data(sockeye)
taglife=sockeye[,"days"]
vit_mod=fc_fit(taglife,model="vitality.ku")
fc_boot(vit_mod,nrep = 100,times=1:50)

fc_boot(vit_mod,nrep = 100,times=seq(1,50,0.1))



# produces NA
vit_mod=fc_fit(taglife[seq(1,length(taglife),2)],model="vitality.ku")

vit_mod=fc_fit(taglife[seq(1,length(taglife),4)],model="vitality.ku")

vit_mod=fc_fit(taglife[1:10],model="vitality.ku")

try(vit_mod=fc_fit(taglife[30:40],model="vitality.ku"))


tryCatch(fc_fit(taglife[30:40],model="vitality.ku"),
         error = function(c) {
           msg <- conditionMessage(c)
           # invisible(
           # structure(msg, class = "try-error")
           # )
           msg
         })



# the following message
tryCatch(fc_fit(taglife[30:40],model="vitality.ku"),
         error = function(a) {
           msg <- conditionMessage(a)
           msg
         })


# the following message
ha=tryCatch(fc_fit(taglife[30:40],model="all"),
            error = function(a) {
              msg <- conditionMessage(a)
              msg
            },expr = message("model could not be fit"))

# calling handlers
ha=withCallingHandlers(fc_fit(taglife[30:40],model="all"),
                       error = function(a) {
                         msg <- conditionMessage(a)
                         msg
                       })




try2 <- function(code, silent = FALSE) {
  tryCatch(code, error = function(c) {
    msg <- conditionMessage(c)
    if (!silent) message(c)
    invisible(structure(msg, class = "try-error"))
  })
}


mod_ls=fc_fit(taglife,model="all")
mod_ls$par_tab

fc_boot(mod_obj=weib_mod,times = 1:20,nrep = 10)
# library(testthat)
