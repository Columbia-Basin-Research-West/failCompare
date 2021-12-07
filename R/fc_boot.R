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
#' fc_boot(weib_mod,nrep=60,times = 10:20)
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

    try({fit=fc_fit(time = boot_dat$time,model = mod_obj$mod_choice,SEs=F,...)
          par_mat=rbind(par_mat,fit$par_tab[,1])})
  }

  # Evaluate initial sim
  tol_val=nrow(par_mat)/counter
  if(tol_val<tol){stop(paste="proportion of fitted models was :",tol_val,"; below the tolerance level: ",tol)}
  time_calc=(((proc.time()-bt)/50)*(nrep-50))[3]
  t_vals=c(time_calc,time_calc/60,time_calc/3600)
  ind=length(which(c(t_vals>60,t_vals>3600)))
  t_est=paste(round(t_vals[ind+1]),c("seconds","minutes","hours")[ind+1])

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
    try({fit=fc_fit(time = boot_dat$time,model = mod_obj$mod_choice,...)
      par_mat=rbind(par_mat,fit$par_tab[,1])})
  }

  out=par_mat
  if(type=="pred"){
    out=t(apply(par_mat,1,function(x) fc_pred(pars = x,times = times,model=mod_obj$mod_choice)))
    colnames(out)=as.character(times)
    }

  return(out)
  }

