#' @title Failure time predictions
#' @description This generates predictions from failure time model objects.
#'
#' @param times time vector 
#' @param mod_obj model object (class = fc_obj)
#' @param pars parameter estimates, if mod_obj absent
#' @param model survival model name, if mod_obj absent
#'
#' @return numeric vector failure/survival probability
#' @import flexsurv survival
#'
#' @export fc_pred
#'
fc_pred <- function(mod_obj=NULL,times,pars=NULL,model=NULL){
  # stopifnot(is.numeric(times)) # requires times argument
  if(!is.null(mod_obj)){
    stopifnot((class(mod_obj)=="fc_obj")) #checking that it sit a fc_obj
    if(!is.null(model) | !is.null(pars)){warning("'pars' and 'model' arguments overridden by 'mod_obj'")} # notifies user that other arguments are overridden
      #accessing model parameters and choices
      pars=mod_obj$par_tab[,1]
      model=mod_obj$mod_choice}

    if(model=="weibull"){
      tmp=1-stats::pweibull(times,shape = pars[1],scale=pars[2])
    }
    if(model=="weibull3"){
      tmp=exp(-((times-pars[2])/pars[3])^pars[1])
      tmp[is.na(tmp)]=1
    }
    if(model=="gompertz"){
      tmp=1-flexsurv::pgompertz(times,shape = pars[1],rate=pars[2])
    }
    if(model=="gamma"){
      tmp=1-stats::pgamma(times,shape = pars[1],rate=pars[2])
    }
    if(model=="lognormal"){
      tmp=1-stats::plnorm(times,meanlog = pars[1],sdlog=pars[2])
    }
    if(model=="llogis"){
      tmp=1-flexsurv::pllogis(times,shape = pars[1],scale=pars[2])
    }
    if(model=="gengamma"){
      tmp=1-flexsurv::pgengamma(times,mu = pars[1],sigma=pars[2],Q=pars[3])
    }
    if(model=="vitality.ku"){
      tmp=vitality::SurvFn.ku(times,r=pars[1],s=pars[2],k=pars[3],u=pars[4])
    }
    if(model=="vitality.4p"){
      tmp=vitality::SurvFn.4p(times,r=pars[1],s=pars[2],lambda=pars[3],beta=pars[4])
    }
    if(model=="kaplan-meier"){
    tmp=sapply(times,function(x){mod_obj$KM_DF$est[as.numeric(cut(x,mod_obj$KM_DF$time))]})
    }

  return(tmp)
}

