#' @title Failure time predictions
#' @description Geneates predictions from failure time model objects
#'
#' @param times time vector
#' @param pars parameter estimates
#' @param model survival model
#'
#' @return failure probability
#'
#' @importFrom survival Surv
#' @importFrom flexsurv pgompertz
#' @importFrom vitality vitality.ku
#' @importFrom vitality vitality.4p
#' @export fc_pred
#'
fc_pred <- function(times,pars,model="gompertz"){
  if(model=="weibull"){
    tmp=1-stats::pweibull(times,shape = pars[1],scale=pars[2])
  }
  if(model=="weibull3"){
    # from Rich's code
    tmp=exp(-((times-pars[2])/pars[3])^pars[1])
    tmp[is.na(tmp)]=1
    # tmp=1-stats::pweibull(times-pars[2],shape = pars[3],scale=pars[2]) # likely equivalent
  }
  if(model=="gompertz"){
    tmp=1-pgompertz(times,shape = pars[1],rate=pars[2])
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
  return(tmp)
}
