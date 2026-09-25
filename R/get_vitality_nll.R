#' @title Recover the negative log-likelihood from a vitality model fit
#'
#' @description \code{vitality::vitality.ku()} and \code{vitality::vitality.4p()} print
#' their minimum negative log-likelihood value when \code{silent = FALSE}, but do not
#' return it as part of their output. This helper recomputes that value from the fitted
#' parameters using the (non-exported) \code{vitality:::dataPrep()} and
#' \code{vitality:::logLikelihood.ku()}/\code{vitality:::logLikelihood.4p()} functions,
#' which is exactly what \code{vitality.ku()}/\code{vitality.4p()} use internally.
#'
#' @param model either \code{"vitality.ku"} or \code{"vitality.4p"}
#' @param time failure-only times, sorted, as passed to the vitality model's \code{time} argument
#' @param sdata survival fraction corresponding to \code{time}, as passed to \code{sdata}
#' @param rc.data logical right-censoring flag, as passed to \code{rc.data}
#' @param params fitted parameter vector, in the order returned by the vitality model
#'
#' @return numeric negative log-likelihood value at \code{params}
#'
#' 
#' @keywords internal
get_vitality_nll=function(model,time,sdata,rc.data,params){
  dTmp=vitality:::dataPrep(time=time,sdata=sdata,datatype="CUM",rc.data=rc.data)
  x1=dTmp$x1; x2=dTmp$x2; Ni=dTmp$Ni
  # vitality.4p() rotates/drops the first interval internally when time[1] > 0
  if(model=="vitality.4p" && time[1]>0){
    x1=c(x1[-c(1,length(x1))],x1[1])
    x2=c(x2[-c(1,length(x2))],0)
    Ni=Ni[-1]
  }
  ll_fun=if(model=="vitality.ku"){vitality:::logLikelihood.ku}else{vitality:::logLikelihood.4p}
  ll_fun(params,x1,x2,Ni)
}
