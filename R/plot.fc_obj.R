
# function that plots fitted values when a failmod_obj is called
#' Title
#'
#' @param x failmod_obj
#' @param type Plotting survival curve of data ("data") versus difference between Kaplan-Meier estimates and predictions from a parametric model ("resid")
#' @param km Show kaplan-meier estimates
#' @param km.ci Show 95% confidence limits surrounding kaplan-meier estimates
#' @param res Number of evenly space points within the range of the data
#' @param ... ignored
#'
#'
#' @return plot of taglife
#'
#' @importFrom survival Surv
#' @importFrom graphics plot lines legend points
#'
#' @export
plot.fc_obj <- function(x,km=FALSE,km.ci=FALSE,res=100,type="data",...){
  stopifnot(is.logical(km))
  stopifnot(is.logical(km.ci))
  stopifnot(type %in% c("data","resid"))

  if(x$mod_choice=="kaplan-meier"){
    if(type=="resid"){stop("residual plot not available for kaplan-meier model")}
    t_rng=x$times$time
    ts=seq(max(min(t_rng*.95),0),(max(t_rng)*1.05),length.out = res)
    plot(surv_frac~time,x$times,pch=3,col=NA,xlab="t",ylab=expression(hat("F")(t)),xlim=c(min(ts),max(ts)),...)
    lines(est~time,x$KM_DF,type="s",col=8,lty=2)
    legend(legend=c(x$mod_choice,"kaplan-meier (Est)"),"bottomleft",col=c(2,8),lwd=c(4,1),lty=c(1,2))
    if(km.ci){
      lines(lcl~time,x$KM_DF,type="s",col=8,lty=3)
      lines(ucl~time,x$KM_DF,type="s",col=8,lty=3)
      legend(legend=c(x$mod_choice,"kaplan-meier (Est)","kaplan-meier (95% CI)"),"bottomleft",col=c(2,8,8),lwd=c(4,1,1),lty=c(1,2,3))
    }
      points(surv_frac~time,x$times,pch=3,col=1)
    # lines(est~time,x$KM_DF,type="s",col=8,lty=2)
    # lines(lcl~time,x$KM_DF,type="s",col=8,lty=3)
    # lines(ucl~time,x$KM_DF,type="s",col=8,lty=3)
    # points(surv_frac~time,x$times,pch=3,...)
    # legend(legend=c("kaplan-meier (Est)","kaplan-meier (95% CI)"),"bottomleft",col=8,lty=c(2,3))
  }
  else{  # Any model other than KM
  # time increment def.
  t_rng=x$fit_vals$time
  ts=seq(max(min(t_rng*.95),0),(max(t_rng)*1.05),length.out = res)
  # survival preds.
  spred=fc_pred(times=ts,model=x$mod_choice,pars = x$par_tab[,1])

  # Data plot
  if(type=="data"){
  plot(surv_frac~time,x$times,pch=3,col=NA,xlab="t",ylab=expression(hat("F")(t)),xlim=c(min(ts),max(ts)),...) #empty plot
  legend(legend=x$mod_choice,"bottomleft",col=2,lwd=4) # legend covered up if km==T
  lines(ts,spred,col=2,lwd=4)
  if(km){
    lines(est~time,x$KM_DF,type="s",col=8,lty=2)
    legend(legend=c(x$mod_choice,"kaplan-meier (Est)"),"bottomleft",col=c(2,8),lwd=c(4,1),lty=c(1,2))
    if(km.ci){
      lines(lcl~time,x$KM_DF,type="s",col=8,lty=3)
      lines(ucl~time,x$KM_DF,type="s",col=8,lty=3)
      legend(legend=c(x$mod_choice,"kaplan-meier (Est)","kaplan-meier (95% CI)"),"bottomleft",col=c(2,8,8),lwd=c(4,1,1),lty=c(1,2,3))
    }}
  points(surv_frac~time,x$times,pch=3,col=1,xlab="t",ylab=expression(hat("F")(t)))
  }
  # Residual plot
  if(type=="resid"){
  est_offset=x$fit_vals$est[match(x$KM_DF$time,x$fit_vals$time)]
  est_lcl=x$KM_DF$lcl-est_offset
  est_ucl=x$KM_DF$ucl-est_offset
  tmpDF=data.frame(x$KM_DF$time,est_offset,x$KM_DF$lcl,x$KM_DF$ucl,est_lcl,est_ucl)
  bnds=unlist(tmpDF[,c("est_lcl","est_ucl")])
  min_bnds=min(bnds);max_bnds=max(bnds)
  plot(x$fit_vals$time,c(1,x$times$surv_frac)-x$fit_vals$est,col=NA,ylim=c(min_bnds,max_bnds),xlim=c(min(ts),max(ts)),xlab="t",ylab="Residual(KM-Fitted)",...)
  lines(y=rep(0,length(ts)),x=ts,col=2,lwd=4)
  points(x$fit_vals$time,c(1,x$times$surv_frac)-x$fit_vals$est,pch=3)
  legend(legend=c(x$mod_choice,"kaplan-meier (Est)"),"bottomleft",col=c(2,8),lwd=c(4,1),lty=c(1,2))
  lines(x$fit_vals$time,c(1,x$times$surv_frac)-x$fit_vals$est,lty=2,col=8)
  if(km.ci){
    lines(x=x$KM_DF$time,est_lcl,col=8,lty=3)
    lines(x=x$KM_DF$time,est_ucl,col=8,lty=3)
    legend(legend=c(x$mod_choice,"kaplan-meier (Est)","kaplan-meier (95% CI)"),"bottomleft",col=c(2,8,8),lwd=c(4,1,1),lty=c(1,2,3))
    }
  }
  }
}

