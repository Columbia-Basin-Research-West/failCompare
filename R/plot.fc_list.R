# function that plots fitted values when a failmod_list is called
#' Title
#'
#' @param x failmod_list objects (ranked or not)
#' @param model surival model
#' @param res fineness of survival function preds
#' @param km Show kaplan-meier estimates
#' @param ... ignore
#'
#' @return plot and a message
#' @export
plot.fc_list <- function(x,model=NULL,km=F,res=100,...){
  # validation
  stopifnot(all(model %in% c("weibull", "gompertz", "gamma", "lognormal", "llogis", "gengamma","vitality.ku","vitality.4p","weibull3")))
  # time increment def.
  t_rng=x$fit_vals$time
  ts=seq(max(min(t_rng*.95),0),(max(t_rng)*1.05),length.out = res)

  # if fail_rank() has not been run yet
  if(is.null(x$GOF_tab)){
    # if the model argument is empty
    if(is.null(model)){
      mod_plts=x$mod_choice}
    # model arguement
    else{mod_plts=model}
    # if >3 total models
    if(length(mod_plts)>3){
      mod_plts=mod_plts[1:3]
      message(paste("First 3 models in the list shown\n\nAdditional models in the list:",paste(x$mod_choice[!x$mod_choice %in% mod_plts],collapse=";")))
    }
  }
  else{  # if ranked
    mod_plts=x$GOF_tab$model
    # if the model argument is empty
    if(is.null(model)){
      mod_plts=x$mod_choice
    }
    # model arguement
    else{mod_plts=model
    }
    if(length(mod_plts)>3){
      mod_plts=x$GOF_tab[,"model"][1:3]
      unused_mod=x$mod_choice[!x$mod_choice %in% mod_plts]
      rnk=which(as.character(x$GOF_tab[,"model"]) %in% unused_mod)
      unused_mod=paste(unused_mod,"(",rnk,")",sep = "")
      message(paste("\nAdditional models with rankings:",paste(unused_mod,collapse="; ")))
    }
    # not more than 3
    else{
      unused_mod=x$mod_choice[!x$mod_choice %in% mod_plts]
      rnk=which(as.character(x$GOF_tab[,"model"]) %in% unused_mod)
      unused_mod=paste(unused_mod,"(",rnk,")",sep = "")
      message(paste("\nAdditional models with rankings:",paste(unused_mod,collapse="; ")))
    }
  }

  # plotting k-m
  plot(surv_frac~time,x$times,pch=3,col=NA,xlab="t",ylab=S(t),xlim=c(min(ts),max(ts)),...) #empty plot

  spred=list()
  for(i in 1:length(mod_plts)){
    # survival preds.
    tmp=x$par_tab[x$par_tab$model==mod_plts[i],"est"]
    spred[[i]]=fc_pred(times=ts,model=mod_plts[i],pars = tmp)
    lines(ts,spred[[i]],col=i+1,lwd=3,lty=i)
  }

  #plotting survival functions
  if(!is.null(x$GOF_tab)){
    if(is.null(model)){
      legend(legend=paste(mod_plts," (",1:length(mod_plts),")",sep=""),"bottomleft",col=(1:length(mod_plts))+1,lwd=3,lty=1:length(mod_plts),title = "Ranked models",bty = "n")
      if(km){
        lines(est~time,x$KM_DF,type="s",col=8,lty=2)
        legend(legend=c("kaplan-meier (Est)",paste(mod_plts," (",1:length(mod_plts),")",sep="")),"bottomleft",col=c(8,(1:length(mod_plts))+1),lwd=3,lty=1:length(mod_plts),title = "Ranked models",bty = "n")}
    }
    else{
      sel_rnk=which(as.character(x$GOF_tab[,"model"]) %in% mod_plts)
      legend(legend=paste(mod_plts," (",sel_rnk,")",sep=""),"bottomleft",col=(1:length(mod_plts))+1,lwd=3,lty=1:length(mod_plts),title = "Ranked models",bty = "n")
    }
  }
  else{
    legend(legend=mod_plts,"bottomleft",col=(1:length(mod_plts))+1,lwd=3,lty=1:length(mod_plts),title = "Models",bty = "n")
    if(km){
      lines(est~time,x$KM_DF,type="s",col=8,lty=2)
      legend(legend=c("kaplan-meier (Est)",mod_plts),"bottomleft",col=c(8,(1:length(mod_plts))+1),lwd=3,lty=1:length(mod_plts),title = "Models",bty = "n")}
  }
  points(surv_frac~time,x$times,pch=3)
}
