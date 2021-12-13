#' @title Fitting a single failure time model 
#'
#' @details function for fitting an indiviudal failure time model assuming inputs have been vetted by primary function fc_fit().
#'
#' @param time failure time
#' @param model failure time model
#' @param SEs calculating standard errors
#' @param censorID logical whether an observation is non-censored
#' @param rc.value value after which
#' @param ... additional arguments passed to fc_tryfit()
#'
#' @return 
#'
fc_fit_single=function(time,model,SEs=TRUE,censorID=NULL,rc.value=NULL,...){
  
  # non_cen
  # fit=list()
  fit_vals=NULL
  for (i in 1:length(model)){
    # FITTING DISTRIBUTIONS IN THE FLEXSURV PACKAGE
    if(model %in% names(fc_mod_ls)[names(fc_mod_ls) %in% names(flexsurv::flexsurv.dists)]){
      flex_mod=quote(model)
      q_e=quote(flexsurv::flexsurvreg(survival::Surv(time=y,event=non_cen) ~ 1,
                                      dist = model,hessian = Hess))
      fit[[model]] <- fc_tryfit(y = y,non_cen = non_cen,fit_call = q_e,model = model,Hess=Hess)
      preds=summary(fit[[i]])[[1]]
      fit_vals=rbind(fit_vals,
                     data.frame(model=model,time=0,est=1,lcl=1,ucl=1),
                     data.frame(model=model,preds[match(y,preds$time),])) # match() used for duplicate rows
    }
    else{
      if(model=="vitality.ku"){
        # Defines function call depending on right censoring or not
        if(rc){
          dTmp=vitality::dataPrep(c(0,y_cen),(n_cen:(n_cen-length(y_cen)))/n_cen,datatype="CUM",rc.data=(n_cen>length(y_cen)))
          q_e=vitality::vitality.ku(dTmp[,"time"],sdata = dTmp[,"sfract"],rc.data = T,pplot =F,silent=T,se=Hess,...)
        }
        else{q_e=vitality::vitality.4p(time = sort(y),sdata = y_sfrac,pplot =F,silent = T,se=Hess)}
        fit[[model]] = fc_tryfit(fit_call = q_e,y = y,model="vitality.ku")
        pars_tmp=fit[[model]]
        fit_vals=rbind(fit_vals,
                       data.frame(model="vitality.ku",
                                  time=c(0,y),
                                  est=c(1,vitality::SurvFn.ku(y,pars_tmp[1],pars_tmp[2],pars_tmp[3],pars_tmp[4])),
                                  lcl=0,ucl=0))}

      if(model=="vitality.4p"){
        # Defines function call depending on right censoring or not
        if(rc){
          dTmp=vitality::dataPrep(c(0,y_cen),(n_cen:(n_cen-length(y_cen)))/n_cen,datatype="CUM",rc.data=(n_cen>length(y_cen)))
          q_e=vitality::vitality.4p(dTmp[,"time"],sdata = dTmp[,"sfract"],rc.data = T,pplot =F,silent=T,se=Hess,...)
        }
        else{q_e=vitality::vitality.4p(time = sort(y),sdata = y_sfrac,pplot =F,silent = T,se=Hess)}
        fit[[model]] = fc_tryfit(fit_call = q_e,y = y,model="vitality.4p")
        pars_tmp=fit[[model]]
        fit_vals=rbind(fit_vals,
                       data.frame(model="vitality.4p",
                                  time=c(0,y),
                                  est=c(1,vitality::SurvFn.4p(y,pars_tmp[1],pars_tmp[2],pars_tmp[3],pars_tmp[4])),
                                  lcl=0,ucl=0))
      }
      if(model=="weibull3"){
        y=sort(time)
        q_e=quote(taglife.fn_weib3(y,model.in = "weibull",tag.se=eval(Hess)))
        tmp=fc_tryfit(y=y,fit_call=q_e,model="weibull3",Hess = eval(Hess))
        if(is.vector(tmp$par_tab)){Hess=F} # switch to FALSE if hessian caused error
        fit[[model]] = tmp$mod_obj
        pars_tmp=tmp$par_tab
        fit_vals=rbind(fit_vals,tmp$fit_vals)
      }
    }
  }
  # if only K-M specified
  if(model=="kaplan-meier"){
    par_tab=fit_vals=NULL
    fit[[1]]=KM_mod
    fit_vals=KM_DF
  }
  # Otherwise
  else{
    # if a non-flexsurv model
    if(model=="vitality.ku" | model=="vitality.4p" | model =="weibull3"){
      if(Hess){par_tab=fit[[1]][,c("params","std")]}
      else{par_tab=matrix(c(fit[[1]],rep(NA,length(fc_mod_ls[[model]]))),ncol=2)}
      pnms=fc_mod_ls[[model]] # number parameters by default
    }
    else{
      par_tab=fit[[1]]$res[,c("est","se")]
    }
  }
  out_ls=list("mod_choice"=model,
              "times"=data.frame(time=y,surv_frac=y_sfrac,non_cen=non_cen),
              "fit_vals"=fit_vals,
              "mod_objs"=fit[[1]],
              "par_tab"=par_tab,
              "KM_DF"=KM_DF,
              "KM_mod"=KM_mod,
              "censored"=rc)
  rownames(out_ls[["par_tab"]])=NULL
  out=structure(out_ls,class="fc_obj")
  
  return(out)
}
