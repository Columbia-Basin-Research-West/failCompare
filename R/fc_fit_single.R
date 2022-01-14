#' @title Fitting a single failure time model 
#'
#' @details function for fitting an indiviudal failure time model assuming inputs have been vetted by user-facing function fc_fit().
#'
#' @param y failure time
#' @param model failure time model
#' @param Hess calculating standard errors
#' @param non_cen logical of length(y)
#' @param y_sfrac survival fraction
#' @param KM_DF K-M model predictions
#' @param KM_mod K-M model object
#'
#' @return "fc_obj" if successful NULL if otherwise
#'
fc_fit_single=function(y,y_sfrac,model,Hess,non_cen,KM_DF,KM_mod){
  KM_mod=get("KM_mod",inherits = T)
  rc=ifelse(all(non_cen),FALSE,TRUE)
    # FITTING DISTRIBUTIONS IN THE FLEXSURV PACKAGE
    if(model %in% names(fc_mod_ls)[names(fc_mod_ls) %in% names(flexsurv::flexsurv.dists)]){
      flex_mod=quote(model)
      q_e=quote(flexsurv::flexsurvreg(survival::Surv(time=y,event=non_cen) ~ 1,
                                      dist = model,hessian = eval(Hess)))
      fit <- fc_tryfit(y = y,non_cen = non_cen,fit_call = q_e,model = model,Hess=Hess)
      preds=summary(fit)[[1]]
      # fit_vals=rbind(fit_vals,
      fit_vals=rbind(data.frame(model=model,time=0,est=1,lcl=1,ucl=1),
                     data.frame(model=model,preds[match(y,preds$time),]))
                     #) # match() used for duplicate rows
    }
    else{
      if(model=="vitality.ku"){
        # Defines function call depending on right censoring or not
        if(rc){
          q_e=quote(vitality::vitality.ku(time=sort(y),sdata = y_sfrac,rc.data = T,pplot =F,silent=T,se=Hess))
        }
        else{
          q_e=quote(vitality::vitality.ku(time = sort(y),sdata = y_sfrac,pplot =F,silent = T,se=Hess,rc.data=F))
        }
        fit = fc_tryfit(fit_call = q_e,y = y,y_sfrac=y_sfrac,model="vitality.ku",Hess=Hess)
        pars_tmp=fit
        fit_vals=data.frame(model="vitality.ku",
                                  time=c(0,y),
                                  est=c(1,vitality::SurvFn.ku(y,pars_tmp[1],pars_tmp[2],pars_tmp[3],pars_tmp[4])),
                                  lcl=0,ucl=0)
      }
      
      if(model=="vitality.4p"){
        # Defines function call depending on right censoring or not
        if(rc){
          q_e=quote(vitality::vitality.4p(time = y,sdata = y_sfrac,pplot =F,silent = T,se = Hess))
        }
        else{
          q_e=quote(vitality::vitality.4p(time = y,sdata = y_sfrac,pplot =F,silent = T,se = Hess))
          }
        fit = fc_tryfit(fit_call = q_e,y = y,y_sfrac=y_sfrac,model="vitality.4p",Hess=Hess)
        pars_tmp=fit
        fit_vals=data.frame(model="vitality.4p",
                                  time=c(0,y),
                                  est=c(1,vitality::SurvFn.4p(y,pars_tmp[1],pars_tmp[2],pars_tmp[3],pars_tmp[4])),
                                  lcl=0,ucl=0)
      }
      if(model=="weibull3"){

        q_e=quote(taglife.fn_weib3(y,model.in = "weibull",tag.se=eval(Hess)))
        tmp=fc_tryfit(y=y,fit_call=q_e,model="weibull3",Hess = Hess)
        if(is.vector(tmp$par_tab)){Hess=F} # switch to FALSE if hessian caused error
        fit = tmp$mod_obj
        pars_tmp=tmp$par_tab
        fit_vals=tmp$fit_vals
      }
    }
  
  # if a non-flexsurv model
  if(model=="vitality.ku" | model=="vitality.4p" | model =="weibull3"){
      if(Hess){par_tab=fit[,c("params","std")]}
      else{par_tab=matrix(c(fit,rep(NA,length(fc_mod_ls[[model]]))),ncol=2)}
      pnms=fc_mod_ls[[model]] # number parameters by default
    }
    else{
      par_tab=fit$res[,c("est","se")]
    }
  mod=list("mod_choice"=model,
              "times"=data.frame(time=y,surv_frac=y_sfrac,non_cen=non_cen),
              "fit_vals"=fit_vals,
              "mod_objs"=fit,
              "par_tab"=par_tab,
              "KM_DF"=KM_DF,
              "KM_mod"=KM_mod,
              "censored"=rc)
  out=structure(mod,class="fc_obj")
  
  return(out)
}
