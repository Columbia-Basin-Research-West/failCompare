#' @title Fitting one or a set of failure time models
#' @description Code for fitting a particular tag-life model or models
#'
#'
#' @param time numeric vector of failure times
#' @param model character vector of specified model(s)
#' @param rc.value right-censoring cutoff value; only observations with times > rc.value are considered to have failed prematurely.
#' @param censorID censored observation vector (T/F)
#' @param ... additional arguments passed
#'
#' @return Returns failure model object of class \code{"fc_obj"} if one model specified OR 
#' a failure model list object of class \code{"fc_list"} if multiple models are specified
#'
#' @details
#' Model fitting routine used to fit one or a set of failure time models.
#' 
#' The following nine model options include
#' `weibull` = 2-parameter Weibull 
#' `weibull3` = 3-parameter Weibull
#' `gompertz` = Gompertz
#' `lognormal` = lognormal
#' `llogis` = log-logistic
#' `gamma` = gamma (2-parameter)
#' `gengamma` = generalized gamma (3-parameter)
#' `vitality.ku` = vitality (2009)
#' `vitality.4p` = vitality (2013)
#' `kaplan-meier` = Kaplan-Meier nonparametric estimate (NOTE: this model cannot be specified in a list with any other model)
#'
#' If a single model is specified a "fc_obj" is created, which can be
#' used to adjust a CJS model in the "ATLAS" package.
#'
#' If multiple models are specified a "fc_list" is created containing
#' output from all model fits. This object serves as an input in the
#' fc_rank() function, which ranks the performance of the model using
#' the \href{http://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-020-00213-z}{Skalski and Whitlock (2020)} GOF measure.
#'
#'   Parameter estimates and accompanying standard errors(if available) are printed. 
#'   The \code{"fc_obJ"} is a list of the following extractable objects:
#'  \itemize{
#'     \item"mod_choice"  = official model name
#'     \item     "times"  = dataframe of failure time, survival fraction, and censoring binary var.
#'     \item  "fit_vals"  = failure times and predicted survival under the model, 95% LCL an UCL if available
#'     \item  "mod_objs"  = actual model object created by "flexsurvdist" or "vitality package"-- much more to extract from "flexsurvdist
#'     \item   "par_tab"  = table of parameter estimates and SE in failCompare recognized order
#'     \item     "KM_DF"  = table of K-M estimates for plotting
#'     \item    "KM_mod"  = survival package K-M model estimates
#'     \item  'censored'  = T/F of censored dataset or not
#'  }
#'
#' @references Skalski, J. R., and S. L. Whitlock. 2020. Vitality models found useful in modeling tag-failure times in acoustic-tag survival studies. Animal Biotelemetry 8:1–10.DOI:10.1186/s40317-020-00213-z
#'
#'
#' @importFrom survival Surv
#' @importFrom flexsurv flexsurvreg
#' @importFrom vitality vitality.ku
#' @importFrom vitality vitality.4p
#' @export fc_fit 
#'
#'
fc_fit=function(time,model,censorID=NULL,rc.value=NULL,...){
  rc=FALSE #temp def
  rt=FALSE #temp def
  if(!is.vector(time)|!is.numeric(time)){stop("A numeric vector is expected for the 'time' argument")}
  if(model[1]=="all"){model=c("weibull",'weibull3', "gompertz", "gamma", "lognormal", "llogis", "gengamma","vitality.ku","vitality.4p")
  message("Fitting all available parametric survival models")}

  # WARNINGS AND VALIDATION
  if(any(is.na(time))){
    message(paste(length(which(is.na(time))),"NA times removed"))}

  test=!model %in% c("weibull",'weibull3', "gompertz", "gamma", "lognormal", "llogis", "gengamma","vitality.ku","vitality.4p","kaplan-meier")
  if(any(test)){
    message(paste("Model names not recognized:",paste(model[which(test)],collapse=";")," \n Default model names = {'weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma','vitality.ku','vitality.4p','kaplan-meier'}",sep=""))
    model=model[!test]
    if(length(model)==0){stop("No valid names in the 'model' argument")}}

  ord=order(time)
  y=sort(time)  # sorted data necessary for Vitality package functions
  y_sfrac=sapply(y,function(x){1-length(which(y<=x))/length(y)}) # survival fraction calc
  non_cen=rep(TRUE,length(y))
  n_cen=length(y)

  if(!is.null(rc.value)){
    rc=TRUE # change this value for later if statement
    y=sort(y)  # sorted data necessary for Vitality package functions
    non_cen=ifelse(y<rc.value,TRUE,FALSE) # vector used by "flexsurv"
    y_sfrac=sapply(y,function(x){1-length(which(y<=x))/length(y)}) # survival fraction calc

    # For vitality model
    y_cen=y[y<rc.value]
    y_cen_sfrac=y_sfrac[y<rc.value]
  }
  
  # censorID
  if(!is.null(censorID)){
    rc=TRUE # change this value for later if statement
    censorID=censorID[ord]
    y_cen=y[censorID]
    stopifnot(length(time)==length(censorID)) # censorID length should match
    if(any(sapply(censorID,function(x){!(x %in% c(0,1) | is.logical(x))}))){stop("1/0 or TRUE/FALSE expected for censorID")}
    if(!is.null(rc.value)){warning("censorID overrides rc.value argument")}
    non_cen=as.logical(censorID)
  }
  

  KM_mod=survival::survfit(survival::Surv(time=y,event=non_cen)~1)
  KM_sls=summary(KM_mod)
  KM_DF=rbind(data.frame(model="kaplan-meier",time=0,est=1,lcl=1,ucl=1),
              data.frame(model="kaplan-meier",time=KM_sls$time,est=KM_sls$surv,lcl=KM_sls$lower,ucl=KM_sls$upper))
  KM_DF[nrow(KM_DF),c("lcl","ucl")]=c(0,0) # replacing last two NAs wih 0s

  if("kaplan-meier" %in% model & length(model)>1){
    model=model[model!="kaplan-meier"]
    message("'kaplan-meier' cannot be included in a fc_list with parametric models, and will be omitted from the string")
  }
  
  ### WEIBULL3 not supported
  if(rc & "weibull3" %in% model){
    model=model[model!="weibull3"]
    message("The right-censored 3-parameter Weibull model ('weibull3') is not available ")
    if(length(model)==0){stop("Cannot fit right-censored 3-parameter Weibull model")}
  }


  fit=list()
  fit_vals=NULL
  for (i in 1:length(model)){
    # FITTING DISTRIBUTIONS IN THE FLEXSURV PACKAGE
    if(model[i] %in% c("weibull", "gompertz", "gamma", "lognormal", "llogis", "gengamma")){
      fit[[model[i]]] <- flexsurv::flexsurvreg(survival::Surv(time=y,event=non_cen) ~ 1, dist = model[i])
      preds=summary(fit[[i]])[[1]]
      # preds[match(time,preds$time),]
      fit_vals=rbind(fit_vals,
                     data.frame(model=model[i],time=0,est=1,lcl=1,ucl=1),
                     data.frame(model=model[i],preds[match(y,preds$time),])) # match() used for duplicate rows
      }
    else{
      if(model[i]=="vitality.ku"){
        if(rc){
          dTmp=vitality::dataPrep(c(0,y_cen),(n_cen:(n_cen-length(y_cen)))/n_cen,datatype="CUM",rc.data=(n_cen>length(y_cen)))
          fit[[model[i]]]=vitality::vitality.ku(dTmp[,"time"],sdata = dTmp[,"sfract"],rc.data = T,pplot =F,lplot = F, silent=T,
                                        se=T,...)#init.params = c(0.0118535685, 0.0070384229, 0.0001868951, 0.0443095393),,...) # added ... here to get the init.params argument to pass through
          pars_tmp=fit[[model[i]]][,"params"]
          fit_vals=rbind(fit_vals,
                         data.frame(model="vitality.ku",
                                    time=c(0,y), # survival proportions
                                    est=c(1,vitality::SurvFn.ku(y,pars_tmp[1],pars_tmp[2],pars_tmp[3],pars_tmp[4])),
                                    lcl=0,ucl=0))
        }
        else{
        fit[[model[i]]] = vitality::vitality.ku(time = sort(y),sdata = y_sfrac,se=T,pplot =F,silent=T,lplot = F,...)
        pars_tmp=fit[[model[i]]][,"params"]
        fit_vals=rbind(fit_vals,
                       data.frame(model="vitality.ku",
                                  time=c(0,y), # survival proportions
                                  est=c(1,vitality::SurvFn.ku(y,pars_tmp[1],pars_tmp[2],pars_tmp[3],pars_tmp[4])),
                                  lcl=0,ucl=0))
        }
      }
      if(model[i]=="vitality.4p"){
        if(rc){
          dTmp=vitality::dataPrep(c(0,y_cen),(n_cen:(n_cen-length(y_cen)))/n_cen,datatype="CUM",rc.data=(n_cen>length(y_cen)))
          fit[[model[i]]]=vitality::vitality.4p(dTmp[,"time"],sdata = dTmp[,"sfract"],rc.data = T,pplot =F,silent=T,
                                                se=T,init.params=c(0.012, 0.01, 0.1, 0.1),...) # added ... here to get the init.params argument to pass through
          pars_tmp=fit[[model[i]]][,"params"]
          fit_vals=rbind(fit_vals,
                         data.frame(model="vitality.4p",
                                    time=c(0,y), # survival proportions
                                    est=c(1,vitality::SurvFn.4p(y,pars_tmp[1],pars_tmp[2],pars_tmp[3],pars_tmp[4])),
                                    lcl=0,ucl=0))
        }
        else{
          fit[[model[i]]] = vitality::vitality.4p(time = sort(y),sdata = y_sfrac,se=T,pplot =F,silent = T,init.params=c(0.012, 0.01, 0.1, 0.1),...)
        pars_tmp=fit[[model[i]]][,"params"]
        fit_vals=rbind(fit_vals,
                       data.frame(model="vitality.4p",
                                  time=c(0,y), # survival proportions
                                  est=c(1,vitality::SurvFn.4p(y,pars_tmp[1],pars_tmp[2],pars_tmp[3],pars_tmp[4])),
                                  lcl=0,ucl=0))
        }
      }
      if(model[i]=="weibull3"){
        tmp=taglife.fn_weib3(sort(y),model.in = "weibull")
        fit[[model[i]]] = tmp$mod_obj
        pars_tmp=tmp$par_tab
        fit_vals=rbind(fit_vals,tmp$fit_vals)
      }
    }
  }

  # Combined table of parameter estimates for all fitted models (Associated with fc_list)
  if(length(model)>1){
    par_ls=sapply(seq_along(fit),function(i){
      tmp_mod_nm=names(fit)[i]
      if(class(fit[[i]])[1]=="flexsurvreg"){tmp=fit[[i]]$res[,c("est","se")]}
      else{tmp=fit[[i]][,c("params","std")]
      pnms=as.character(1:dim(tmp)[1]) # number parameters by default
      # substitute vitality parameter names
      if(tmp_mod_nm=="vitality.ku"){pnms=c("r","s","k","u")}
      if(tmp_mod_nm=="vitality.4p"){pnms=c("r","s","lambda","beta")}
      if(tmp_mod_nm=="weibull3"){pnms=c("shape","thrsh","scale")}
      dimnames(tmp)=list(pnms,c("est","se"))
      tmp=data.frame(tmp)
      }
      return(tmp)
    },simplify = FALSE)
    tmp_tab=do.call(rbind,par_ls)
    param=unlist(sapply(par_ls,rownames,simplify = F))
    par_tab=data.frame(model=rep(model,sapply(par_ls,nrow)),param,tmp_tab,row.names = NULL)

    out_ls=list("mod_choice"=model,
                "times"=data.frame(time=y,surv_frac=y_sfrac,non_cen=non_cen),
                "fit_vals"=fit_vals,
                "mod_objs"=fit,
                "par_tab"=par_tab,
                "KM_DF"=KM_DF,
                "KM_mod"=KM_mod,
                "censored"=rc)
    out=structure(out_ls,class="fc_list")

    }

  # Table of parameter estimates for the single-model case (Associated with fc_obj)
  else{
    # if only K-M specified
    if(model=="kaplan-meier"){
      par_tab=fit_vals=NULL
      fit[[1]]=KM_mod
      fit_vals=KM_DF
    }
    # Otherwise
    else{
      if(model=="vitality.ku" | model=="vitality.4p" | model =="weibull3"){
        par_tab=fit[[1]][,c("params","std")]
        pnms=as.character(1:dim(par_tab)[1]) # number parameters by default
      if(model=="vitality.ku"){ pnms=c("r","s","k","u")}
      if(model=="vitality.4p"){pnms=c("r","s","lambda","beta")}
      dimnames(par_tab)=list(pnms,c("est","se"))
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

    out=structure(out_ls,class="fc_obj")
  }

return(out)
  # DEV NOTES
  # can use a tapply() to execute function for multiple LotIDs
  # Consider trying alternative optimization methods
}

# function that limits the amount of output displayed when
# a fc_obj or a fc_list is called
#' @title print.fc_list()
#'
#' @description Printed output for class "fc_list"
#'
#' @param x an object of class "fc_list"
#' @param ... ignored
#'
#' @return description of list of models
#'
#' @export
#'
print.fc_list <- function(x,...){
  cat("Failure model list object\n\n")
  cat("Contains the following",length(x[["mod_choice"]]),"models: \n",paste(x[["mod_choice"]],collapse = " ; "))
  if(is.null(x$"GOF_tab")){ cat("\n\n*use this object to compare models using the function: fc_rank()\n")}
  if(!is.null(x$"GOF_tab")){cat("\n\nRanked list\n")
    print(x$"GOF_tab")} # and comment too!
  invisible(x)
}

#' @title print.fc_obj()
#'
#' @description Printed output for class "fc_obj"
#'
#' @param x an object of class "fc_obj"
#' @param ... ignored
#'
#' @export
#'
print.fc_obj <- function(x,...){
  if(x[["mod_choice"]]=="kaplan-meier"){
    cat("Kaplan-Meier estimates for increments between failure times\n")
    print(x[["KM_DF"]][,-c(1)])
  }
  else{
  cat(paste(x[["mod_choice"]],"failure model object \n\n"))
  cat("Parameter estimates:\n")
  print(x[["par_tab"]])
  }
  #cat("\n*This object can be used to adjust survival estimates using the 'ATLAS' package\n")

  invisible(x)
}


#' @title Summary of fc object
#' @details Detailed summary of model objects in fc object
#'
#' @param object object of class fc_obj
#' @param ... ignore
#'
#' @return Summary of fc_obj model fit calls
#' @export
#'
summary.fc_obj <- function(object,...){
  cat("Summary of",paste(object[["mod_choice"]],"failure model object \n\n"))
  print(object$"mod_obj")
  cat("\n*This object can be used to adjust survival estimates using the 'ATLAS' package\n")
  invisible(object)
}

#' @title Sumary of fc list
#' @details Detailed summary of model objects in fc list object
#'
#' @param object object of class fc_list
#' @param ... ignore
#'
#' @return Summary of model fitting calls and GOF rankings (if available)
#' @export
summary.fc_list <- function(object,...){
  cat("Summary failure model list \n\n")
  cat("Contains the following",length(object[["mod_choice"]]),"models: \n",paste(object[["mod_choice"]],collapse = " ; "),"\n\n")
  print(object$"mod_obj")
  if(is.null(object$"GOF_tab")){ cat("\n\n*this object can used to compare model fit using the function: fc_rank()\n")}
  if(!is.null(object$"GOF_tab")){cat("\n\nRanked list\n")
    print(object$"GOF_tab")}
  invisible(object)
}

