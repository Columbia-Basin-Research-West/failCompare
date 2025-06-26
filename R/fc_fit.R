#' @title Fitting one or a set of failure time models
#'
#' @description Routines for fitting a common failure time model or models 
#'
#' @param time numeric vector of failure times
#' @param model character  string specififying the model(s) to be fit
#' @param rc.value rc.value right-censoring cutoff value (i.e.,only observations with times > rc.value are censored due to termination of the experiment or study)
#' @param censorID binary or logical variable the same length as \code{time} indicating censored observations, with zeros or FALSE indicating a cenosored observation
#' @param SEs logical for whether standard errors should be estimated
#' @param ... additional arguments passed to optimizer
#'
#' @return Returns failure model object of class \code{"fc_obj"} if one model specified OR 
#' a failure model list object of class \code{"fc_list"} if multiple models are specified.
#'
#' @details
#' This is a model fitting routine used to fit one or a set of failure time models:
#' 
#'  \itemize{
#'     \item  \code{"weibull"}  = 2-parameter Weibull
#'     \item  \code{"weibull3"} = 3-parameter Weibull
#'     \item  \code{"gompertz"}  = Gompertz Model
#'     \item  \code{"gamma"}  = Gamma distribution (2-parameter)
#'     \item  \code{"lognormal"}  = Log-Normal distribution
#'     \item  \code{"llogis"} = Log-Logistic distribution
#'     \item  \code{"gengamma"}  = Generalized Gamma Distribution (3-parameter; Prentice 1974 parameterization)
#'     \item  \code{"vitality.ku"}  = 4-parameter vitality model from Li and Anderson (2009)
#'     \item  \code{"vitality.4p"}  = 4-parameter vitality model from Li and Anderson (2013)
#'     \item  \code{"kaplan-meier"} = Kaplan-Meier nonparametric estimate (NOTE: this model cannot be specified in a list with any other model
#'  }
#'
#' Details on the parameterization of these distributions can be found in the appendix of the 
#' \href{http://www.cbr.washington.edu/analysis/apps/failcompare}{failCompare user manual} .
#' If a single model is specified, a \code{"fc_obj"} is created, which can be
#' used to adjust a CJS model in the "cbrATLAS" package.
#'
#' If multiple models are specified, a \code{"fc_list"} is created containing
#' output from all models that could be fit with default optimizer settings. 
#' A warning will appear if any of the models could not be fit, in which case 
#' the user should either remove the model from consideration or specifiy initial parameter values.
#' 
#' Objects of class \code{fc_list} may serve as an input in the
#' fc_rank() function, which ranks the performance of the model using
#' the \href{http://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-020-00213-z}{Skalski and Whitlock (2020)} GOF measure.
#'
#' Printing a \code{fc_obJ} will display 
#' parameter estimates and accompanying standard errors, if available.
#'
#'   Each \code{fc_obJ} is a list of the following extractable objects:
#'  \itemize{
#'     \item"mod_choice"  = model name
#'     \item     "times"  = dataframe of failure time, survival fraction, and censoring binary var
#'     \item  "fit_vals"  = failure times and predicted survival under the model, 95% LCL an UCL if available
#'     \item  "mod_objs"  = actual model object created by "flexsurvdist" or "vitality package"-- much more to extract from "flexsurvdist
#'     \item   "par_tab"  = table of parameter estimates and SE in failCompare recognized order
#'     \item     "KM_DF"  = table of product limit (Kaplan-Meier)  estimates for plotting (Kaplan and Meier 1954)
#'     \item    "KM_mod"  = survival package K-M model estimates
#'     \item  'censored'  = binary/logical variable the length of the data describing individual observations that are censored
#'  }
#'
#' @references 
#'
#' Kaplan, E.L., and Meier, P. 1958. Nonparametric estimation from incomplete observations. Journal of the American Statistical Association 53(282):457-481.
#'
#' Li, T., and Anderson, J.J. 2009. The vitality model: a way to understand population survival and demographic heterogeneity. Theoretical Population Biology 76(2):118-131.
#'
#' Li, T., and Anderson, J.J. 2013. Shaping human mortality patterns through intrinsic and extrinsic vitality processes. Demographic Research 28:341-372.
#'
#' Prentice, R. L. 1974. A Log Gamma Model and Its Maximum Likelihood Estimation. Biometrika: 61(3):539-544. 
#'
#' Skalski, J. R., and S. L. Whitlock. 2020. Vitality models found useful in modeling tag-failure times in acoustic-tag survival studies. Animal Biotelemetry 8(1):1-10.DOI:10.1186/s40317-020-00213-z.
#'
#'
#' @importFrom survival Surv
#' @importFrom flexsurv flexsurvreg
#' @importFrom vitality vitality.ku
#' @importFrom vitality vitality.4p
#' @export fc_fit 
#'
#'
fc_fit=function(time,model,SEs=TRUE,censorID=NULL,rc.value=NULL,...){
  rc=FALSE #temp def
  if(!is.vector(time)|!is.numeric(time)|any(time<=0)){stop("Expects postive numeric vector the 'time' argument")}
  if(model[1]=="all"){model=c("weibull",'weibull3', "gompertz", "gamma", "lognormal", "llogis", "gengamma","vitality.ku","vitality.4p")
  message("attempting to fit all available parametric survival models")}
  
  if(!is.logical(SEs)){stop:"Expects a logical for the 'SEs' arguement"}
  Hess=SEs

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
  fc_mod_ls=fc_mod_ls

  # HANDLING CENSORING
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
  
  # kaplan-meier estimates
  KM_mod=survival::survfit(survival::Surv(time=y,event=non_cen)~1)
  KM_sls=summary(KM_mod)
  KM_DF=rbind(data.frame(model="kaplan-meier",time=0,est=1,lcl=1,ucl=1),
              data.frame(model="kaplan-meier",time=KM_sls$time,est=KM_sls$surv,lcl=KM_sls$lower,ucl=KM_sls$upper))
  KM_DF[nrow(KM_DF),c("lcl","ucl")]=c(0,0) # replacing last two NAs wih 0s

  if("kaplan-meier" %in% model & length(model)>1){
    model=model[model!="kaplan-meier"]
    message("'kaplan-meier' cannot be included in a fc_list with parametric models, and will be omitted from the string")
  }
  
  ### WEIBULL3 censoring not supported message
  # if(rc & "weibull3" %in% model ){
  if(suppressWarnings(!all(censorID))){
    model=model[model!="weibull3"]
    message("The right-censored 3-parameter Weibull model ('weibull3') is not available ")
    if(length(model)==0){stop("Cannot fit right-censored 3-parameter Weibull model")}
  }
  
  
  # Combined table of parameter estimates for all fitted models (Associated with fc_list)
  if(length(model)==1){
    # if only K-M specified
    if(model=="kaplan-meier"){
      par_tab=fit_vals=NULL
      fit=KM_mod
      fit_vals=KM_DF
      out=list("mod_choice"=model,
                  "times"=data.frame(time=y,surv_frac=y_sfrac,non_cen=non_cen),
                  "fit_vals"=fit_vals,
                  "mod_objs"=fit,
                  "par_tab"=par_tab,
                  "KM_DF"=KM_DF,
                  "KM_mod"=KM_mod,
                  "censored"=rc)
      rownames(out[["par_tab"]])=NULL
      out=structure(out,class="fc_obj")
      return(out) # stops execution here if K-M and only K-M specified
    }
    mc=names(match.call(expand.dots = T))

    out=tryCatch(fc_fit_single(y,y_sfrac,model,Hess,non_cen,KM_DF,KM_mod,...),
                 error = function(x){
                   stop(paste(c(model," model could not be fit\n"),collapse = ""))
                 })
    return(out) # stops execution here if only one parametric modle specified
  }
  else{
    fit=list()
    for (i in 1:length(model)){
      fit[[i]] <- tryCatch(fc_fit_single(y=y,y_sfrac = y_sfrac,non_cen = non_cen,
                    Hess = Hess,KM_DF = KM_DF,KM_mod = KM_mod,
                    model = model[i],...),
                    error = function(e){
                      message(paste(c(model[i]," model could not be fit\n"),collapse = ""))
                    },finally = NULL)
    }
  }
  nonnull_mod_ind <- which(!sapply(fit,is.null))

  # combining models that could be fit (i.e., not returnning an NA)
  if(length(nonnull_mod_ind)>1) {
    out_ls <- fc_combine(fit[!sapply(fit,is.null)]) 
    return(out_ls)}
  
  if(length(nonnull_mod_ind)==1) {
    out_mod <- fit[[nonnull_mod_ind]] 
    return(out_mod)}
  
  stop("no model(s) could be fitted")

}


#' @title Generic function that limits the amount of output displayed when
#' an fc_list is called
#'
#' @description Printed output for class "fc_list"
#'
#' @param x an object of class "fc_list"
#' @param ... additional arguments to print()
#'
#' @return description of list of models
#'
#' @export
#'
print.fc_list <- function(x,...){
  if(is.null(x$"GOF_tab")){
    cat("Failure model list object\n\n")
    cat("Contains the following",length(x[["mod_choice"]]),"models: \n",paste(x[["mod_choice"]],collapse = " ; "))
    cat("\n\n*use this object to compare models using the function: fc_rank()\n")
  }
  if(!is.null(x$"GOF_tab")){
    cat("Failure model list object (Ranked by GOF)\n")
    cat("\nRanked list\n\n")
    print(x$"GOF_tab")} # and comment too!
  invisible(x)
}

#' @title Generic function that limits the amount of output displayed when
#' an fc_obj is called
#'
#' @param x an object of class "fc_obj"
#' @param ... additional arguments to print()
#'
#' @export
#'
print.fc_obj <- function(x,...){
  if(x[["mod_choice"]]=="kaplan-meier"){
    cat("Kaplan-Meier estimates for increments between failure times\n")
    print(x[["KM_DF"]][,-c(1)])
  }
  else{
  cat(paste("Failure model object \n\nType:",x[["mod_choice"]],"\n\n"))
  cat("Parameter estimates:\n")
  print(x[["par_tab"]])
  }
  invisible(x)
}


#' @title Generic function for summarizing an object of class "fc_obj"
#'
#' @param object object of class fc_obj
#' @param ... additional arguments to summary
#'
#' @return Summary of fc_obj model of calls to model fitting functions.
#' @export
#'
summary.fc_obj <- function(object,...){
  cat("Summary of",paste(object[["mod_choice"]],"failure model object \n\n"))
  print(object$"mod_obj")
  cat("\n*This object can be used to adjust survival estimates using the 'ATLAS' package\n")
  invisible(object)
}

#' @title Generic function for summarizing an object of class "fc_list"
#'
#' @param object object of class fc_list
#' @param ... additional arguments to summary
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

