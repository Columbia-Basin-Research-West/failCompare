#' @title Error handling for fitting failCompare models
#' 
#' @param fit_call call to dependent model fitting functions.
#' @param model model argument passed from fc_fit()
#' @param non_cen logical indicating censored variables for use by flexsurv and vitality models
#' @param y numeric time argument of failure times carried through
#' @param y_sfrac survival fraction
#' @param Hess logical argument to fc_fit() carried through
#' @param inits initial parameter values for optimization of a single model
#' @param ... additional arguments passed to dependent functions
#'
#' @return model fitting output for internal use by fc_fit
#'
#' @details Prevents errors from interrupting single- and multi-model runs using fc_fit
#' 
#'
fc_tryfit=function(y=y,y_sfrac=NULL,fit_call,model="weibull3",non_cen=NULL,Hess=NULL,inits=NULL,...){
  frst_ft=tryCatch(eval(fit_call),
                   error = function(e) {
                     disp=paste(c("Encountered error(s) when trying to fit ",model," model:\n"),collapse = "")
                     msg=conditionMessage(e)
                     message(paste(disp,unique(msg)))
                     err=invisible(structure(msg, class = "try-error"))
                     return(err)},
                   warning = function(w) {
                     disp=paste(c("Encountered warning(s) when trying to fit ",model," model:\n"),collapse = "")
                     msg=conditionMessage(w)
                     message(paste(disp,unique(msg)))
                     suppressWarnings(eval(fit_call))
                     # return()
                   }
                   )
  # If there's an optim error replacing call and reporting warnings
  if(any(class(frst_ft)=="try-error")){
    # if(Hess){message("Switching 'SEs' argument to FALSE to avoid error")}
    if(model=="weibull3"){fit_call$tag.se=F}
    if(model %in% names(flexsurv::flexsurv.dists)){fit_call$hessian=F}
    if(length(grep(model,pattern = "vitality"))==1){fit_call$se=F}

  # Second try at optimizing
    sec_ft=tryCatch(eval(fit_call),
                   error = function(e) {
                     # disp=paste(c("Error(s) in ",model," model fitting:\n"),collapse = "")
                     # msg=conditionMessage(e)
                     # message(paste(disp,unique(msg),"\nTry changing initial values or optimizer settings"))
                     message("consider setting initial values or altering optimizer settings")
                     err=invisible(structure(msg, class = "try-error"))
                     return(err)},
                    warning = function(w) {
                      disp=paste(c("Warning(s) in ",model," model fitting:\n"),collapse = "")
                      msg=conditionMessage(w)
                      message(paste(disp,unique(msg)))
                      }
                   )
    lst_ft=sec_ft}
  else{lst_ft=frst_ft}
  return(lst_ft)
}

