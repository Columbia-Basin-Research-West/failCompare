#' @title Error handling for fitting failCompare models
#' 
#' @details Prevents errors from interupting single- and multi-model runs using fc_fit
#'
#' @param fit_call call to dependent model fitting functions.
#' @param model model argument passed from fc_fit()
#' @param non_cen logical indicating censored variables for use by flexsurv and vitality models
#'
#' @return model fitting output for internal use by fc_fit
#'
fc_tryfit=function(y=y,fit_call,model="weibull3",non_cen=NULL,Hess=NULL){
  # seeing if an error in generated
  frst_ft=tryCatch(eval(fit_call),
                   error = function(e) {
                     disp=paste(c("Error(s) in ",model," model fitting:\n"),collapse = "")
                     msg=conditionMessage(e)
                     message(paste(disp,unique(msg)))
                     err=invisible(structure(msg, class = "try-error"))
                     return(err)})
  # If theres a optim error replacing call and reporting warnings
  if(class(frst_ft)=="try-error"){
    if(model=="weibull3"){fit_call$tag.se=F}
    if(model %in% names(flexsurv::flexsurv.dists)){fit_call$hessian=F}
    if(length(grep(model,pattern = "vitality"))==1){fit_call$se=F}
    message("Switching 'SEs' argument to FALSE to avoid error")
  # Second try at optimizing
    sec_ft=tryCatch(eval(fit_call),
                   error = function(e) {
                     disp=paste(c("Error(s) in ",model,"model fitting:\n"),collapse = "")
                     msg=conditionMessage(e)
                     message(paste(disp,unique(msg),"\nTry changing initial values or optimizer settings"))
                     err=invisible(structure(msg, class = "try-error"))
                     return(err)},
                    warning = function(w) {
                      disp=paste(c("Warnings(s) in ",model," model fitting:\n"),collapse = "")
                      msg=conditionMessage(w)
                      message(paste(disp,unique(msg)))
                      formals(taglife.fn_weib3)[["tag.se"]]
                      finally = eval(fit_call)})
    
    lst_ft=sec_ft}
  else{lst_ft=frst_ft}
  return(lst_ft)
}

