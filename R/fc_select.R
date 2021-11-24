#' @title Selecting a failure time model from a list
#' 
#' @description Select a failure time model from predefined list of candidate models produced by the function \code{fc_fit()}
#'
#' @param mod_ls failure model list object (i.e., class \code{fc_list})
#' @param model model selected from list of those available. Possible options include:
#'  \itemize{
#'     \item   "weibull"  = 2-parameter Weibull
#'     \item   "weibull3" = 3-parameter Weibull
#'     \item  "gompertz"  = Gompertz Model
#'     \item     "gamma"  = Gamma distribution (2-parameter)
#'     \item "lognormal"  = Log-Normal distribution
#'     \item     "llogis" = Log-Logistic distribution
#'     \item  "gengamma"  = Generalized Gamma Distribution (3-parameter; Prentice 1974 parameterization)
#'     \item  'vitality.ku'  = 4-parameter vitality model from Li and Anderson (2009)
#'     \item  'vitality.4p'  = 4-parameter vitality model from Li and Anderson (2013)
#'  }
#' 
#' @return Failure time model object of class \code{fc_obj}
#' 
#' @export
#'
fc_select <- function(mod_ls,model){
  # validation
  stopifnot(is.character(model) & length(model)==1)
  if(!model %in% c(mod_ls$mod_choice,"kaplan-meier")){
    stop(paste("model is not in the list","\n select one of the following: ",
               paste(c(mod_ls$mod_choice,"kaplan-meier"),collapse = ","),sep=""))}
  # subsetting list object
  
  par_tab=mod_ls$parm_tab[mod_ls$parm_tab$model==model,]
  rownames(par_tab)=par_tab[,2]
  
  if(model=="kaplan-meier"){
    out_ls = list(mod_choice = model,
                  times = mod_ls$times,
                  fit_vals = mod_ls$KM_DF,
                  mod_objs =mod_ls$KM_mod,
                  par_tab = NULL,
                  KM_DF=mod_ls$KM_DF,
                  KM_mod=mod_ls$KM_mod)
    out = structure(out_ls, class = "fc_obj")
  }
  else{
    out_ls = list(mod_choice = model,
                  times = mod_ls$times,
                  fit_vals = mod_ls$fit_vals[mod_ls$fit_vals$model==model,],
                  mod_objs = mod_ls$mod_objs[[model]],
                  par_tab = mod_ls$par_tab[mod_ls$par_tab$model==model,-c(1:2)],
                  KM_DF=mod_ls$KM_DF,
                  KM_mod=mod_ls$KM_mod)
    # reorder rownames from 1
    rownames(out_ls$fit_vals)=1:nrow(out_ls$fit_vals)
    rownames(out_ls$par_tab)=1:nrow(out_ls$par_tab)
    out = structure(out_ls, class = "fc_obj")}
  
  return(out)
}
