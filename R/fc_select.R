#' @title Selecting a failure time model from a list
#' 
#' @description select failure time model from predefined list of candidate models produced by the function \code{fc_fit()}. Kaplan-Meier nonparametric model is selectable from any list.
#'
#' @param mod_ls failure model list object (i.e., class \code{fc_list})
#' @param model model selected from list of those available. Options include:
#'  \itemize{
#'     \item  \code{"weibull"}  = 2-parameter Weibull
#'     \item  \code{"weibull3"} = 3-parameter Weibull
#'     \item  \code{"gompertz"}  = Gompertz Model
#'     \item  \code{"gamma"}  = Gamma distribution (2-parameter)
#'     \item  \code{"lognormal"}  = Log-Normal distribution
#'     \item  \code{"llogis"} = Log-Logistic distribution
#'     \item  \code{"gengamma"}  = Generalized Gamma Distribution 
#'     \item  \code{"vitality.ku"}  = 4-parameter vitality model 
#'     \item  \code{"vitality.4p"}  = 4-parameter vitality model 
#'     \item  \code{"kaplan-meier"} = Kaplan-Meier nonparametric estimate (always selectable)
#'  }
#' 
#' @return Returns a failure time model object of class \code{fc_obj} .
#' 
#' @seealso
#' \code{\link{fc_fit}}
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
