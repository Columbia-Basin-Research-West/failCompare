#' @title Combination of multiple failure time model objects into a list of models
#'
#' @description A combination of multiple failure time model objects into a list of models.
#'
#' @details A convenience function for combining model failure time model objecs \code{fc_obj} into a failure model
#' list object \code{fc_list}. Lists that include the "Kaplan-Meier" model or duplicates are not allowed.
#'
#' @param mod_ls list of fc_mods
#'
#' @return fc_list object
#' 
#' @seealso \code{fc_select} and \code{fc_fit}
#' 
#' @references 
#' Li, T., and Anderson, J.J. 2009. The vitality model: a way to understand population survival and demographic heterogeneity. Theoretical Population Biology 76(2):118–131.
#'
#' Li, T., and Anderson, J.J. 2013. Shaping human mortality patterns through intrinsic and extrinsic vitality processes. Demographic Research 28:341–372.
#'
#' @examples 
#' 
#' ### Load example dataframe
#' data(sockeye)
#' taglife=sockeye[,"days"] #define vector of times
#' 
#' ### Fit a 2-parameter Weibull model
#' weib_mod=fc_fit(time=taglife,model="weibull")
#' 
#' ### Fit a 4-parameter Vitality 2013 model
#' vit_mod=fc_fit(time=taglife,model="vitality.4p")
#' 
#' # Combine two "fc_obj" objects into a model list of class "fc_list" 
#' fc_combine(mod_ls = list(weib_mod,vit_mod))
#' 
#' @export
#' 
fc_combine <- function(mod_ls){
  
  ls_nms=as.character(sapply(mod_ls,function(x){x$mod_choice}))
  
  #checking for too many models in the list
  if(length(mod_ls)<2){stop("expects more than 1 model object in a list object")}
  if(length(mod_ls)>9){
    stop("list of models exceeds the total number of default failCompare models 
    (combining and rank models with different sample sizes/censoring schemes is not advised)")}
  #checking for only fc_obj
  if(!all(sapply(mod_ls,function(x){class(x)=="fc_obj"})) | any(duplicated(ls_nms)))
    {stop("expecting a list of unique model objects")}
  modnms=names(fc_mod_ls)
  #checking that only default parametric models are used
  if(!all(sapply(mod_ls,function(x){x$mod_choice %in% modnms}))){stop(paste("expecting parametric model objects of type:",paste(modnms,collapse="; ")))}
  
  out_ls=list(
              "mod_choice"=ls_nms,
              "times"=mod_ls[[1]]$times,
              "fit_vals"=do.call(rbind,lapply(mod_ls,function(x){x$fit_vals})),
              "mod_objs"=do.call(c,lapply(mod_ls,function(x){list(x$mod_objs)})),
              "par_tab"=do.call(rbind,lapply(mod_ls,function(x){x$par_tab})),
              "KM_DF"=mod_ls[[1]]$KM_DF,
              "KM_mod"=mod_ls[[1]]$KM_mod) # advisable to check that all km_mods are the same
  out_ls[["par_tab"]]=data.frame(
    param=unlist(sapply(ls_nms,function(x){fc_mod_ls[[x]]})),
    model=rep(x = ls_nms,times=sapply(ls_nms,function(x){length(fc_mod_ls[[x]])})),
    out_ls[["par_tab"]])

  # making rownames comparable
  rownames(out_ls[["par_tab"]])=NULL
  rownames(out_ls[["fit_vals"]])=1:nrow(out_ls[["fit_vals"]]) # renumbering
  
  out=structure(out_ls,class="fc_list") # adds fc_list class designation
  
  return(out)
}

