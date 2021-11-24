#' @title Combination of multiple failure time model objects into a list of models
#'
#' @details A convenience function for combining model failure time model objecs \code{fc_obj} into a failure model
#' list object \code{fc_list}. Lists that include the "kaplan-meier" model or duplicates are not allowed.
#'
#' @param mod_ls list of fc_mods
#'
#' @return fc_list object
#' 
#' @seealso \code{fc_select} and \code{fc_fit}
#' 
#' @references 
#' Li, T., and Anderson, J.J. 2009. The vitality model: A way to understand population survival and demographic heterogeneity. Theoretical Population Biology 76(2): 118–131.
#'
#' Li, T., and Anderson, J.J. 2013. Shaping human mortality patterns through intrinsic and extrinsic vitality processes. Demographic research 28: 341–372.
#'
#' @export
#' 
fc_combine <- function(mod_ls){
  #checking for too many models in the list
  if(length(mod_ls)>8){
    stop("list of models exceeds the total number of default failCompare models 
    (combining and rank models with different sample sizes/censoring schemes is not advised)")}
  #checking for only fc_obj
  if(!all(sapply(mod_ls,function(x){class(x)=="fc_obj"}))){stop("Expecting a list of model objects")}
  modnms=c('weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma','vitality.ku','vitality.4p')
  #checking that only default parametric models are used
  if(!all(sapply(mod_ls,function(x){x$mod_choice %in% modnms}))){stop(paste("Expecting parametric model objects    of type:",paste(modnms,collapse="; ")))}
  
  fit_vals=do.call(rbind,lapply(mod_ls,function(x){x$fit_vals}))
  
  out_ls=list("mod_choice"=sapply(mod_ls,function(x){x$mod_choice}),
              "times"=mod_ls[[1]]$times,
              "fit_vals"=do.call(rbind,lapply(mod_ls,function(x){x$fit_vals})),
              "mod_objs"=do.call(c,lapply(mod_ls,function(x){x$mod_objs})),
              "par_tab"=do.call(rbind,lapply(mod_ls,function(x){x$par_tab})),
              "KM_DF"=mod_ls[[1]]$KM_DF,
              "KM_mod"=mod_ls[[1]]$KM_mod) # advisable to check that all km_mods are the same
  
  if(!any(names(out_ls[["par_tab"]]) %in% "model")){
    out_ls[["par_tab"]]=data.frame(model=rep(out_ls[["mod_choice"]],sapply(mod_ls,function(x){nrow(x$par_tab)})),
                                   param=unlist(sapply(mod_ls,function(x){rownames(x$par_tab)})),#rownames(out_ls[["par_tab"]]),
                                   out_ls[["par_tab"]])
    rownames(out_ls[["par_tab"]])=NULL}
  
  out=structure(out_ls,class="fc_list") # adds fc_list class designation
  
  return(out)
}
