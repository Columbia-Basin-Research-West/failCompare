#' @title Look-up parameter names for gompertz model
#' 
#' @param mod_nm character vector of acceptable failCompare model
#'
#' @return dataframe of with model name "model" and "param" names in order reported by optimizer
#' 
#' 
get_param_nm <- function(mod_nm){
  pos_mods=names(fc_mod_ls) #accepted models
  if(!is.character(mod_nm)){
    stop("Expects a character vector as an input")}
  if(any(duplicated(mod_nm)) | any(!mod_nm %in% pos_mods)){
    stop(paste("Expects a model name or nonrepeating vector  of following inside of c():\n'",
               paste(pos_mods,collapse= "', '"),"'",sep=""))}
  # fc_mod_ls is a named list automatically loaded with the package 
  tab_ls=sapply(mod_nm,function(x){
    param=as.character(fc_mod_ls[[x]])
    model=rep(x,length(par))
    cbind(model,param)}
    ,simplify = F,USE.NAMES = F)
  out=tab_ls
  if(length(mod_nm)==1){out=data.frame(tab_ls)}
  else{out=data.frame(do.call(rbind,tab_ls))}
  return(out)
}

