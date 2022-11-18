#' @title Compute sample survival function
#' 
#' @description Computes a sample survival function.
#' 
#' @details Calculates a sample survival function accounting for right censoring. In the absence of censoring, it uses the basic survival function estimator, or otherwise uses the Kaplan-Meier product limit estimate.
#'
#' @param time failure or censoring time
#' @param censorID binary or logical variable the same length as \code{time} indicating censored observations, with zeros or FALSE indicating a censored observation
#' @param rc.value time after which all values are censored
#'
#' @return a numeric vector of survival fraction estimates sample survival function
#' 
#' @export
#'
fc_surv <- function(time,censorID=NULL,rc.value=NULL){
  if(is.unsorted(time)){message("Warning: times are not sorted")}
  org_time_ord=order(time)
  y=time
  y_sfrac=sapply(y,function(x){1-length(which(y<=x))/length(y)}) # survival fraction calc
  
  # Checking censoring values
  if(is.null(censorID)){non_cen=rep(1,length(y))}
  else{
    if(any(sapply(censorID,function(x){!(x %in% c(0,1) | is.logical(x))}))){stop("1/0 or TRUE/FALSE expected for censorID")}
    if(!is.null(rc.value)){stop("'rc.value' will override 'censorID' argument")}
    non_cen=censorID}
  
  if(!is.null(rc.value)){
    rc=TRUE # change this value for later if statement
    non_cen=ifelse(time<rc.value,1,0) # vector used by "flexsurv"
    y=sort(y)  # sorted data necessary for Vitality package functions
    y_sfrac=sapply(y,function(x){1-length(which(y<=x))/length(y)}) # survival fraction calc
    
    # For vitality model
    y_cen=y[y<rc.value]
    y_cen_sfrac=y_sfrac[y<rc.value]
    n_cen=length(y)
  }
  
  KM_mod=survival::survfit(survival::Surv(time=y,event=non_cen)~1)
  KM_sls=summary(KM_mod)
  S_est=KM_sls$surv[match(y,KM_sls$time)]

  return(S_est)
}
