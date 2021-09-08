#' @title Compute sample survival function
#' 
#' @details Calculates sample survival function accounting for right censoring. In the absence of censoring, uses the basic survival function estimator, otherwise uses the Kaplan-Meier product limit estimate.
#'
#' @param time failure or censoring time
#' @param non_cen logical vector the same length as `time`, with TRUE indicating censoring time
#'
#' @return sample survival function
#' 
#' @export
#'
fc_surv <- function(time,non_cen=NULL,rt.value=NULL,rc.value=NULL){
  org_time_ord=order(time)
  y=time
  y_sfrac=sapply(y,function(x){1-length(which(y<=x))/length(y)}) # survival fraction calc
  
  if(is.null(non_cen)){non_cen=rep(1,length(y))}
  
  if(!is.null(rt.value)){
    non_trunc=time<rt.value # vector used by "flexsurv"
    y=y[non_trunc]
    non_cen=non_cen[non_trunc]
    # recalculate the survival fraction after truncating values beyond the threshold
    y_sfrac=sapply(y,function(x){1-length(which(y<=x))/length(y)}) # survival fraction calc
  }
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
  # KM_DF=rbind(data.frame(model="kaplan-meier",time=0,est=1,lcl=1,ucl=1),
  #             data.frame(model="kaplan-meier",time=KM_sls$time,est=KM_sls$surv,lcl=KM_sls$lower,ucl=KM_sls$upper))
  # KM_DF[nrow(KM_DF),c("lcl","ucl")]=c(0,0) # replacing last two NAs wih 0s
  
  S_est=KM_sls$surv[match(y,KM_sls$time)]

  return(S_est)
}
