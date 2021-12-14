#' @title Log-rank test of two data sets
#'
#' @description A log-rank test of two data sets using the "survival" package
#'
#' @param data dataframe containing all variables
#' @param time numeric failure times
#' @param group character or factor grouping variable
#' @param censorID logical vector the same length as "time" indicating censored observations 
#'
#' @return Returns the results of a log-rank test for comparing two survival distributions.
#' @export
fc_diff <- function(data,time,group,censorID=NULL){
  if(is.null(censorID)){non_cen=rep(1,length(time))
  data$non_cen=non_cen
  }
  Surv=survival::Surv
  survdiff=survival::survdiff
  f1=stats::as.formula(paste("Surv(",time,",event=non_cen)~",group,collapse=""))
  survdiff(f1,data)
}
