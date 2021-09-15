
#' @title Log-rank test of two data sets
#'
#' @param data dataframe containing all variables
#' @param time failure times
#' @param group grouping variable
#' @param censorID censored observation vector (T/F)
#'
#' @return Results of a log-rank test for comparing two survival distributions
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
