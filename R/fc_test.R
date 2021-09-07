
#' @title Log-rank test of two data sets
#'
#' @param data 
#' @param time 
#' @param group 
#'
#' @return
#' @export
fc_test <- function(data,time,group,non_cen=NULL){
  if(is.null(non_cen)){non_cen=rep(1,length(time))
  data$non_cen=non_cen
  }
 # print(paste("Surv(",time,",event=non_cen)~",group,collapse=""))
  f1=as.formula(paste("Surv(",time,",event=non_cen)~",group,collapse=""))
  survdiff(f1,data)
}
