
#' @title Ranking failure time models based on goodness of fit (GOF)
#'
#' @param x an object of class "fc_list"
#'
#' @details Goodness of fit measure from Skalski and Whitlock (2020)
#'
#' @return Table of models ranked according to GOF measure
#' @export
#'
fc_rank <- function(x){
  stopifnot(class(x)=="fc_list")

  # Estimating Kaplan Meier values for the data
  KM_mod=survival::survfit(survival::Surv(x[["times"]]$time)~1)

  surv_predDF=data.frame(matrix(NA,ncol=6,nrow=length(KM_mod$time)+1)) # adding one for the first part of step function
  names(surv_predDF)=c("model","time","est","lcl","ucl","npars")
  surv_predDF$time=c(0,KM_mod$time);
  surv_predDF$model="Kaplan-Meier"
  surv_predDF$est=c(1,KM_mod$surv);surv_predDF$lcl=c(1,KM_mod$lower)
  surv_predDF$ucl=c(1,KM_mod$upper);surv_predDF$npars=NA

  fit_DF=x[["fit_vals"]]
  # print(match(fit_DF$time,surv_predDF$time))
  # print(str(surv_predDF$est[match(fit_DF$time,surv_predDF$time)]))
  # # fitted values and KM (nonparametric method)
  fit_DF$km_est=surv_predDF$est[match(fit_DF$time,surv_predDF$time)]
  # # # resid
  fit_DF$resid=fit_DF$est-fit_DF$km_est
  # # # squared difference
  fit_DF$sqr_err=(fit_DF$est-fit_DF$km_est)^2

  # summary table
  sumDF=data.frame(
    SSE_KM=tapply(fit_DF$sqr_err[-1],fit_DF$model[-1],sum),
    n=length(x[["times"]]$time),
    npars=tapply(x[["par_tab"]]$param,x[["par_tab"]]$model,length))

  sumDF$denom=(sumDF$n - sumDF$npars -1)
  sumDF$GOF=sumDF$SSE_KM/sumDF$denom

  sumDF=data.frame(model=rownames(sumDF),sumDF); rownames(sumDF)=NULL; sumDF[order(sumDF$GOF,decreasing = F),]
  sumDF$model=factor(sumDF$model,levels = sumDF$model[order(sumDF$GOF,decreasing = T)])
  sumDF$model=factor(as.character(sumDF$model),levels=unique(as.character(sumDF$model)))
  sumDF=sumDF[order(sumDF$GOF,decreasing = F),]
  sumDF$GOF=round(sumDF$GOF,4)
  rownames(sumDF)=NULL

  # out=c(x,structure(list("GOF_tab"=sumDF),class="failmod_list"))
  # out=c(x,structure(list("GOF_tab"=sumDF),class="failmod_list"))

  if(is.null(x$"GOF_tab")){
    out1=c(x,list("GOF_tab"=sumDF))
    out2=structure(out1,class="failmod_list")}
  cat("Candidate models ranked by goodness of fit measure:\n\n")
  return(out2)
}
