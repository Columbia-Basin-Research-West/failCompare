#' @title Ranking failure time models 
#'
#' @description This provides a ranking of failure time models.
#'
#' @param x an object of class \code{"fc_list"}
#'
#' @details This uses the Skalski and Whitlock (2020) goodness-of-fit measure to rank parametric failure time models. 
#' The statistic is based on the squared difference between a given parametric model and the product-limit estimate
#' of the survival estimate described by Kaplan and Meier (1954).
#'
#' @return Creates a table of models ranked in ascending order according to a goodness-of-fit measure.
#' 
#' @export
#'
#' @references 
#' 
#' Kaplan, E.L., and Meier, P. 1958. Nonparametric estimation from incomplete observations. Journal of the American Statistical Association 53(282):457-481.
#' 
#' Skalski, J. R., and S. L. Whitlock. 2020. Vitality models found useful in modeling tag-failure times in acoustic-tag survival studies. Animal Biotelemetry 8(1):1-10. doi:10.1186/s40317-020-00213-z#'
#'
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
  # # fitted values and KM (nonparametric method)
  fit_DF$km_est=surv_predDF$est[match(fit_DF$time,surv_predDF$time)]
  
  # Checking for censored columns
  fit_DF$non_cen=FALSE
  fit_DF$non_cen[(fit_DF$time!=0)]=rep(x[["times"]]$non_cen,length(x$mod_choice))

  # SUBSET to non-zeros and uncensored obs.
  fit_DF=fit_DF[fit_DF$non_cen,]
  
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

  if(is.null(x$"GOF_tab")){
    out1=c(x,list("GOF_tab"=sumDF))
    out2=structure(out1,class="fc_list")}
  # cat("Candidate models ranked by goodness of fit measure:\n\n")
  # print(sumDF)
  return(out2)
}

