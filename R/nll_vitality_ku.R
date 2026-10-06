#' @title Obtain negative log-likelihood from a vitality model fit
#'
#' @param params vitality parameters (r,s,k,u)
#' @param data_time failure times
#' @param non_cen binary indicator variable for non-censored obs (non-censored=1; censored=0)
#' @param indiv_nll logical; indicates whether individual likelihood contribution (defaults to FALSE)
#'
#' @returns negative log likelihood based on individual time
#' @export
#'
#' @examples
#' 
#' data("sockeye")
#' vit_mod=failCompare::fc_fit(time=sockeye$days,model="vitality.ku")
#' 
#' # extract negative-log likelihood value from model object
#' vit_mod$nll
#' 
#' nll_vitality_ku(params = vit_mod$par_tab[, "params"],
#'                 data_time = vit_mod$times$time,
#'                 non_cen = vit_mod$times$non_cen)
#' 
#' # table showing likelihood contributions
#' nll_vitality_ku(params = vit_mod$par_tab[, "params"],
#'                 data_time = vit_mod$times$time,
#'                 non_cen = vit_mod$times$non_cen,
#'                 indiv_nll=TRUE)
#'                 
#' # with censoring
#'                 
#' vit_mod_w_censor=failCompare::fc_fit(time=sockeye$days,model="vitality.ku",rc.value=17)
#' 
#' vit_mod_w_censor$nll
#' 
#' nll_vitality_ku(params = vit_mod$par_tab[, "params"],
#'                 data_time = vit_mod$times$time,
#'                 non_cen = vit_mod$times$non_cen)
#' 
#' # table showing individual likelihood contributions and censoring status
#' nll_vitality_ku(params = vit_mod$par_tab[, "params"],
#'                 data_time = vit_mod$times$time,
#'                 non_cen = vit_mod$times$non_cen,
#'                 indiv_nll=TRUE)
#' 
#'                 
#' 
#' 
nll_vitality_ku <- function(params, data_time, non_cen,indiv_nll=FALSE) {
  r <- params[1]; s <- params[2]; k <- params[3]; u <- params[4]
  
  # Parameter constraints to keep boundaries stable
  if (r <= 0 || s <= 0 || u <= 0 || k < 0) return(Inf)
  
  # Calculate continuous survival probability S(t)
  # St <- s_vitality_ku(data_time, r, s, u, k)
  St <- fc_pred(times = data_time,pars =params,model="vitality.ku")
  
  # Fast numerical derivative to find PDF f(t)
  dt <- 1e-5
  # St_dt <- s_vitality_ku(data_time + dt, r, s, u, k)
  St_dt <- fc_pred(times = data_time + dt,pars =params,model="vitality.ku")
  ft <- -(St_dt - St) / dt
  
  # Safety floors to avoid log(0) numeric failures
  ft <- pmax(ft, 1e-10)
  St <- pmax(St, 1e-10)
  
  ll <- non_cen * log(ft) + (1 - non_cen) * log(St)
  
  if(indiv_nll){
    indiv_DF <- data.frame(data_time,nll=-ll,non_cen)
    return(indiv_DF)
  }
  
  
  # Negative log-likelihood formulation, summed across ALL observations
  jll <- sum(ll)
  return(-jll)
}
