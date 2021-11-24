#' @title Nonparametric bootstrap of failure time model object
#'
#' @param mod_obj failure time model object of class "fc_obj"
#' @param nrep number of replicate samples
#'
#' @return matrix of bootstrap parameter estimates with, dimensions (rows= nrep, columns= estimated pars in mod_obj)
#'
fc_boot=function(mod_obj,nrep){
  stopifnot(class(mod_obj)=="fc_obj")
  obs_dat=mod_obj[["times"]]
  par_mat=NULL
  for(i in 1:nrep){
    ind=sample(1:nrow(obs_dat),replace = T)
    t_raw=obs_dat[ind,1]
    
    t_ord=order(obs_dat[ind,1]) # sorting times
    t_new=fc_surv(t_raw[t_ord])
    boot_dat=data.frame(ind=ind[t_ord],
                        time=t_raw[t_ord],
                        surv_frac=t_new,
                        non_cen=(obs_dat$non_cen[ind])[t_ord])
    # ensure that error-catching works here
    fit=fc_fit(time = boot_dat$time,model = mod_obj$mod_choice)
    # plot(fit)
    par_mat=rbind(par_mat,fit$par_tab[,1])}
  return(par_mat)}
