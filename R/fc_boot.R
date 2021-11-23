library(failCompare)

# loading ATLAS example data set
time=read.csv("http://www.cbr.washington.edu/sites/default/files/Taglife%20ATLAS%20Practice%20File.csv")
taglife=time$tag_life_days

######################################## #
#  multiple model fitting and ranking
######################################## #

# fitting and ranking example (excluding vitality models, which need initial values)
mod_ls=failCompare::fc_fit(time=taglife,model=c('weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma'))
mod_ls_ranked=fc_rank(mod_ls)
# "mod_ls" and "mod_ls_ranked" are "fc_list" objects used for batch-fitting and rank models in failcompare
class(mod_ls)

########################### #

# fitting a 2-parameter weibull by itself
TL_weib2=failCompare::fc_fit(time=taglife,model="weibull")

# selecting the gengamma  from the model list
TL_gengam=fc_select(mod_ls_ranked,model = "gengamma")

# fc_obj objects CAN be used to correct taglife
class(TL_weib2)=="fc_obj"
class(TL_gengam)=="fc_obj"

# weibull3 failure model object 
TL_weib2
#
# Printing the objects displays this
# Parameter estimates:
#   est        se
# 1 12.0898 1.3689615
# 3 16.3512 0.1986752

# Many different objects hidden inside
names(TL_weib2)

# "mod_choice"  = official model name
#      "times"  = DF of failure time, survival fraction, and censoring binary var.
#   "fit_vals"  = failure times and predicted survival under the model, 95% LCL an UCL if available
#   "mod_objs"  = actual model object created by "flexsurvdist" or "vitality package"-- much more to extract from "flexsurvdist
#    "par_tab"  = table of parameter estimates and SE in failCompare recognized order
#      "KM_DF"  = table of K-M estimates for plotting
#     "KM_mod"  = survival package K-M model estimates
#   'censored'  = T/F of censored dataset or not


plot(TL_gengam)

# survival times
time_preds=c(18:27)

# predicted survival probability over time
tag_surv=failCompare::fc_pred(mod_obj = TL_weib2,times = time_preds)

tag_surv





TL_gengam$times$non_cen
length(TL_gengam$times$non_cen)
TL_gengam$times$non_cen=paste(rep(letters,4),rep(1:4,each=26))[1:100]

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


par(mfrow=c(3,3))
boot_out=fc_boot(mod_obj=TL_gengam,nrep = 100)

t(apply(boot_out,1,function(x){fc_pred(times = time_preds,pars=x,model = "gengamma")}))

t(apply(apply(boot_out,1,function(x){fc_pred(times = time_preds,pars=x,model = "gengamma")}),
        1,function(x){c(mean(x),quantile(x,probs = c(0.025,0.975)))}))

# plot(fc_fit(time = time_preds,model = "gengamma"))
