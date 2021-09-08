#' @title Simulated Kolmogorov-Smirnov Test
#'
#' @param times numeric vector of failure times
#' @param reps replicates for bootstrap (default to 50k)
#' @param dist distribution
#' @param label optional argument for labeling plots
#'
#' @return P-value and plot of sample distribution of D statistic.
#' 
#' @export
#' 

ks_boot <- function(
  times, # vectors of values
  reps=50000,
  dist="gompertz",
  label=""
){
  n=length(times)
  
  # TWO-PARAMETER MODELS
  if(dist=="gompertz"){
    fit=flexsurv::flexsurvreg(survival::Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pgompertz",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rgompertz(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pgompertz",shape=est_pars[1],rate=est_pars[2])$statistic})}
  
  if(dist=="llogis"){
    fit=flexsurv::flexsurvreg(survival::Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pllogis",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rllogis(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pllogis",est_pars[1],est_pars[2])$statistic})}
  
  if(dist=="lnorm"){
    fit=flexsurv::flexsurvreg(survival::Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"plnorm",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rlnorm(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"plnorm",est_pars[1],est_pars[2])$statistic})}
  
  if(dist=="gamma"){
    fit=flexsurv::flexsurvreg(survival::Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pgamma",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rgamma(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pgamma",est_pars[1],est_pars[2])$statistic})}
  
  if(dist=="weibull(2)"){
    fit=flexsurv::flexsurvreg(survival::Surv(times)~1,dist = "weibull")
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pweibull",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rweibull(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pweibull",est_pars[1],est_pars[2])$statistic})}
  
  # THREE-PARAMETER MODELS
  
  if(dist=="gengamma"){
    fit=flexsurv::flexsurvreg(survival::Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pgengamma",est_pars[1],est_pars[2],est_pars[3])$statistic
    MAT=replicate(reps,rgengamma(n,est_pars[1],est_pars[2],est_pars[3]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pgengamma",est_pars[1],est_pars[2],est_pars[3])$statistic})}
  
  if(dist=="weibull(3)"){
    weib3_res=weibull3_NOSE(times,plots=F)
    est_pars=weib3_res[[1]][,2]
    lines(x=sort(times),y=sort(weib3_res[[3]]$S,decreasing = T),col=2)
    D0=ks.test(times,"pweibull3",est_pars[1],est_pars[2],est_pars[3])$statistic
    MAT=replicate(reps,rweibull3(n,est_pars[1],est_pars[2],est_pars[3]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pweibull3",est_pars[1],est_pars[2],est_pars[3])$statistic})
  }
  
  # FOUR-PARAMETER MODELS
  if(dist=="Vitality09"){
    s_y=sort(times) #sorting taglife values
    y_sfrac=sapply(s_y,function(x){1-length(which(s_y<=x))/length(s_y)})
    est_pars=vitality::vitality.ku(time = sort(s_y),sdata = y_sfrac,se=F,pplot =F,lplot = T, silent = T)
    lines(sort(times),vitality::SurvFn.ku(sort(times),est_pars[1],est_pars[2],est_pars[3],est_pars[4]),col=2)
    D0=ks.test(times,"pvit09",est_pars[1],est_pars[2],est_pars[3],est_pars[4])$statistic
    MAT=replicate(reps,rvitality(
      parms=est_pars, # four vitality parameters
      times_dat=times,  # survival times used for determining # samples to generate and range of slices
      t_seq_fineness=0.001, # time increments to with which to slice up the survival curve
      quant_seq=seq(0,1,0.005),
      model="Vitality09"))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pvit09",est_pars[1],est_pars[2],est_pars[3],est_pars[4])$statistic})
    # print(str(Dsim))
  }
  
  if(dist=="Vitality13"){
    s_y=sort(times) #sorting taglife values
    y_sfrac=sapply(s_y,function(x){1-length(which(s_y<=x))/length(s_y)})
    # est_pars=vitality.4p(time = sort(s_y),sdata = y_sfrac,se=F,pplot =F, silent = T)
    
    est_pars=vitality::vitality.4p(time = s_y,sdata =  y_sfrac,se=F,init.params=c(0.012, 0.01, 0.1, 0.1),
                                   lower = c(0, 0, 0, 0), upper = c(100,50,1,50),rc.data = F,
                                   datatype = "CUM",ttol = 1e-06,silent = T,pplot = F,Iplot = F,Mplot = F)
    
    lines(sort(times),SurvFn.4p(sort(times),est_pars[1],est_pars[2],est_pars[3],est_pars[4]),col=2)
    # print(ks.test(times,"pvit13",est_pars[1],est_pars[2],est_pars[3],est_pars[4]))
    D0=ks.test(times,"pvit13",exact=T,est_pars[1],est_pars[2],est_pars[3],est_pars[4])$statistic
    MAT=replicate(reps,rvitality(
      parms=est_pars, # four vitality parameters
      times_dat=times,  # survival times used for determining # samples to generate and range of slices
      t_seq_fineness=0.001, # time increments to with which to slice up the survival curve
      quant_seq=seq(0,1,0.005),
      model="Vitality13"))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pvit13",exact=T,est_pars[1],est_pars[2],est_pars[3],est_pars[4])$statistic})
    # print(str(MAT))
  }
  
  # P-VAL CALCULATION
  # pval=(table(Dsim>D0)/reps)[2]
  pval=length(which(Dsim>D0))/reps
  
  hist(Dsim,border = NA,col=8,breaks=seq(0,1,0.01),main="",
       col.main=ifelse(pval<0.05,"red","black"),probability=T)
  mtext(side=3,text = paste("P-value =",signif(pval,4)),line = 4,adj = 0.8,col=2)
  abline(v=D0,col=2,lwd=2)
  
  out=list(pval,Dsim,D0)
  names(out)=c("pval","Dsim","D0")
  
  cat("Results of a one-sample Kolmogorov-Smirnov test based on a simulation\n\n")
  cat("model = ",dist,"\n\n")
  cat("iterations = ",reps,"\n\n")
  cat("observed test statistic\n D[obs] =",out$D0,"\n\n")
  cat("p-value = ",pval)
}

#' @title Generating  samples from 2009 and 2013 Vitality models
#'
#' @param parms vector of parameters, Vitality 2009 (r,s,k,u), Vitality 2013 (r,s,lambda,beta)
#' @param times_dat  survival times used for determining # samples to generate and range of slices
#' @param t_seq_fineness time increments to with which to slice up the survival curve
#' @param quant_seq bins in which to place simulated times
#' @param model either "Vitality09" ot "Vitality13"
#'
#' @return
#' @export
rvitality=function(
  parms, # four vitality parameters
  times_dat,  # survival times used for determining # samples to generate and range of slices
  t_seq_fineness=0.005, # time increments to with which to slice up the survival curve
  quant_seq=seq(0,1,0.005), # bins in which to place simulated times
  model="Vitality09"
){
  out=list()
  stopifnot(any(model %in% c("Vitality09","Vitality13")))
  # fineness of slices accross x axis (time)
  # vit_pred_seq=seq(ifelse(min(times_dat)*.8>0,min(times_dat)*.8,0),max(times_dat)*1.2,t_seq_fineness) # span of slices across the x axis
  vit_pred_seq=seq(min(times_dat)*0.8,max(times_dat)*1.2,t_seq_fineness) # span of slices across the x axis
  # evaluates survival curve at slice
  ts=seq(min(times_dat),min(times_dat),.5)
  
  if(model=="Vitality09"){
    pred_survs=vitality::SurvFn.ku(vit_pred_seq,parms[1],parms[2],parms[3],parms[4])}
  
  if(model=="Vitality13"){
    pred_survs=vitality::SurvFn.4p(vit_pred_seq,parms[1],parms[2],parms[3],parms[4])}
  
  # place sequence of times and predicted survival in a dataframe
  vit_sliceDF=data.frame(vit_pred_seq,pred_survs)
  # identify the survival increment that each time interval is in
  vit_sliceDF$bin=cut(vit_sliceDF$pred_survs,breaks=quant_seq)
  # convert bin same to an index
  vit_sliceDF$binID=as.numeric(vit_sliceDF$bin)
  vit_sliceDF$binID[is.na(vit_sliceDF$binID)]=length(quant_seq)+1
  nsamp=length(times_dat)
  bin_samp=sample(unique(vit_sliceDF$binID),nsamp,replace = T)
  bin_samp
  vit_sliceDF
  # sample(subset(vit_sliceDF,binID==bin_samp[1])$vit_pred_seq,1,replace = T)
  # ind_t_samp=sapply(bin_samp,function(x){sample(subset(vit_sliceDF,binID==x)$vit_pred_seq,1,replace = T)})
  # ind_t_samp
}


#'@title Cumulative distribution function of Vitality 2009 model
#'
#' @param x time
#' @param par1 r
#' @param par2 s
#' @param par3 k
#' @param par4 u
#'
#' @return
pvit09=function(x,par1,par2,par3,par4){1-vitality::SurvFn.ku(x,par1,par2,par3,par4)}

#'@title Cumulative distribution function of Vitality 2013 model
#'
#' @param x time
#' @param par1 r
#' @param par2 s
#' @param par3 lambda
#' @param par4 beta
#'
#' @return
pvit13=function(x,par1,par2,par3,par4){1-vitality::SurvFn.4p(x,par1,par2,par3,par4)}

pllogis=flexsurv::pllogis
rllogis=flexsurv::rllogis
