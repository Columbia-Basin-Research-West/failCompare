#' @title Simulated Kolmogorov-Smirnov Test
#' 
#' 
#' @param times numeric vector of failure times
#' @param iters replicates for bootstrap (default to 50k)
#' @param model distribution
#' @param label optional argument for labeling plots
#' @param plot optional argument for creating histogram
#'
#' @details performs a a simulation-based Kolmogorov-Smirnov test.
#'
#' @return p-value and optionallyhistogram of based on a Monte Carlo estimate of the sampling distribution of the D statistic.
#'
#'
#' @seealso \code{\link[stats]{ks.test}}.
#'  \code{\link[stats]{rweibull}}.
#' 
#' 
#' @importFrom stats ks.test
#' @importFrom stats rlnorm
#' @importFrom stats rgamma
#' @importFrom stats rweibull
#'
#' @importFrom flexsurv flexsurvreg
#' @importFrom flexsurv pllogis
#' @importFrom flexsurv rllogis
#' @importFrom flexsurv rgompertz
#' @importFrom flexsurv rgengamma
#'
#'
#' @importFrom stats rweibull
#' @import graphics
#'
#' @export
#'  
fc_test <- function(
  times, # vectors of values
  iters=50000,
  model="gompertz",
  label="",
  plot=FALSE
){
  test=!model %in% c("weibull",'weibull3', "gompertz", "gamma", "lognormal", "llogis", "gengamma","vitality.ku","vitality.4p","kaplan-meier")
  if(any(test)){
    message(paste("Model names not recognized:",paste(model[which(test)],collapse=";")," \n Default model names = {'weibull','weibull3','gompertz','gamma','lognormal','llogis','gengamma','vitality.ku','vitality.4p','kaplan-meier'}",sep=""))
    model=model[!test]
    if(length(model)==0){stop("No valid names in the 'model' argument")}}
  n=length(times)
  
  # TWO-PARAMETER MODELS
  if(model=="gompertz"){
    fit=flexsurvreg(survival::Surv(times)~1,dist = model)
    est_pars=fit$res[,1]
    D0=ks.test_fc(times,"pgompertz",est_pars[1],est_pars[2])$statistic
    MAT=replicate(iters,rgompertz(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"pgompertz",shape=est_pars[1],rate=est_pars[2])$statistic})}
  
  if(model=="llogis"){
    fit=flexsurvreg(survival::Surv(times)~1,dist = model)
    est_pars=fit$res[,1]
    D0=ks.test_fc(times,"pllogis",est_pars[1],est_pars[2])$statistic
    MAT=replicate(iters,rllogis(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"pllogis",est_pars[1],est_pars[2])$statistic})}
  
  if(model=="lnorm"){
    fit=flexsurvreg(survival::Surv(times)~1,dist = model)
    est_pars=fit$res[,1]
    D0=ks.test_fc(times,"plnorm",est_pars[1],est_pars[2])$statistic
    MAT=replicate(iters,rlnorm(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"plnorm",est_pars[1],est_pars[2])$statistic})}
  
  if(model=="gamma"){
    fit=flexsurvreg(survival::Surv(times)~1,dist = model)
    est_pars=fit$res[,1]
    D0=ks.test_fc(times,"pgamma",est_pars[1],est_pars[2])$statistic
    MAT=replicate(iters,rgamma(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"pgamma",est_pars[1],est_pars[2])$statistic})}
  
  if(model=="weibull"){
    fit=flexsurvreg(survival::Surv(times)~1,dist = "weibull")
    est_pars=fit$res[,1]
    D0=ks.test_fc(times,"pweibull",est_pars[1],est_pars[2])$statistic
    MAT=replicate(iters,rweibull(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"pweibull",est_pars[1],est_pars[2])$statistic})}
  
  # THREE-PARAMETER MODELS
  
  if(model=="gengamma"){
    fit=flexsurvreg(survival::Surv(times)~1,dist = model)
    est_pars=fit$res[,1]
    D0=ks.test_fc(times,"pgengamma",est_pars[1],est_pars[2],est_pars[3])$statistic
    MAT=replicate(iters,rgengamma(n,est_pars[1],est_pars[2],est_pars[3]))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"pgengamma",est_pars[1],est_pars[2],est_pars[3])$statistic})}
  
  if(model=="weibull3"){
    # weib3_res=weibull3_NOSE(times,plots=F)
    weib3_res=taglife.fn_weib3(tags.in = times)
    est_pars=weib3_res[[1]][,2]
    D0=ks.test_fc(times,"pweibull3",est_pars[1],est_pars[2],est_pars[3])$statistic
    MAT=replicate(iters,rweibull3(n,est_pars[1],est_pars[2],est_pars[3]))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"pweibull3",est_pars[1],est_pars[2],est_pars[3])$statistic})
  }
  
  # FOUR-PARAMETER MODELS
  if(model=="Vitality09"){
    s_y=sort(times) #sorting taglife values
    y_sfrac=sapply(s_y,function(x){1-length(which(s_y<=x))/length(s_y)})
    est_pars=vitality::vitality.ku(time = sort(s_y),sdata = y_sfrac,se=F,pplot =F,lplot = T, silent = T)
    D0=ks.test_fc(times,"pvit09",est_pars[1],est_pars[2],est_pars[3],est_pars[4])$statistic
    MAT=replicate(iters,rvitality(
      parms=est_pars, # four vitality parameters
      times_dat=times,  # survival times used for determining # samples to generate and range of slices
      t_seq_fineness=0.001, # time increments to with which to slice up the survival curve
      quant_seq=seq(0,1,0.005),
      model="Vitality09"))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"pvit09",est_pars[1],est_pars[2],est_pars[3],est_pars[4])$statistic})
    # print(str(Dsim))
  }
  
  if(model=="Vitality13"){
    s_y=sort(times) #sorting taglife values
    y_sfrac=sapply(s_y,function(x){1-length(which(s_y<=x))/length(s_y)})
    # est_pars=vitality.4p(time = sort(s_y),sdata = y_sfrac,se=F,pplot =F, silent = T)
    
    est_pars=vitality::vitality.4p(time = s_y,sdata =  y_sfrac,se=F,init.params=c(0.012, 0.01, 0.1, 0.1),
                                   lower = c(0, 0, 0, 0), upper = c(100,50,1,50),rc.data = F,
                                   datatype = "CUM",ttol = 1e-06,silent = T,pplot = F,Iplot = F,Mplot = F)
    
    D0=ks.test_fc(times,"pvit13",exact=T,est_pars[1],est_pars[2],est_pars[3],est_pars[4])$statistic
    MAT=replicate(iters,rvitality(
      parms=est_pars, # four vitality parameters
      times_dat=times,  # survival times used for determining # samples to generate and range of slices
      t_seq_fineness=0.001, # time increments to with which to slice up the survival curve
      quant_seq=seq(0,1,0.005),
      model="Vitality13"))
    Dsim=apply(MAT,2,function(x){ks.test_fc(x,"pvit13",exact=T,est_pars[1],est_pars[2],est_pars[3],est_pars[4])$statistic})
    # print(str(MAT))
  }
  
  # P-VAL CALCULATION
  pval=length(which(Dsim>D0))/iters
  
  if(plot){
    hist(Dsim,col=8,breaks=seq(0,1,0.01),main="",xlab="D",
         col.main=ifelse(pval<0.05,"red","black"),probability=T)
    mtext(side=3,text = paste("P-value =",signif(pval,4)),line = -2,adj = 0.8,col=2)
    abline(v=D0,col=2,lwd=2)
  }
  
  out=list(pval,Dsim,D0)
  names(out)=c("pval","Dsim","D0")
  
  cat("Results of a one-sample Kolmogorov-Smirnov test based on a simulation\n\n")
  cat("model = ",model,"\n\n")
  cat("iterations = ",iters,"\n\n")
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
#' @return random values
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
#' @return cumulative probability
pvit09=function(x,par1,par2,par3,par4){1-vitality::SurvFn.ku(x,par1,par2,par3,par4)}

#'@title Cumulative distribution function of Vitality 2013 model
#'
#' @param x time
#' @param par1 r
#' @param par2 s
#' @param par3 lambda
#' @param par4 beta
#'
#' @return cumulative probability
pvit13=function(x,par1,par2,par3,par4){1-vitality::SurvFn.4p(x,par1,par2,par3,par4)}

#' @title random number generation for 3-parameter weibull
#'
#' @param n sample size
#' @param shape beta
#' @param scale lambda
#' @param thres gamma
#'
#' @return vector of random values from the 3-parameter weibull model
rweibull3=function(n, shape, scale = 1, thres = 0){
  thres + rweibull(n, shape, scale)
}

#' @title ks.test with suppressed warnings
#'
#' @param ... inputs to stats::ks.test() function
#'
#' @return expected output from ks.test
#'
ks.test_fc=function(...){
  suppressWarnings(ks.test(...))
}


  