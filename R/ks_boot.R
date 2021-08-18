ks_boot
#' Title
#'
#' @param times numeric vector of failure times
#' @param reps replicates for bootstrap (default to 50k)
#' @param dist distribution
#' @param label optional argument for labeling plots
#'
#' @return
#' @export
#'
#' @examples
function(
  times, # vectors of values
  reps=50000,
  dist="gompertz",
  label=""
){
  n=length(times)
  
  plot(Surv(times),xlim=c(min(times)*0.8,max(times)*1.2),main=paste(dist,label,sep="--"))
  
  # TWO-PARAMETER MODELS
  if(dist=="gompertz"){
    fit=flexsurvreg(Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pgompertz",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rgompertz(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pgompertz",shape=est_pars[1],rate=est_pars[2])$statistic})}
  
  if(dist=="llogis"){
    fit=flexsurvreg(Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pllogis",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rllogis(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pllogis",est_pars[1],est_pars[2])$statistic})}
  
  if(dist=="lnorm"){
    fit=flexsurvreg(Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"plnorm",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rlnorm(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"plnorm",est_pars[1],est_pars[2])$statistic})}
  
  if(dist=="gamma"){
    fit=flexsurvreg(Surv(times)~1,dist = dist)
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pgamma",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rgamma(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pgamma",est_pars[1],est_pars[2])$statistic})}
  
  if(dist=="weibull(2)"){
    fit=flexsurvreg(Surv(times)~1,dist = "weibull")
    est_pars=fit$res[,1]
    lines(x=summary(fit)[[1]][,1],y=summary(fit)[[1]][,2],col=2)
    D0=ks.test(times,"pweibull",est_pars[1],est_pars[2])$statistic
    MAT=replicate(reps,rweibull(n,est_pars[1],est_pars[2]))
    Dsim=apply(MAT,2,function(x){ks.test(x,"pweibull",est_pars[1],est_pars[2])$statistic})}
  
  # THREE-PARAMETER MODELS
  
  if(dist=="gengamma"){
    fit=flexsurvreg(Surv(times)~1,dist = dist)
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
    est_pars=vitality.sa(time = sort(s_y),sdata = y_sfrac,se=F,pplot =F,lplot = T, silent = T)
    lines(sort(times),SurvFn(sort(times),est_pars[1],est_pars[2],est_pars[3],est_pars[4]),col=2)
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
    
    est_pars=vitality.4p(time = s_y,sdata =  y_sfrac,se=F,init.params=c(0.012, 0.01, 0.1, 0.1),
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
  
  hist(Dsim,border = NA,col=8,breaks=seq(0,1,0.01),main=paste("P-value =",signif(pval,4)),
       col.main=ifelse(pval<0.05,"red","black"),probability=T)
  abline(v=D0,col=2,lwd=2)
  
  out=list(pval,Dsim,D0)
  names(out)=c("pval","Dsim","D0")
  out
}