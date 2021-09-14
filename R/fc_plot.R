#' @title Ploting failure time and sample survival function
#'
#' @param time 
#' @param surv 
#' @param group 
#'
#' @return histogram of failure times and scatter plot of sample survival function
#' @export
fc_plot=function(time,surv,group=NULL,...){
  
  tmx=max(round(time))
  tmn=min(round(time))
  inc=(tmx-tmn)*0.15
  xdim=c(max(tmn-inc,0),tmx+inc)
  t_brk=seq(xdim[1],xdim[2],length.out = 20)
  
  if(is.null(group)){
    hist(x=time,
         col=8,
         # xlim=xdim,   # x-axis scale
         breaks=20,
         main="Failure time distribution", # plot title
         xlab="Days") # axis label
    
    # Sample survival function
    plot(x=time,y=surv,
         xlim=xdim,  
         main="Survival function", 
         xlab="Days",ylab="S(t)",...)
  }
  else{
    stopifnot(is.character(group)|is.factor(group))
    n_grps=length(unique(group))
    if(length(unique(group))>3){stop("Only 3 groups may be plotted at once")}
    
    # histogram of all failure times
    hist(x = time,breaks=t_brk,xlab="Day",xlim=c(0,60),main="",col=8)
    
    # overlay histogram of spring failure is red
    for(i in 2:n_grps){
      hist(time[group==unique(group)[i]],breaks=t_brk,add=T,col=i)
    legend("topleft",legend=unique(group),fill=c(8,2:(n_grps)),bty="n")
      
    plot(x=time,y=surv,
           xlim=xdim,  
           main="Survival function", 
           xlab="Days",ylab="S(t)",cex=1.1)
    for(i in 2:n_grps){
      points(time[group==unique(group)[i]],y=surv[group==unique(group)[i]],col=i,cex=1.1)}
    legend("bottomleft",legend=unique(group),col=c(1,2:(n_grps)),pch=1,bty="n")
    
    }
  }
}

