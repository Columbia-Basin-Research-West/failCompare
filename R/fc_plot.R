#' @title Plotting failure time and sample survival function
#'
#' @details plot of type "data" shown by default. For "residual" type plot showing (Kaplan-Meier estimates - parametric model fit), plot a singular model of class = fc_obj.
#'
#' @param time failure time (x axis)
#' @param surv survival function (y axis)
#' @param group grouping variable, limit of 3
#' @param censorID binary or logical variable the same length as \code{time} indicating censored observations, with zeros or FALSE indicating a censored observation
#' @param hist logical show histogram of failure times 
#' @param surv_curv logical show histogram of failure times 
#' @param xlim x axis limits for survival plot
#' @param ylim y axis limits for survival plot, used to override default of c(0,1)
#' @param main title for scatterplot
#' @param xlab x axis limits for survival plot
#' @param ylab y axis limits for survival plot, used to override default of c(0,1)
#' @param ... additional argument passed to plot()
#'
#' @return histogram of failure times and/or scatter plot of sample survival function.
#' @export
fc_plot=function(time,surv,censorID,group=NULL,hist=T,surv_curv=T,main,ylim,xlim,ylab,xlab,...){
  if(all(!surv_curv,!hist)){stop("At least 'surv_curv' or 'hist' must be TRUE")}
  if(missing(censorID)){censorID=rep(FALSE,length(time))}
  else{if(!is.null(group)){warning("group symbols override censorIDs")}}
    # Checking censoring
  if(!is.null(censorID)){
    if(any(sapply(censorID,function(x){!(x %in% c(0,1) | is.logical(x))}))){stop("1/0 or TRUE/FALSE expected for censorID")}
    # if(!(censorID %in% c(0,1) | is.logical(censorID))){stop("1/0 or TRUE/FALSE expected for censorID")}
    non_cen=censorID}

  if(missing(main)){plt_title="Survival function"}
  else{plt_title=main}
  
  if(missing(xlab)){x_title="t"}
  else{x_title=xlab}
  
  if(missing(ylab)){y_title="S(t)"}
  else{y_title=ylab}
  
  tmx=max(round(time))
  tmn=min(round(time))
  inc=(tmx-tmn)*0.15
  h_xdim=c(max(tmn-inc,0),tmx+inc)
  t_brk=seq(h_xdim[1],h_xdim[2],length.out = 20)
  
  #override xlim
  if(missing(xlim)){
    xdim=c(max(tmn-inc,0),tmx+inc)}
  else{xdim=xlim}
  
  #override ylim
  if(missing(ylim)){
    ydim=c(0,1)}
  else{ydim=ylim}
  
  if(is.null(group)){
    if(hist){
         hist(x=time,
         col=8,
         # xlim=xdim,   # x-axis scale
         breaks=20,
         main="Failure time distribution", # plot title
         xlab=x_title)} # axis label
    
    # Sample survival function
    if(surv_curv){
      plot(x=time,y=surv,
         xlim=xdim, 
         ylim = ydim,
         main=plt_title,pch=1,col=c("black","darkgray")[non_cen+1],
         xlab=x_title,ylab=y_title,...)}
  }
  else{
    stopifnot(is.character(group)|is.factor(group))
    n_grps=length(unique(group))
    if(length(unique(group))>3){stop("Only 3 groups may be plotted at once")}
    
    if(hist){
      # histogram of all failure times
      hist(x = time,breaks=t_brk,xlab=x_title,xlim=c(0,60),main="",col=8)
      # overlay histogram of spring failure is red
      for(i in 2:n_grps){
        hist(time[group==unique(group)[i]],breaks=t_brk,add=T,col=i)
      legend("topleft",legend=unique(group),fill=c(8,2:(n_grps)),bty="n")
      }}
      
    if(surv_curv){
      plot(x=time,y=surv,
             xlim=xdim,  
             main=plt_title, 
             xlab=x_title,ylab=y_title,cex=1.1,pch=1)
      for(i in 2:n_grps){
        points(time[group==unique(group)[i]],y=surv[group==unique(group)[i]],col=i,cex=1.1)}
      legend("bottomleft",legend=unique(group),col=c(1,2:(n_grps)),pch=1,bty="n")
      }
    }
}

