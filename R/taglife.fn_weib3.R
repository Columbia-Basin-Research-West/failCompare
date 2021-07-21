
#' Fitting 3-parameter Weibull model to failure time data (adapted from Rich's code)
#'
#' @param tags.in vector of observed time to failure (days)
#' @param model.in name of model to use.  Current option is "weibull" for 3-parameter weibull
#' @param tag.se logical for whether to compute SEs
#'
#' @return list with model objects (mod_obj), fitted values (fit_vals) and table of parameter estimates (par_tab)
#' @export
#'
taglife.fn_weib3=function(tags.in,model.in="weibull",tag.se=T){
  par_out_simp=NULL # added by SLW
  lik.in=NULL

  ######################### start loop for fitting tag-life curve after changing model or censoring tags from tail-end of tag-life study

  #### currently set to default to weibull, need to add prompt to choose model here
  if(model.in=="weibull"){lik.in=logweib3.lik; par.in=c(min(tags.in),.1,max(tags.in)) }
  stopifnot('Tag life model name not valid' =!is.null(lik.in)) #if no model picked

  # likelihood
  answer <- stats::optim(par=par.in,fn=lik.in,hessian=F,tags.in=tags.in,control = list(reltol = 1e-10,maxit=500,ndeps=rep(.01,3)))
  if(tag.se){
    answer <- stats::optim(par=answer$par,fn=lik.in,method="BFGS",hessian=T,tags.in=tags.in,control = list(parscale = abs(1/answer$par), reltol = 1e-20,maxit=500,ndeps=abs(1/answer$par)^2)) # estimate Hessian if se(tag life params wanted)
  }
  if(model.in=="weibull"){
    params.out=data.frame(matrix(NA,nrow=3,ncol=3))

    params.out[,1]=c("Beta","Gamma","Eta")
    params.out[,2]=round(abs(answer$par),4)
    names(params.out)=c("Parameter","Estimate","s.e.")

  }

  if(tag.se){
    # calculate standard errors on parameter estimates, based on Hessian
    temp = answer$hessian
    singular = colSums(temp) == 0
    var.vec = rep(0, 3)
    var.vec[!singular] = diag(solve(temp[!singular, !singular]))
    params.out[,3]=apply(cbind("(",round(sqrt(var.vec),4),")"),1,paste,collapse="")

    par_out_simp=cbind(params.out[,2],sqrt(var.vec))
    dimnames(par_out_simp)=list(c("shape","thrsh","scale"),c("params","std")) # needs to be c("params","std") to match with vitality output in fail_fit()

  }else{params.out[,3]=NA}

  ##### end of fitting loop.  need to add prompt that fit is accepted.  If not, go back to model choice
  ##### fn will return last model fit

  # fit_vals=data.frame(model="weibull3",time=c(0,tags.in),est=fail.fn_weib3(c(0,tags.in),params.out[1,2],params.out[2,2],params.out[3,2]),lcl=NA,ucl=NA)
  fit_vals=data.frame(model="weibull3",time=c(0,tags.in),est=fail_pred(c(0,tags.in),pars = params.out[,2],model = "weibull3"),lcl=NA,ucl=NA)
  mod_obj=par_out_simp
  out=list(mod_obj=mod_obj,fit_vals=fit_vals,par_tab=par_out_simp)

  return(out)
}


#' 3-parameter weibull likelihood
#'
#' @param x estimated parameters (beta, gamma, eta) or (shape,thrsh,scale)
#' @param tags.in observed time to failure
#'
#' @return log likelihood
#' @export
logweib3.lik=function(x,tags.in){
  # log mle for 3-parameter Weibull
  # INPUT
  # x: estimated parameters (beta, gamma, eta)
  # tags.in: observed time to failure (days)
  # OUTPUT
  # negative log-likelihood
  x=abs(x)
  b.in=x[1] # beta
  g.in=x[2] # gamma
  n.in=x[3] # eta
  if (sum(g.in>tags.in)>0){g.in=min(tags.in)-.0000001}

  l.out=log(b.in)-(b.in*log(n.in))+((b.in-1)*log(tags.in-g.in))-(((tags.in-g.in)/n.in)^b.in)
  l.out=-sum(l.out)
  return(l.out)
}

