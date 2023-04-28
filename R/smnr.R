##initial value
#' Semiparametric Mixture Data Generator
#'
#'Generate data from two-component semiparametric mixture of regression models with
#' m(x) functions:\cr \eqn{m_1(x) = 4 - sin(2\pi x)} and \eqn{m_2(x) = 1.5 + cos(3\pi x)}, used in the examples for
#' backfitlocal,backfitglobal etc.
#' @param n number of observations
#' @param p probability of an observation in first component
#' @param var variance of the observations
#' @param u grid point
#'
#' @return list of x,y, true mean functions and true interpolated mean functions.
#' @export
#'
#' @examples
#' n=100;
#' u=seq(from=0,to=1,length=100);
#' true_p=c(0.3,0.7);
#' true_var=c(0.09,0.16);
#' out=gen_mixreg(n,true_p[1],true_var,u);
gen_mixreg<-function(n,p,var,u){

  n1=rbinom(1,n,p);#rbinom(n,size,prob)，n生成随机数的数量，size伯努利试验的次数，prob一次试验成功的概率。
  x=runif(n,0,1);

  err=numeric(n);
  err[1:n1]=rnorm(n1,0,1)*sqrt(var[1]);
  err[(n1+1):n]=rnorm(n-n1,0,1)*sqrt(var[2]);
  true=numeric(n);
  true[1:n1]=4-sin(2*pi*x[1:n1]);
  true[(n1+1):n]=1.5+cos(3*pi*x[(n1+1):n]);

  y=true+err;
  true_mu=matrix(c(4-sin(2*pi*x),1.5+cos(3*pi*x)),nrow=n);
  true_mu_u=matrix(c(4-sin(2*pi*u),1.5+cos(3*pi*u)),nrow=length(u));

  out=list(x=x,y=y,true_mu=true_mu,true_mu_u=true_mu_u)
  out

}
gen_mixreg_boot<-function(x,phat,muhat,varhat){
  n=dim(muhat)[1];
  k=dim(muhat)[2];
  p=round(phat[1]*100)/100;
  n0=rbinom(1,n,p[1])
  n1=c(n0,n-n0);
  pos=1;
  err=matrix(rep(0,n),nrow=n);
  true=err;

  for(j in 1:k){
    ind=pos:sum(n1[1:j]);
    err[ind]=rnorm(n1[j],mean=0,sd=1)*sqrt(varhat[j]);
    true[ind]=muhat[ind,j];
    pos=sum(n1[1:j])+1;
  }

  y_boot=true+err;
}



#true=list(true_p=true_p,true_mu=true_mu,true_var=true_var)
#ini=list(est_p=est_p,est_mu=est_mu,est_var=est_var)


#' Semiparametric Mixtures of Nonparametric Regressions with Local EM-type Algorithm
#'
#'A new class of semiparametric mixture of regression models,
#'where the mixing proportions and variances are constants,but the component regression functions are smooth functions of a covariate.
#'A onestep backfitting estimate and local EM-type algorithms have been proposed to achieve
#'the optimal convergence rate for both the global parameters and the nonparametric regression functions.
#'
#' @param x independent variable x
#' @param y dependent variable y
#' @param u  vector of grid points
#' @param h bandwidth for nonparametric regression; If NULL, will be calculated by reliable
#' and extremely fast kernel density estimator for one-dimensional data (ref: Z.I.Botev,J.F.Grotowski and D.P.Kroese,"KERNEL DENSITY ESTIMATION VIA DIFFUSION" ,Submitted to the Annals of Statistics, 2009)
#' @param ini initial values of parameters;ini=list(est_p=est_p,est_mu=est_mu,est_var=est_var); If NULL, estimated automatically by regression spline approximation.
#' @param true true value of the parameters;\cr true=list(true_p=true_p,true_mu=true_mu,true_var=true_var); See details.
#'
#' @return
#' est_p: estimated proportions; est_mu: estimated mean functions ,est_mu_u: estimated mean functions by linear interpolating,est_var: estimated variance,lh:likelihood;
#' initialtrue: which estimation was returned: see details.
#'@details
#'True parameter: If specified, likelihood of the estimated parameters will be calculated and compare
#' with likelihood calculated from initial value. If likelihood from true value is greater than that from the initial values
#' then return estimations from true value with initialtrue=1, otherwise return estimations from initial values with initialtrue=0.
#' @export
#'
#' @examples
#'
#'#produce data that matches the description using gen_mixreg function
#'#true_mu=(4-sin(2*pi*x),1.5+cos(3*pi*x))
#'n=100;
#'u=seq(from=0,to=1,length=100);
#'true_p=c(0.3,0.7);
#'true_var=c(0.09,0.16);
#'out=gen_mixreg(n,true_p[1],true_var,u);
#'
#'x=out$x;
#'y=out$y;
#'true_mu=out$true_mu;
#'true=list(true_p=true_p,true_mu=true_mu,true_var=true_var)
#'
#'#estimate parameters using backfitlocal function.
#'\donttest{esttrue=backfitlocal(x,y,u,h=0.08,true=true)}
#'\donttest{est=backfitlocal(x,y,u,h=NULL)}
backfitlocal<-function(x,y,u = NULL,h = NULL,ini = NULL,true = NULL){
  n=length(y)#added sep18 2022 Xin Shen to fix the error message of 'no visible binding for global variable ‘n’'
  if(is.null(u)){
    u = seq(from=min(x),to=max(x),length=100);
  }
  if(is.null(h)){
    h = kdebw(x,2^14);
  }

  if(is.null(ini)){
    IDX=kmeans(cbind(x,y),2)$cluster;
    ini_p = rbind(2-IDX,IDX-1);
    ini = list(p=t(ini_p))
    ini = mixbspline(u,x,y,ini);
  }
  est_mu_x = ini$est_mu;
  est_p = ini$est_p; est_p=cbind(rep(est_p[1],n),rep(est_p[2],n));
  est_var = ini$est_var; est_var=cbind(rep(est_var[1],n),rep(est_var[2],n));

  est1=backfitting_step1(x,y,u,est_p,est_mu_x,est_var,h);
  est_p=est1$est_p;
  est_var=est1$est_var;
  est_mu_x=est1$est_mu_x;
  est2=backfitting_step2(x,y,u,est_p,est_mu_x,est_var);
  est_p=est2$est_p;
  est_var=est2$est_var;
  est3=backfitting_step3(x,y,u,est_p,est_mu_x,est_var,h);
  est_mu_x=est3$est_mu_x;
  est_mu_u=est3$est_mu_u;
  likelihood=est3$lh;

  res=list(est_p=est_p,est_mu_x=est_mu_x,est_var=est_var,est_mu_u=est_mu_u,lh=likelihood)
  out=list(est_p=est_p,est_mu=est_mu_x,est_var=est_var,est_mu_u=est_mu_u,lh=likelihood)

  if(is.null(true)){
    out = out
  }else{
    est_mu_x = true$true_mu;
    est_p = true$true_p; est_p = cbind(rep(est_p[1],n),rep(est_p[2],n));
    est_var = true$true_var; est_var = cbind(rep(est_var[1],n),rep(est_var[2],n));

    est1=backfitting_step1(x,y,u,est_p,est_mu_x,est_var,h);
    est_p=est1$est_p;
    est_var=est1$est_var;
    est_mu_x=est1$est_mu_x;
    est2=backfitting_step2(x,y,u,est_p,est_mu_x,est_var);
    est_p=est2$est_p;
    est_var=est2$est_var;
    est3=backfitting_step3(x,y,u,est_p,est_mu_x,est_var,h);
    est_mu_x=est3$est_mu_x;
    est_mu_u=est3$est_mu_u;
    likelihood=est3$lh;

    if(res$lh<likelihood){
      out=list(est_p=res$est_p,est_mu=res$est_mu_x,est_mu_u=res$est_mu_u,est_var=res$est_var,lh=res$lh,initialtrue=0)
    }else{
      out=list(est_p=est_p,est_mu=est_mu_x,est_mu_u=est_mu_u,est_var=est_var,lh=likelihood,initialtrue=1)
    }
  }

  return(out)
}



backfitting_step1<-function(x,y,u,est_p,est_mu_x,est_var,h){

  n = dim(est_mu_x)[1];
  numc = dim(est_mu_x)[2];
  m = length(u);

  r = matrix(numeric(n*numc),nrow=n);
  f = r;
  est_p_u = matrix(numeric(m*numc),nrow=m);
  est_mu_u = matrix(numeric(m*numc),nrow=m);
  est_var_u = matrix(numeric(m*numc),nrow=m);

  rsd=1;likelihood=1;iter=0;acc=10^(-3);

  repeat
  {
    iter = iter+1;

    #E-step
    for(i in 1:numc){
      f[,i] = dnorm(y,est_mu_x[,i],est_var[,i]^0.5);
    }

    #M-step
    for(i in 1:numc){
      r[,i] = est_p[,i]*f[,i]/apply(f*est_p,1,sum);
      # W = exp(-(((matrix(rep(t(x),m),nrow=m)-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));
      W = exp(-(((t(matrix(rep(x,m),nrow=m))-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));#Updated sep20 2022 after meeting with Yan Ge.
      R = t(matrix(rep(r[,i],m),ncol=m));
      RY = t(matrix(rep(r[,i]*y,m),ncol=m));
      est_p_u[,i] = apply(W*R,1,sum)/apply(W,1,sum);
      est_p[,i] = approx(u,est_p_u[,i],x)$y;#approx:线性插值
      est_mu_u[,i] = apply(W*RY,1,sum)/apply(W*R,1,sum);
      est_mu_x[,i] = approx(u,est_mu_u[,i],x)$y;
      RE = (t(matrix(rep(y,m),ncol=m))-matrix(rep(est_mu_u[,i],n),ncol=n))^2;
      est_var_u[,i] = apply(W*RE*R,1,sum)/apply(W*R,1,sum);
      est_var[,i] = approx(u,est_var_u[,i],x)$y;
    }

    rsd = likelihood-sum(log(apply(f*est_p,1,sum)));
    likelihood = sum(log(apply(f*est_p,1,sum)));
    if(abs(rsd)<acc|iter>200){break}
  }

  out=list(est_p=est_p,est_mu_x=est_mu_x,est_var=est_var)
  out

}

backfitting_step2<-function(x,y,u,est_p,est_mu_x,est_var){

  n = dim(est_mu_x)[1];
  numc = dim(est_mu_x)[2];
  m = length(u);

  r = matrix(numeric(n*numc),nrow=n);
  f = r;

  rsd=1;likelihood=1;iter=0;acc=10^(-3);
  ###
  repeat
  {
    iter = iter+1;

    #E-step
    for(i in 1:numc){
      f[,i] = dnorm(y,est_mu_x[,i],est_var[,i]^0.5);
    }

    #M-step
    for(i in 1:numc){
      r[,i] = est_p[,i]*f[,i]/apply(f*est_p,1,sum);
      est_p[,i] = matrix(rep(sum(r[,i])/n,n),nrow=n);
      est_var[,i] = matrix(rep(sum(r[,i]*(y-est_mu_x[,i])^2)/sum(r[,i]),n),nrow=n);
    }

    rsd = likelihood-sum(log(apply(f*est_p,1,sum)));
    likelihood = sum(log(apply(f*est_p,1,sum)));
    if(abs(rsd)<acc|iter>200){break}

  }

  out=list(est_p=est_p,est_var=est_var)
  out

}

backfitting_step3<-function(x,y,u,est_p,est_mu_x,est_var,h){

  n = dim(est_mu_x)[1];
  numc = dim(est_mu_x)[2];
  m = length(u);

  r = matrix(numeric(n*numc),nrow=n);
  f = r;est_mu_u = matrix(numeric(m*numc),nrow=m);

  rsd=1;likelihood=1;iter=0;acc=10^(-3);

  repeat
  {
    iter = iter+1;

    #E-step
    for(i in 1:numc){
      f[,i] = dnorm(y,est_mu_x[,i],est_var[,i]^0.5);

    }

    #M-step
    for(i in 1:numc){
      r[,i] = est_p[,i]*f[,i]/apply(f*est_p,1,sum);
      # W = exp(-(((matrix(rep(t(x),m),nrow=m)-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));
      W = exp(-(((t(matrix(rep(x,m),nrow=m))-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));#Updated sep20 2022 after meeting with Yan Ge.
      R = t(matrix(rep(r[,i],m),ncol=m));
      RY = t(matrix(rep(r[,i]*y,m),ncol=m));
      est_mu_u[,i] = apply(W*RY,1,sum)/apply(W*R,1,sum);
      est_mu_x[,i] = approx(u,est_mu_u[,i],x)$y;
    }

    rsd = likelihood-sum(log(apply(f*est_p,1,sum)));
    likelihood = sum(log(apply(f*est_p,1,sum)));
    if(abs(rsd)<acc|iter>200){break}
  }

  out=list(est_mu_x=est_mu_x,est_mu_u=est_mu_u,lh=likelihood)
  out

}




#' Semiparametric Mixtures of Nonparametric Regressions with Global EM-type Algorithm
#'
#'
#'A new class of semiparametric mixture of regression models,
#'where the mixing proportions and variances are constants,but the component regression functions are smooth functions of a covariate.
#'A onestep backfitting estimate and global EM-type algorithms have been proposed to achieve
#'the optimal convergence rate for both the global parameters and the nonparametric regression functions.
#'
#' @param x independent variable x
#' @param y dependent variable y
#' @param u  vector of grid points
#' @param h bandwidth for nonparametric regression; If NULL, will be calculated by reliable
#' and extremely fast kernel density estimator for one-dimensional data (ref: Z.I.Botev,J.F.Grotowski and D.P.Kroese,"KERNEL DENSITY ESTIMATION VIA DIFFUSION" ,Submitted to the Annals of Statistics, 2009)
#' @param ini initial values of parameters;ini=list(est_p=est_p,est_mu=est_mu,est_var=est_var); If NULL, estimated automatically by regression spline approximation.
#' @param true true value of the parameters; \cr true=list(true_p=true_p,true_mu=true_mu,true_var=true_var); See details.
#'
#' @return
#' est_p: estimated proportions; est_mu: estimated mean functions,est_mu_u: estimated mean functions with linear interpolating,est_var: estimated variance,lh:likelihood;
#' initialtrue:which estimation was returned: see details.
#'@details
#'True parameter: If specified, likelihood of the estimated parameters will be calculated and compare
#' with likelihood calculated from initial value. If likelihood from true value is greater than that from the initial values
#' then return estimations from true value with initialtrue=1, otherwise return estimations from initial values with initialtrue=0.
#' @export
#' @examples
#'#produce data that matches the description using gen_mixreg function
#'#true_mu=c(4-sin(2*pi*x),1.5+cos(3*pi*x))
#'n=100;
#'u=seq(from=0,to=1,length=100);
#'true_p=c(0.3,0.7);
#'true_var=c(0.09,0.16);
#'out=gen_mixreg(n,true_p[1],true_var,u);
#'
#'x=out$x;
#'y=out$y;
#'true_mu=out$true_mu;
#'true=list(true_p=true_p,true_mu=true_mu,true_var=true_var)
#'
#'#estimate parameters using backfitglobal function.
#'esttrue=backfitglobal(x,y,u,h=0.08,true=true)
#'est=backfitglobal(x,y,u,h=NULL)
backfitglobal<-function(x,y,u = NULL,h = NULL,ini = NULL,true = NULL){
  if(is.null(u)){
    u = seq(from=min(x),to=max(x),length=100);
  }
  if(is.null(h)){
    h = kdebw(x,2^14);
  }

  if(is.null(ini)){
    IDX=kmeans(cbind(x,y),2)$cluster;
    ini_p = rbind(2-IDX,IDX-1);
    ini = list(p=t(ini_p))
    ini = mixbspline(u,x,y,ini);
  }
  n=length(y)#added sep18 2022 Xin Shen to fix the error message of 'no visible binding for global variable ‘n’'
  est_mu_x = ini$est_mu;
  est_p = ini$est_p; est_p=cbind(rep(est_p[1],n),rep(est_p[2],n));
  est_var = ini$est_var; est_var=cbind(rep(est_var[1],n),rep(est_var[2],n));

  n = dim(est_mu_x)[1];
  numc = dim(est_mu_x)[2];
  m = length(u);

  r = matrix(numeric(n*numc),nrow=n);
  f = r;est_mu_u = matrix(numeric(m*numc),nrow=m);

  rsd=1;likelihood=1;iter=0;acc=10^(-3);

  repeat
  {
    iter = iter+1;

    #E-step
    for(i in 1:numc){
      f[,i] = dnorm(y,est_mu_x[,i],est_var[,i]^0.5);
    }

    #M-step
    for(i in 1:numc){
      r[,i] = est_p[,i]*f[,i]/apply(f*est_p,1,sum);
      est_p[,i] = matrix(rep(sum(r[,i])/n,n),nrow=n);
      est_var[,i] = matrix(rep(sum(r[,i]*(y-est_mu_x[,i])^2)/sum(r[,i]),n),nrow=n);
      # W = exp(-(((matrix(rep(t(x),m),nrow=m)-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));
      W = exp(-(((t(matrix(rep(x,m),nrow=m))-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));#Updated sep20 2022 after meeting with Yan Ge.
      R = t(matrix(rep(r[,i],m),ncol=m));
      RY = t(matrix(rep(r[,i]*y,m),ncol=m));
      est_mu_u[,i] = apply(W*RY,1,sum)/apply(W*R,1,sum);
      est_mu_x[,i] = approx(u,est_mu_u[,i],x)$y;
    }

    rsd = likelihood-sum(log(apply(f*est_p,1,sum)));
    likelihood = sum(log(apply(f*est_p,1,sum)));
    if(abs(rsd)<acc|iter>200){break}
  }

  res=list(est_p=est_p,est_mu_x=est_mu_x,est_var=est_var,est_mu_u=est_mu_u,lh=likelihood)
  out=list(est_p=est_p,est_mu=est_mu_x,est_var=est_var,est_mu_u=est_mu_u,lh=likelihood)

  if(is.null(true)){
    out = out
  }else{
    est_mu_x = true$true_mu;
    est_p = true$true_p; est_p = cbind(rep(est_p[1],n),rep(est_p[2],n));
    est_var = true$true_var; est_var = cbind(rep(est_var[1],n),rep(est_var[2],n));

    n = dim(est_mu_x)[1];
    numc = dim(est_mu_x)[2];
    m = length(u);

    r = matrix(numeric(n*numc),nrow=n);
    f = r;est_mu_u = r;

    rsd=1;likelihood=1;iter=0;acc=10^(-3);

    repeat
    {
      iter = iter+1;

      #E-step
      for(i in 1:numc){
        f[,i] = dnorm(y,est_mu_x[,i],est_var[,i]^0.5);
      }

      #M-step
      for(i in 1:numc){
        r[,i] = est_p[,i]*f[,i]/apply(f*est_p,1,sum);
        est_p[,i] = matrix(rep(sum(r[,i])/n,n),nrow=n);
        est_var[,i] = matrix(rep(sum(r[,i]*(y-est_mu_x[,i])^2)/sum(r[,i]),n),nrow=n);
        W = exp(-(((t(matrix(rep(t(x),m),nrow=m))-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));
        R = t(matrix(rep(r[,i],m),ncol=m));
        RY = t(matrix(rep(r[,i]*y,m),ncol=m));
        est_mu_u[,i] = apply(W*RY,1,sum)/apply(W*R,1,sum);
        est_mu_x[,i] = approx(u,est_mu_u[,i],x)$y;
      }

      rsd = likelihood-sum(log(apply(f*est_p,1,sum)));
      likelihood = sum(log(apply(f*est_p,1,sum)));
      if(abs(rsd)<acc|iter>200){break}
    }

    if(res$lh<likelihood){
      out=list(est_p=res$est_p,est_mu=res$est_mu_x,est_mu_u=res$est_mu_u,est_var=res$est_var,lh=res$lh,initialtrue=0)
    }else{
      out=list(est_p=est_p,est_mu=est_mu_x,est_mu_u=est_mu_u,est_var=est_var,lh=likelihood,initialtrue=1)
    }
  }

  return(out)
}



#'  Mixture of Nonparametric Regressions Using B-spline
#'
#'This function is used to fit mixture of nonparametric regressions using b-spline
#'assuming the error is normal and the variance of the errors are equal
#' @param u vector of grid points.
#' @param x the predictors we observed; needs to be standardized;
#' @param y the observation we observed
#' @param ini initial values of proportion. If left as NULL will be calculated from mixture of the linear regression
#' with normal error assumption.
#' @param k number of component, the default value is 2;
#' @param mm order of spline, default is cubic spline with mm=4;
#' @param tt the knots; default is quantile(x,seq(0.05,0.95,0.2))
#' @param a the interval that contains x; default is c(mean(x)-2*sqrt(var(x)),mean(x)+2*sqrt(var(x)));
#'
#' @return mu_u:estimated mean functions with interpolate;est_mu: estimated mean functions,est_p: estimated proportion;est_var: estimated variance
#'
mixbspline<-function(u,x,y,ini = NULL,k = NULL,mm = NULL,tt = NULL,a = NULL){
  if(is.null(k)){
    k = 2;
  }
  if(is.null(mm)){
    mm = 4;
  }
  if(is.null(tt)){
    tt = quantile(x,seq(0.05,0.95,0.2));
  }
  if(is.null(a)){
    a = c(mean(x)-2*sqrt(var(x)),mean(x)+2*sqrt(var(x)));
  }
  if(is.null(ini)){
    ini = mixlin_smnr(x,y,k);
  }

  t = c(rep(a[1],mm),rep(a[2],mm));
  K = length(tt);
  BX = bspline_basis(1,mm,t,x);
  BXu = bspline_basis(1,mm,t,u);
  for(i in 2:(mm+K)){
    BX = cbind(BX,bspline_basis(1,mm,t,x));
    BXu = cbind(BXu,bspline_basis(1,mm,t,u));
  }

  fit = mixlin_smnr(BX,y,k,ini,0);

  out=list(mu_u=BXu%*%fit$beta,est_mu=BX%*%fit$beta,est_p=fit$pi,est_var=(fit$sigma)^2)
  return(out)
}



#is used to fit mixture of the linear regression
#assuming the error is normal

mixlin_smnr<-function(x,y,k = NULL,ini = NULL,intercept = NULL){
  n = length(y); x = as.matrix(x);
  if(dim(x)[1]<dim(x)[2]){
    x = t(x); y = t(y);
  }

  if(is.null(k)){
    k = 2;
  }
  if(is.null(intercept)){
    intercept = 1;
  }
  if(intercept==1){
    X = cbind(rep(1,n),x);
  }else{
    X = x;
  }
  if(is.null(ini)){
    ini = mixlin_smnrveq(x,y,k);
  }
  p = ini$p;

  stop = 0; numiter = 0; dif = 10; difml = 10; ml = 10^8; dif1=1;
  beta = matrix(rep(0,dim(X)[2]*k),ncol=k); e =  matrix(rep(0,n*k),nrow=n); sig = numeric();

  for(j in 1:k){
    temp = t(X)%*%diag(p[,j]);
    beta[,j] = MASS::ginv(temp%*%X+10^(-5)*diag(dim(X)[2]))%*%temp%*%y;
    e[,j] = y-X%*%beta[,j];
    sig[j] = sqrt(p[,j]%*%(e[,j]^2)/sum(p[,j]));
  }
  prop = apply(p,2,sum)/n;

  acc=10^(-3)*mean(mean(abs(beta)));

  while(stop==0){
    denm =  matrix(rep(0,n*k),nrow=n);
    for(j in 1:k){
      denm[,j] = dnorm(e[,j],0,sig[j]);
    }
    f = apply(prop*denm,1,sum); prebeta=beta;

    numiter = numiter+1; lh[numiter] = sum(log(f));
    if(numiter > 2){
      temp = lh[numiter-1]-lh[numiter-2];
      if(temp == 0){
        dif = 0;
      }else{
        c = (lh[numiter]-lh[numiter-1])/temp;#convergence rate of Em algorithm
        if(c>0 && c<1){
          preml = ml; ml=lh[numiter-1]+(lh[numiter]-lh[numiter-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
          dif = ml-lh[numiter]; difml = abs(preml-ml);
        }
      }
    }

    #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
    if(dif<0.001 && difml<10^(-6) || numiter>5000 || dif1<acc){
      stop = 1;
    }else{
      p = prop*denm/matrix(rep(f,k),ncol=k);
      for(j in 1:k){
        temp = t(X)%*%diag(p[,j]);
        beta[,j] = MASS::ginv(temp%*%X)%*%temp%*%y;
        e[,j] = y-X%*%beta[,j];
        sig[j] = sqrt(t(p[,j])%*%(e[,j]^2)/sum(p[,j]));
      }
      prop = apply(p,2,sum)/n; dif1 = max(abs(beta-prebeta));
    }
    if(min(sig)/max(sig)<0.001 || min(prop)<1/n || is.nan(min(beta[,1]))){
      ini=list(beta=beta,sig=sig,prop=prop)
      res = mixlin_smnrveq(x,y,k,ini,0);
      out=list(beta=res$beta,sig=res$sig,pi=res$pi,p=p,id=0,lh=lh[numiter])
    }else{
      out=list(beta=beta,sigma=sig,pi=prop,p=p,id=1,lh=lh[numiter])
    }

  }
  #out=list(beta=beta,sigma=sig,pi=prop,p=p,id=1,lh=lh[numiter])

  return(out)
}



#is used to fit mixture of the linear regression
#assuming the error is normal and the variance of the errors are equal
mixlin_smnrveq<-function(x,y,k = NULL,ini = NULL,intercept = NULL){
  n = length(y);
  if(dim(x)[1]<dim(x)[2]){
    x = t(x); y = t(y);
  }

  if(is.null(k)){
    k = 2;
  }
  if(is.null(intercept)){
    intercept = 1;
  }
  if(intercept==1){
    X = cbind(rep(1,n),x);
  }else{
    X = x;
  }
  if(is.null(ini)){
    ini = mixmnorm(cbind(x,y),k,0);
    prop = ini$pi; p = t(ini$p);
    beta = matrix(rep(0,dim(X)[2]*k),ncol=k); e =  matrix(rep(0,n*k),nrow=n); sig = numeric();

    for(j in 1:k){
      temp = t(X)%*%diag(p[,j]);
      beta[,j] = MASS::ginv(temp%*%X)%*%temp%*%y;
      e[,j] = y-X%*%beta[,j];
      sig[j] = sqrt(t(p[,j])%*%(e[,j]^2)/sum(p[,j]));
    }
    sig = rep(sqrt(sum(e^2*p)/n),k);
    ini$beta = beta; ini$sigma = sig;
  }else{
    beta = ini$beta; sig = ini$sig; prop = ini$prop;
    denm = matrix(rep(0,n*k),nrow=n);
    for(j in 1:k){
      e[,j] = y-X%*%beta[,j];
      denm[,j] = dnorm(e[,j],0,sig[j]);
    }
    f = apply(prop*denm,1,sum);
    p = prop*denm/matrix(rep(f,k),ncol=k);
  }

  stop = 0; numiter = 0; dif = 10; difml = 10; ml = 10^8;

  while(stop==0){
    denm = matrix(rep(0,n*k),nrow=n);
    for(j in 1:k){
      denm[,j] = dnorm(e[,j],0,sig[j]);
    }
    f = apply(prop*denm,1,sum);
    numiter = numiter+1; lh[numiter] = sum(log(f));
    if(numiter > 2){
      temp = lh[numiter-1]-lh[numiter-2];
      if(temp == 0){
        dif = 0;
      }else{
        c = (lh[numiter]-lh[numiter-1])/temp;#convergence rate of Em algorithm
        if(c>0 && c<1){
          preml = ml; ml=lh[numiter-1]+(lh[numiter]-lh[numiter-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
          dif = ml-lh[numiter]; difml = abs(preml-ml);
        }
      }
    }

    #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
    if(dif<0.001 && difml<10^(-6) || numiter>5000){
      stop = 1;
    }else{
      p = prop*denm/matrix(rep(f,k),ncol=k);
      for(j in 1:k){
        temp = t(X)%*%diag(p[,j]);
        beta[,j] = MASS::ginv(temp%*%X)%*%temp%*%y;
        e[,j] = y-X%*%beta[,j];
      }
      prop = apply(p,2,sum)/n; sig = rep(sqrt(sum(e^2*p/n)),k);
    }
    if(min(prop)<1/n || is.nan(min(beta[,1]))){
      out = ini;
    }else{
      out=list(beta=beta,sigma=sig,pi=prop,p=p,ind=1)
    }
  }

  #out=list(beta=beta,sigma=sig,pi=prop,p=p,ind=1)

  return(out)
}


bspline_basis<-function(i,m,t,x){
  x = as.matrix(x,nrow=length(x));
  B = matrix(rep(0,dim(x)[1]*dim(x)[2]),nrow=dim(x)[1]);
  if(m>1){
    b = bspline_basis(i,m-1,t,x);
    dn = x-t[i];
    dd = t[i+m-1]-t[i];
    if(dd!=0){
      B = B+b*(dn/dd);
    }
    b = bspline_basis(i+1,m-1,t,x);
    dn = t[i+m]-x;
    dd = t[i+m]-t[i+1];
    if(dd!=0){
      B = B+b*(dn/dd);
    }
  }else{
    # B = t[i]<=x[1]&&x[1]<t[i+1]
    B = (t[i]<=x)&(x<t[i+1]) #Xin Shen updated on sep23 2022 after checking with Yan Ge
  }
  return(B)
}






kdebw<-function(y,n = NULL,MIN = NULL,MAX = NULL){
  if (is.null(n)) {
    n = 2^12;
  }
  n = 2^ceiling(log2(n));

  if(is.null(MIN) || is.null(MAX)){
    minimum = min(y); maximum = max(y);
    Range = maximum-minimum;
    MIN = minimum-Range/10; MAX = maximum+Range/10;
  }

  #set up the grid over which the density estimate is computed;
  R = MAX-MIN; dx = R/(n-1); xmesh = MIN+seq(from=0,to=R,by=dx); N = length(y);

  #bin the data uniformly using the grid define above;
  initial_data = graphics::hist(y,xmesh,plot=FALSE)$counts/N;
  a = dct1d(initial_data);#discrete cosine transform of initial data
  I = as.matrix(c(1:(n-1))^2); a2 = as.matrix((a[2:length(a)]/2)^2);

  #ft1=function(t){
  #  y=fixed_point(abs(t),I,a2,N)-abs(t);
  #}
  #t_star=abs(stats::uniroot(ft1,c(-0.1,1)));
  ft2<-function(t){
    y = abs(fixed_point(t,I,a2,N)-t);
  }
  #if(is.null(t_star) || is.infinite(t_star)){
  t_star = optimize(ft2,c(0,0.1))$minimum;
  #}

  bandwidth = sqrt(t_star)*R;

  return(bandwidth)
}

fixed_point<-function(t,I,a2,N){
  functional<-function(df,s){
    K0 = prod(seq(1,2*s-1,2))/sqrt(2*pi);
    const = (1+(1/2)^(s+1/2))/3;
    t = (2*const*K0/N/df)^(2/(3+2*s));
    f = 2*pi^(2*s)*sum(I^s*a2*exp(-I*pi^2*t));
    return(f)
  }

  f=numeric()
  f[5] = 2*pi^10*sum(I^5*a2*exp(-I*pi^2*t));
  for(s in seq(4,2,-1)){
    f[s] = functional(f[s+1],s);
  }
  time = (2*N*sqrt(pi)*f[2])^(-2/5);
  out = (t-time)/time;

  return(out)
}

dct1d<-function(data){
  data = c(data,0);
  n = length(data);
  weight = c(1,2*exp(-1i*seq(1,(n-1),1)*pi/(2*n)));
  data = c(data[seq(1,n,2)],data[seq(n,2,-2)]);
  out = Re(weight*fft(data));

  return(out)
}


