#' Varpop function Data Generator
#'
#' Used to generate data for MixReg_Pvary function examples.
#' @param n number of observations
#' @param numc number of components
#'
#' @return
#' list of x,y, true mean functions and true p.
#'
#'@export
#' @examples
#' gen_mixreg1(n=100,numc=2)
gen_mixreg1<-function(n,numc){
  x=runif(n,0,1);
  true_p=matrix(c(0.1+0.8*sin(pi*x),1-0.1-0.8*sin(pi*x)),ncol=2);

  cum_p=matrix(c(0.1+0.8*sin(pi*x),rep(1,n)),ncol=2);
  gen_p=runif(n,0,1);
  bernulli=matrix(rep(0,n*numc),nrow=n);
  bernulli[,1]=(gen_p<cum_p[,1]);
  for(i in 2:numc){
    bernulli[,i]=(cum_p[,i-1]<gen_p)*(gen_p<cum_p[,i]);
  }

  e1=0.3;e2=0.4;
  func1=4-2*x+e1*rnorm(n,0,1);
  func2=3*x+e2*rnorm(n,0,1);
  func=matrix(c(func1,func2),ncol=2);
  y=apply(func*bernulli,1,sum);

  out=list(x=x,y=y,true_p=true_p);
  return(out)
}



#' Estimate Mixture of Regression Models with Varying Mixing Proportions
#'
#'Return EM algorithm output for mixture of regression models with varying mixing proportions.
#'The varying proportions are estimated by local constant method.
#' @param x is an nxp matrix of predictors, intercept will be automatically added to x.
#' @param y n vector dimensional of response values
#' @param numc number of components, if NULL,numc = 2
#' @param z is a variable that proportion varying with. It is a univariate vector and can be part of x
#' @param u grid points from z, default is NULL, and 100 equally spaced grid points will
#' automatically generated using the minimum and maximum of z.
#' @param h bandwidth for kernel
#' @param kernel kernel for local constant method, 1 = gaussian kernel, 2 = epanechnikov
#' @param ini list of initial values:
#'  list(est_p=est_p,est_beta=est_beta,est_var=est_var),
#'  est_p: numc values, est_beta: numc*(p+1), est_var: numc
#' @param true list of true values: list(true_p=true_p,true_beta=true_beta,true_var=true_var) true value
#'
#' @return
#' est_p_z : the estimated mixing proportions
#' est_beta : the estimated regression coefficients
#' est_var : the estimated global variance
#' lh : the estimated log-likelihood
#' initialtrue : whether there is a true value, 1 true, 0 false
#' @export
#'
#'@importFrom mixtools regmixEM
#'
#' @examples
#' n=100;
#' numc=2;
#' u=seq(from=0,to=1,length=100);
#' true_beta=cbind(c(4,-2),c(0,3));
#' true_var=c(0.09,0.16);
#' data = gen_mixreg1(n,numc);
#' x = data$x;
#' y = data$y;
#' true_p = data$true_p;
#' true=list(true_p=true_p,true_beta=true_beta,true_var=true_var)
#' est=MixReg_Pvary(x,y,numc,z=x,u,h=0.08,true=true)
#'
#'
#'
MixReg_Pvary<-function(x,y,numc = NULL,z = NULL,u = NULL,h = NULL,kernel = NULL,ini = NULL,true = NULL){
  n=length(y)
  X = cbind(rep(1,n),x);

  if(is.null(numc)){
    numc = 2;
  }
  if(is.null(z)){
    z = X[,2];
  }
  if(is.null(u)){
    u = seq(from=min(z),to=max(z),length=100);
  }
  if(is.null(h)){
    h = kdebw(x,2^14);
  }
  if(is.null(kernel)){
    kernel = 1;
  }

  #if(is.null(ini)){
  #  ini = emest_linear(x,y,numc)
  #}

  if(is.null(ini)){
    ini = regmixEM(y,x)
    est_p_z = ini$lambda; est_p_z = cbind(rep(est_p_z[1],n),rep(est_p_z[2],n));
    est_beta = ini$beta;
    est_var = ini$sigma; est_var=est_var^2;
  }else{
    est_p_z = ini$est_p; est_p_z = cbind(rep(est_p_z[1],n),rep(est_p_z[2],n));
    est_beta = ini$est_beta;
    est_var = ini$est_var;
  }


  n = length(y);
  m = length(u);

  est_p_u = matrix(rep(0,m*numc),nrow=m);
  f = matrix(rep(0,n*numc),nrow=n); r = f;

  rsd=1;likelihood=1;iter=0;acc=10^(-5);
  repeat
  {
    iter=iter+1;

    for(i in 1:numc){
      f[,i] = dnorm(y,X%*%est_beta[,i],est_var[i]^0.5*rep(1,n));
    }

    for(i in 1:numc){
      r[,i] = est_p_z[,i]*f[,i]/apply(f*est_p_z,1,sum);
      s = sqrt(r[,i]);
      yy = y*s;xx = X*matrix(rep(s,2),nrow=n);
      b = coefficients(lm(yy~xx));
      #est_beta[,i] = c(b[2],b[3]);
      est_beta[,i] = b[-1];
      est_var[i] = sum(r[,i]*(y-X%*%est_beta[,i])^2)/sum(r[,i]);

      #W = exp(-(((matrix(rep(t(z),m),nrow=m)-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));
      #W = 0.75*(1-(matrix(rep(t(z),m),nrow=m)-matrix(rep(u,n),ncol=n)))/h[i];W[W<0] = 0;
      #R = t(matrix(rep(r[,i],m),ncol=m));
      #est_p_u[,i] = apply(W*R,1,sum)/apply(W,1,sum);

      est_p_u[,i] = local_constant(r[,i],z,u,h,kernel);
      est_p_z[,i] = approx(u,est_p_u[,i],z)$y;
    }

    rsd = likelihood-sum(log(apply(f*est_p_z,1,sum)));
    likelihood = sum(log(apply(f*est_p_z,1,sum)));
    if(abs(rsd)<acc|iter>200){break}
  }

  res=list(est_p_u=est_p_u,est_p_z=est_p_z,est_beta=est_beta,est_var=est_var,lh=likelihood)
  out=list(est_p_u=est_p_u,est_p_z=est_p_z,est_beta=est_beta,est_var=est_var,lh=likelihood)

  if(is.null(true)){
    out = out
  }else{
    est_p_z = true$true_p;
    est_beta = true$true_beta;
    est_var = true$true_var;
    #est_p_z = cbind(rep(est_p_z[1],n),rep(est_p_z[2],n));

    n = length(y);
    m = length(u);

    est_p_u = matrix(rep(0,m*numc),nrow=m);
    f = matrix(rep(0,n*numc),nrow=n); r = f;

    rsd=1;likelihood=1;iter=0;acc=10^(-5);
    repeat
    {
      iter=iter+1;

      for(i in 1:numc){
        f[,i] = dnorm(y,X%*%est_beta[,i],est_var[i]^0.5*rep(1,n));
      }

      for(i in 1:numc){
        r[,i] = est_p_z[,i]*f[,i]/apply(f*est_p_z,1,sum);
        s = sqrt(r[,i]);
        yy = y*s;xx = X*matrix(rep(s,2),nrow=n);
        b = coefficients(lm(yy~xx));
        #est_beta[,i] = c(b[2],b[3]);
        est_beta[,i] = b[-1];
        est_var[i] = sum(r[,i]*(y-X%*%est_beta[,i])^2)/sum(r[,i]);

        if(is.null(kernel)){
          W = exp(-(((matrix(rep(t(z),m),nrow=m)-matrix(rep(u,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));
          #W = 0.75*(1-(matrix(rep(t(z),m),nrow=m)-matrix(rep(u,n),ncol=n)))/h[i];W[W<0] = 0;
          R = t(matrix(rep(r[,i],m),ncol=m));
          est_p_u[,i] = apply(W*R,1,sum)/apply(W,1,sum);
        }else{
          est_p_u[,i] = local_constant(r[,i],z,u,h,kernel);
        }

        est_p_z[,i] = approx(u,est_p_u[,i],z)$y;
      }

      rsd = likelihood-sum(log(apply(f*est_p_z,1,sum)));
      likelihood = sum(log(apply(f*est_p_z,1,sum)));
      if(abs(rsd)<acc|iter>200){break}
    }

    if(res$lh<likelihood){
      out=list(est_p_u=res$est_p_u,est_p_z=res$est_p_z,est_beta=res$est_beta,est_var=res$est_var,lh=res$lh,initialtrue=0)
    }else{
      out=list(est_p_u=est_p_u,est_p_z=est_p_z,est_beta=est_beta,est_var=est_var,lh=likelihood,initialtrue=1)
    }
  }

  return(out)
}





gen_mixreg1_boot<-function(X,beta,var,p,numc){
  n=dim(X)[1];
  d=dim(X)[2];
  cum_p=c(p[,1],rep(1,n));
  cum_p=matrix(cum_p,ncol=2);
  gen_p=runif(n,0,1);
  bernulli=matrix(rep(0,n*numc),nrow=n);
  bernulli[,1]=(gen_p<cum_p[,1]);
  for(i in 2:numc){
    bernulli[,i]=(cum_p[,i-1]<gen_p)*(gen_p<cum_p[,i]);
  }
  e1=sqrt(var[1]);
  e2=sqrt(var[2]);
  func1=X%*%beta[,1]+e1*rnorm(n,0,1);
  func2=X%*%beta[,2]+e2*rnorm(n,0,1);
  func=matrix(c(func1,func2),ncol=2);
  y_boot=apply(func*bernulli,1,sum);
  jieguo_boot=list(y=y_boot);
  jieguo_boot
}



emest_linear<-function(x,y,numc){
  n=length(y)
  X = cbind(rep(1,n),x);
  d = dim(X)[2];

  est_p = (1/numc)*matrix(rep(1,numc),nrow=1);
  est_var = 0.1*matrix(rep(1,numc),nrow=1);
  est_beta = matrix(rnorm(d*numc,0,1),ncol=numc);

  n = length(y);
  f = matrix(numeric(n*numc),nrow=n);r=f;

  rsd=1;likelihood=1;iter=0;acc=10^(-3);

  repeat
  {
    iter = iter+1;

    for(i in 1:numc){
      f[,i] = dnorm(y,X%*%est_beta[,i],est_var[i]^0.5*rep(1,n));
    }

    ######
    for(i in 1:numc){
      r[,i] = est_p[i]*f[,i]/(f%*%t(est_p));#(f%*%est_p);apply(f*est_p,1,sum);
      s = sqrt(r[,i]);
      yy = y*s;xx = X*matrix(rep(s,2),nrow=n);
      b = coefficients(lm(yy~xx));
      #est_beta[,i] = c(b[2],b[3]);
      est_beta[,i] = b[-1];
      est_var[i] = sum(r[,i]*(y-X%*%est_beta[,i])^2)/sum(r[,i]);
    }
    ######
    est_p = cbind(mean(r[,1]),mean(r[,2]));

    rsd = likelihood-sum(log(f%*%t(est_p)));
    likelihood = sum(log(f%*%t(est_p)));
    if(abs(rsd)<acc|iter>20){break}
  }
  out=list(est_p=est_p,est_beta=est_beta,est_var=est_var)
  out
}



#y is the response; x is the preditor; u is the grid points where the unknow function is evaluated; h is the bandwidth
#mu is the output. we then run: plot(u,mu,'-')
#kernel, 1 means gaussian, 2 means epanechnikov
local_constant<-function(y,x,u,h,kernel){
  len1 = length(x);
  len2 = length(u);

  U1 = matrix(rep(u,len1),ncol=len1);
  U2 = matrix(rep(x,len2),ncol=len2);
  U = U1-t(U2);
  t = U/h;

  Y = t(matrix(rep(y,len2),ncol=len2));

  if(kernel == 1){
    W = h^-1*(2*pi)^-0.5*exp(-0.5*t^2);
  }
  if(kernel ==2){
    W = h^-1*0.75*(1-t^2);
    W[W<0] = 0;
  }

  p_u = apply(W*Y,1,sum)/apply(W,1,sum);
  p_u
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
