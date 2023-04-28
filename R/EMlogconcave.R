#' Parameter Estimation of Multivariate Normal Mixture Using EM Algorithm
#'
#'
#'Computes the parameter estimation of a possibly multivariate normal mixture model
#'using EM algorithm. The result is a input for EMlogconc function.
#'
#' @param x Row vector of observations.
#' @param k Gaussian mixture with k components
#'
#' @return List of parameters: logL:best log-likelihood; pis: estimated component proportions;
#' mu:estimated component means; sig: estimated component sd.
#'
#' @details
#' Uses random initialization for 20 restarts;
#' Stopping criterion: relative change in log likelihood < 10^(-8);
#'
#' @export
#'
#' @examples
#'
#'m=2;
#'x=rnorm(100,2,sqrt(2));x=matrix(x,nrow=100);x[1:60]=x[1:60]+5;
#'EMnorm = EMnormal(x,m);
#'
EMnormal <- function(x,k){
  # computes the parameter estimates of a (possibly multivariate)
  # Gaussian mixture with k components using EM.
  # The observations are the row vectors of x.
  # Uses random initializations for 20 restarts.
  # Stopping criterion: relative change in log likelihood < 10^(-8)

  x=as.matrix(x);
  n=dim(x)[1];
  d=dim(x)[2];

  tau=matrix(numeric(n*k),nrow=n);  #membership weight matrix
  lik=tau;
  bestlogL=-Inf; logL=1;

  for(restarts in 1:20){
    #Initialize random starting parameters:
    mu=rmvnorm(k,apply(x,2,mean),cov(x));  # mu is kxd
    sig=array(rep(1,d*d*k),dim=c(d,d,k));
    for(m in 1:k){
      sig[,,m]=cov(x);
    }
    pis=matrix(rep(1/k,k),nrow=1);

    for(iters in 1:1000){  #max. number of iterations
      #E-step:
      for(m in 1:k){
        tau[,m]=dmvnorm(x,mu[m,],as.matrix(sig[,,m]));
      }
      tau=(matrix(rep(1,n),nrow=n)%*%pis)*tau;
      tau=tau/(tau%*%matrix(rep(1,k*k),nrow=k));

      #M-step:
      for(m in 1:k){
        w=tau[,m]/sum(tau[,m]);
        mu[m,]=t(w)%*%x;
        xcenter=x-matrix(rep(1,n),nrow=n)%*%mu[m,];
        sig[,,m]=t(xcenter)%*%(xcenter*(w%*%matrix(rep(1,d),nrow=1)));
        lik[,m]=dmvnorm(x,mu[m,],as.matrix(sig[,,m]));
      }

      pis=matrix(apply(tau,2,sum)/n,nrow=1);
      oldlogL=logL;
      logL=sum(log(lik%*%t(pis)));

      if(iters>1 && (logL-oldlogL)/abs(oldlogL)<10^(-8)){
        break
      }
    }

    if(logL>bestlogL){
      bestlogL=logL;
      bestpis=pis;
      bestmu=mu;
      bestsig=sig;
    }
  }

  parameters=list(logL=bestlogL,pis=bestpis,mu=bestmu,sig=bestsig)
  parameters
}






#' Clustering with Mixtures of Log-concave Distributions using EM Algorithm
#'(Univariate)
#'
#'Does EM by computing the logconcave MLE for each component in the M-step,
#' with initial values from estimations from normal mixture assumption.
#' This function will need the output from EMnorm function.
#'
#' @param x Observations. Given as col vectors of x.
#' @param k the number of components.
#' @param EMnorm Output  of EMnorm function
#'
#' @return List of logL:log-likelihood; pis: estimated component proportions;
#' f:component densities at x.
#' @export
#'@importFrom mvtnorm rmvnorm dmvnorm
#' @examples
#'m=2;
#' set.seed(4)
#'x=rnorm(100,2,sqrt(2));x=matrix(x,nrow=100);x[1:60]=x[1:60]+5;
#'EMnorm = EMnormal(x,m);
#'EMlogc = EMlogconc(x,m,EMnorm);
#'
EMlogconc <- function(x,k,EMnorm){


  x=as.matrix(x);
  n=dim(x)[1];
  tau=matrix(numeric(n*k),nrow=n);  #membership probabilities
  f=tau;  #component densities at the x's

  #Initialize using the membership probs. from EMnormal:
  for(m in 1:k){
    f[,m]=dnorm(x,EMnorm$mu[m],sqrt(EMnorm$sig[m]));
  }
  pis=EMnorm$pis;

  for(iters in 1:5){
    #E-step:
    for(i in 1:n){
      tau[i,]=pis*f[i,]/as.numeric(f[i,]%*%t(pis));
    }

    #M-step:
    for(m in 1:k){
      f[,m]=exp(lcmle(x,tau[,m]));
    }

    pis=matrix(apply(tau,2,sum)/n,nrow=1);
    logL=sum(log(f%*%t(pis)));
  }

  parameters=list(logL=logL,pis=pis,f=f)
  parameters
}






lcmle <- function(x,w){
  # computes logMLE of a log-concave density using ICMA for f=exp(log) piecewise
  # linear, based on data x with weights w.
  # Merges data that are tied or too close (< std(x)/1000).
  # ICMA works on cone of increasing arguments, so work with s(n:-1:1) and
  # undo in calls to loglikelihd, grad loglikelihd. Return grad(n:-1:1) in latter

  # sort data and take care of ties and data that are too close:

  x=as.matrix(x);
  n=max(dim(x));
  x=t(x); w=t(w);
  w=w/sum(w);
  mu=sum(w*x);
  sig=sqrt(sum(w*(x-mu)^2));
  eps=sig/1000;  #eps/1000 is merging threshold
  xsort=sort(x,index.return = T)$x;idx=sort(x,index.return = T)$ix;
  wsort=w[idx]; xx=xsort[1]; ww=wsort[1]; nn=1;

  for(i in 2:n){
    if(xsort[i]-xx[nn]<eps){
      ww[nn]=ww[nn]+wsort[i];
    }else{
      xx=cbind(xx,xsort[i]);
      ww=cbind(ww,wsort[i]);
      nn=nn+1;
    }
  }

  ww=ww/sum(ww);

  #fit normal to initialize:
  xx1=xx[1]-sig/10;  #xx1 is dummy variable
  xx=cbind(xx1,xx);
  d=min(-10,-((xx1-mu)/sig)^2);  #log f at dummy xx(1)=xx1
  ww=cbind(1,ww);  #ww(1) is dummy weight
  cw=cumsum(ww[seq(nn+1,2,-1)]);
  cw=cw[seq(nn,1,-1)]*diff(as.numeric(xx));
  swork=1:nn;  #construct starting values for slopes:
  swork[1]=(log(dnorm(xx[2],mu,sig))-d)/(xx[2]-xx[1]);
  for(i in 2:nn){
    swork[i]=(log(dnorm(xx[i+1],mu,sig)/dnorm(xx[i],mu,sig)))/(xx[i+1]-xx[i]);
  }

  #start ICMA:
  repeatICMA=1;
  while(repeatICMA>0){  #is dummy d = log f(xx1) not small enough?
    d=4*d;
    swork[1]=swork[1]-3/4*d/(xx[2]-xx[1]);  #adjust d and swork(1)
    swork=swork[seq(nn,1,-1)];  #working value is isotonic
    swork=ICMA(xx,d,swork,cw);  #compute MLE
    swork=swork[seq(nn,1,-1)];  # working value swork is isotonic

    if(swork[1]>1.1*swork[2]){
      repeatICMA=0;  #d is small enough
    }
  }

  #evaluate logMLE at xx:
  logf=matrix(rep(0,n),nrow=1); logf[1]=d;
  for(i in 2:(nn+1)){
    logf[i]=logf[i-1]+swork[i-1]*(xx[i]-xx[i-1]);
  }
  logf=logf[2:(nn+1)];  #ignore dummy variable xx1
  xx=xx[2:(nn+1)];

  #evaluate logMLE at xsort:
  logfsort=matrix(rep(0,n),nrow=1); logfsort[1]=logf[1]; xxidx=1;
  for(i in 2:n){
    if(xsort[i]-xx[xxidx]<eps){
      logfsort[i]=logf[xxidx];
    }else{
      xxidx=xxidx+1;
      logfsort[i]=logf[xxidx];
    }
  }

  #evaluate logMLE at x:
  logf=matrix(rep(0,n),nrow=1); logf[idx]=logfsort;  #undo the sorting
  logf
}




ICMA <- function(x,d,s,cw){
  # modified ICMA using Hessian weights, and the following further modifications:
  # also terminates if the gain in phi is < eta, or if more than 500 iterations;
  # bounds weights below at 10^-3 to avoid problems when using w^(-1);
  # performs at most 6 bisections for the line search

  eta <- 0.00001  # accuracy parameter
  eps <- 0.01   # line search parameter
  x=as.matrix(x);
  n=max(dim(x))-1;
  #s: starting value
  grad=gradphi(x,d,s,cw);
  w=mmax(0.001,grgrphi(x,d,s));    # second derivatives as weights
  t1=abs(sum(s*grad));
  t2=abs(sum(grad));
  t3=min(cumsum(grad[seq(n,1,-1)]));
  ctr=0;  gain=1;

  while((t1>eta || t2>eta || t3<eta) && ctr<500 && gain>eta){
    oldphi=phi(x,d,s,cw);
    ytil=MINLOWSET(c(s-(w^(-1))*grad,w));

    if(phi(x,d,ytil,cw)<phi(x,d,s,cw)+eps*sum(grad*(ytil-s))){
      s=ytil;
    }else{
      lam=1; ss=0.5; z=ytil; ctr2=0;
      t4=(phi(x,d,z,cw)<phi(x,d,s,cw)+(1-eps)*sum(grad*(z-s)));
      t5=(phi(x,d,z,cw)>phi(x,d,s,cw)+eps*sum(grad*(z-s)));
      while((t4>0.5 || t5>0.5) && ctr2<6){
        if(t4>0.5){
          lam=lam+ss;
        }else{
          lam=lam-ss;
        }

        z=s+lam*(ytil-s);
        ss=ss/2;
        t4=(phi(x,d,z,cw)<phi(x,d,s,cw)+(1-eps)*sum(grad*(z-s)));
        t5=(phi(x,d,z,cw)>phi(x,d,s,cw)+eps*sum(grad*(z-s)));
        ctr2=ctr2+1;
      }

      s=z;
    }

    grad=gradphi(x,d,s,cw);
    w=mmax(0.001,grgrphi(x,d,s));
    t1=abs(sum(s*grad));
    t2=abs(sum(grad));
    t3=min(cumsum(grad[seq(n,1,-1)]));
    ctr=ctr+1;
    gain=oldphi-phi(x,d,s,cw);
  }

  y=s;
  if(ctr>499){
    warning("Number of iterations in ICMA exceeded")
  }

  y
}



gradphi <- function(x,d,s,cw){
  # computes gradient of phi at s
  # NOTE: switches argument s to s(n:-1:1), consequently switches result(n:-1:1)

  x=as.matrix(x);
  n=max(dim(x))-1;
  s=s[seq(n,1,-1)];
  a=matrix(rep(0,n+1),nrow=1);
  y=matrix(rep(0,n),nrow=1);
  x1sq=x[1]^2;

  for(i in 2:(n+1)){
    a[i]=a[i-1]+s[i-1]*(x[i]-x[i-1]);
  }
  t=(x[n+1]-x[n])*exp(d+a[n+1]);
  y[n]=0.5*t*(x[n+1]-x[n]);

  for(k in seq(n-1,2,-1)){
    t=t+(x[k+2]-x[k])*exp(d+a[k+1]);
    y[k]=0.5*t*(x[k+1]-x[k]);
  }

  t=(x[3]-x[2])*exp(d+a[2]);
  y[1]=0.5*t*(x[2]-x[1]);
  y=y-cw;
  y=y[seq(n,1,-1)];

  y
}




grgrphi <- function(x,d,s){
  # computes diagonal of second derivative of phi at s
  # NOTE: switches argument s to s(n:-1:1), consequently switches result(n:-1:1)

  x=as.matrix(x);
  n=max(dim(x))-1;
  s=s[seq(n,1,-1)];
  a=matrix(rep(0,n+1),nrow=1);
  y=matrix(rep(0,n),nrow=1);
  x1sq=x[1]*x[1];

  for(i in 2:(n+1)){
    a[i]=a[i-1]+s[i-1]*(x[i]-x[i-1]);
  }
  t=(x[n+1]-x[n])*exp(d+a[n+1]);
  y[n]=0.5*t*(x[n+1]-x[n])^2;

  for(k in seq(n-1,2,-1)){
    t=t+(x[k+2]-x[k])*exp(d+a[k+1]);
    y[k]=0.5*t*(x[k+1]-x[k])^2;
  }

  t=(x[3]-x[2])*exp(d+a[2]);
  y[1]=0.5*t*(x[2]-x[1])^2;
  y=y[seq(n,1,-1)];

  y
}




phi <- function(x,d,s,cw){
  # computes -(likelihood - integrated likelihood) for log(concave +ct^2,
  # linearized), for f=exp(log) is piecewise linear
  # NOTE: switches argument s to s(n:-1:1)

  x=as.matrix(x);
  n=max(dim(x))-1;
  s=s[seq(n,1,-1)];
  a=matrix(rep(0,n+1),nrow=1);
  x1sq=x[1]^2;
  a[2]=s[1]*(x[2]-x[1]);
  y=(x[3]-x[2])*exp(d+a[2]);

  for(i in 3:n){
    a[i]=a[i-1]+s[i-1]*(x[i]-x[i-1]);
    y=y+(x[i+1]-x[i-1])*exp(d+a[i]);
  }

  a[n+1]=a[n]+s[n]*(x[n+1]-x[n]);
  y=y+(x[n+1]-x[n])*exp(d+a[n+1]);
  y=0.5*y-sum(s*cw);

  y
}





MINLOWSET <- function(gw){
  # minimum lower sets algorithm, see RWD, p. 24
  # computes isotonic regression of g w.r.t. weights w, i.e. slope of
  # greatest convex minorant (i.e. alternative to PAVA)

  n=length(gw)/2;
  g=gw[1:n];
  w=gw[(n+1):(2*n)];
  curr=1;

  while(curr<n+0.5){
    h=cumsum(g[curr:n]*w[curr:n])/cumsum(w[curr:n]);
    a=min(h[seq(n-curr+1,1,-1)]);
    ind=which.min(h[seq(n-curr+1,1,-1)]); ind=n+1-ind;
    g[curr:ind]=seq(a,a,ind-curr+1);
    curr=ind+1;
  }
  y=g;

  y
}






#' inverse fcdf
#'
#' @param x observations
#'
#' @return pnorm
#'
#'@importFrom GoFKernel inverse
#'
fcdf=function(x){
  y=pnorm(x)
  y
}
finv <- GoFKernel::inverse(fcdf)

#' Clustering with mixtures of log-concave distributions using EM Algorithm
#' (Multivariate)
#'
#' Does EM by computing the logconcave MLE for each marginal in
#' each component in the M-step. Then use the normal copula to
#' model the correlation structure in each component:
#' Transform each marginal to a normal using the estimated cdf.
#'
#' @param x Observations. Given as row vectors of x.
#' @param k the number of components.
#' @param EMnorm Output of EMnorm function
#' @param mixiter maximum iteration the algorithm run for EM-algorithm, defualt is 5.
#'
#' @return List of logL:log-likelihood; pis: estimated component proportions;
#' f:component densities at the observations;Sig: Estimated covariances of the transformed marginals;
#' fh:logconcave marginal component densities at x; Fh:cdfs of the fh;
#'@importFrom mvtnorm rmvnorm dmvnorm
#'@export
#'
#' @examples
#'
#' m=2;
#' x=mvtnorm::rmvnorm(100,c(0,0),matrix(c(2,1,1,2),nrow=2));x=matrix(x,nrow=100);
#' x[1:60,]=x[1:60,]+5;
#' EMnorm = EMnormal(x,m);
#' \donttest{EMlogc = EMlogconcHD(x,m,EMnorm,mixiter=1);}
#'
#'
EMlogconcHD <- function(x,k,EMnorm,mixiter=5){
  # Does EM by computing the logconcave MLE for each marginal in
  # each component in the M-step. Then use the normal copula to
  # model the correlation structure in each component:
  # Transform each marginal to a normal using the estimated cdf.
  # Initializes using the output EMnorm of EM with normal components.
  # Observations are given as row vectors of x. k is the number of components.

  x=as.matrix(x);
  n=dim(x)[1];
  d=dim(x)[2];

  tau=matrix(numeric(n*k),nrow=n);  #membership weight matrix
  f=tau;  #component densities at the x's
  fh=array(rep(1,n*d*k),dim=c(n,d,k));  #logconcave marginal component densities at the x's
  Fh=fh;  #cdfs of the fh
  Sig=array(rep(1,d*d*k),dim=c(d,d,k));  #covariances of the transformed marginals

  #Initialize using the membership probs. from EMnormal:

  for(m in 1:k){
    f[,m]=dmvnorm(x,EMnorm$mu[m,],as.matrix(EMnorm$sig[,,m]));
  }
  pis=EMnorm$pis;

  for(iters in 1:mixiter){
    #E-step:
    tau=(matrix(rep(1,n),nrow=n)%*%pis)*f;
    tau=tau/(tau%*%matrix(rep(1,k*k),nrow=k));

    #M-step:
    for(m in 1:k){
      J=matrix(rep(1,n),nrow=n);  #initialize Jacobian
      Y=matrix(rep(0,n*d),nrow=n);  #Y is transformed data

      for(j  in 1:d){
        fh[,j,m]=exp(lcmle(x[,j],tau[,m]));
        xs=sort(x[,j],index.return = T)$x;idx=sort(x[,j],index.return = T)$ix;
        fhs=fh[idx,j,m];  #xs is sorted data, fhs corresponds to x
        cdfincrements=(xs[2:n]-xs[1:(n-1)])*(fhs[2:n]+fhs[1:(n-1)])/2;
        cdf=cumsum(c(0,cdfincrements));
        #scale so that cdf(x1)=1/(n+1), cdf(xn)=n/(n+1):
        Fh[idx,j,m]=1/(n+1)+cdf/cdf[n]*(n-1)/(n+1);
        for(l in 1:n){
          Y[l,j]=finv(Fh[l,j,m]);
        }
        J=J*fh[,j,m]/dnorm(Y[,j]);
      }

      w=tau[,m]/sum(tau[,m]);
      ycenter=Y-matrix(rep(1,n),nrow=n)%*%(t(w)%*%Y);
      Sig[,,m]=t(ycenter)%*%(ycenter*(w%*%matrix(rep(1,d),nrow=1)));
      f[,m]=dmvnorm(Y,sigma=as.matrix(Sig[,,m]))*J;
    }

    pis=matrix(apply(tau,2,sum)/n,nrow=1);
    logL=sum(log(f%*%t(pis)));
  }

  parameters=list(logL=logL,pis=pis,f=f,Sig=Sig,fh=fh,Fh=Fh)
  parameters
}







