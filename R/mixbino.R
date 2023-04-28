
#' Nonparametric Mixture of Binomial Regression with a Degenerate Component
#'
#'Fits the nonparametric mixture of binomial distribution with one degenerate component
#' \cr \eqn{w(t)*B(N,p(t))+(1-w(t))*B(N,0)}
#' @param tg the grid points on which we want to evaluate the function w(t) and p(t).
#' @param t the time variable associated with x.
#' @param x the integer observation.
#' @param N the number of experiments for Binomial distribution.
#' @param tune the bandwidth will be h*tune. Default will be 1.
#' @param acc1 stopping criteria
#'
#' @return
#' w contains the estimate of first component proportion function
#' pt is the estimate of component probability function for the first component
#' @export
#'
#' @examples
#' n=100;tg=seq(from=0,to=1,length.out=50);t=seq(from=0,to=1,length.out=n);
#' pt=0.5*(1-cos(2*pi*t));
#' b=rbinom(n,1,0.2);
#' y=numeric();
#' for(i in 1:n){
#'  if(b[i]==1){
#'   y[i]=0;
#'   }else{
#'       y[i]=rbinom(1,7,pt[i]);
#'         }
#'         }
#'  ft=mixbino(tg,t,y,7);
#'  mean(ft$w)
#'  p=0.5*(1-cos(2*pi*tg));
mixbino<-function(tg,t,x,N,tune = NULL,acc1 = NULL){
  if(length(x)!=length(t)){
    errorCondition("x and t are in different sizes")
  }
  if(is.null(tune)){
    tune = 1;
  }

  #make x and t row vector
  x = as.matrix(x); t = as.matrix(t);
  if(dim(x)[1]==1){
    x = x; t = t;
  }else{
    x = t(x); t=t(t);
  }

  #Create initial values for w and pt, where pt records the estimate of p(t)
  #evaluated at all observed t values.
  #p recrods the classification probability that each observation from the first component.
  lx = length(x);
  temp = 1:lx;
  ind1 = temp[x==0]; ind2 = temp[x!=0];
  w = length(ind1)/lx;
  pt = mean(x[ind2])/N;
  p = numeric();
  p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt^x[ind1]*(1-pt)^(N-x[ind1]));
  p[ind2] = 0;
  ltg = length(tg);
  acc = 0.01;

  if(tune>0.2){#the final bandwidth will be h*tune
    subx = x[ind2]/N;
    subt = t[ind2];#subx and subt must from the second component
    lsx = length(subx);
    hx = median(abs(subx-median(subx)))/0.6745*(4/3/lsx)^0.2;
    hy = median(abs(subt-median(subt)))/0.6745*(4/3/lsx)^0.2;
    h = sqrt(hy*hx);
    if(is.null(acc1)){
      acc1 = min(0.01,0.01*h);
    }
    difh = acc1+1;

    while(difh>acc1){
      preh = h; PT = numeric(); W = numeric();
      for(i in 1:ltg){
      #for(i in 1:50){
        dif = 1;
        while(dif>acc){
          prew = w;
          #M step to find pt and w
          temp = dnorm(t,tg[i],h)*(1-p);
          pt = sum(temp*x)/(sum(temp))/N;
          temp1 = dnorm(t,tg[i],h)%*%p;
          w = sum(temp1)/(sum(temp1)+sum(temp));
          dif = abs(w-prew); if(length(dif)>1){dif=min(dif)}
          p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt^x[ind1]*(1-pt)^(N-x[ind1]));
        }
        PT[i] = pt; W[i] = 1-w;
      }
      out=list(pt=PT,w=W)

      w = approx(tg,out$w,t)$y; #问题在这
      pt = approx(tg,out$pt,t)$y;
      p[ind1] = w[ind1]/(w[ind1]+(1-w[ind1])*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt[ind1]^x[ind1]*(1-pt[ind1])^(N-x[ind1]));
      lsx = lx-sum(p);
      mx = (1-p)%*%t(x)/lsx; mx = as.numeric(mx);
      mt = (1-p)%*%t(t)/lsx; mt = as.numeric(mt);
      hx = sqrt((1-p)%*%t(x-mx)^2/lsx)/N*(4/3/lsx)^0.2;
      hy = sqrt((1-p)%*%t(t-mt)^2/lsx)/N*(4/3/lsx)^0.2;
      h = sqrt(hy*hx); difh = abs(h-preh); difh = as.numeric(difh);
    }
  }else{#tune means the percentage of data included for local estimation
    h = tune*(range(t)[2]-range(t)[1]);
    for(i in 1:ltg){
      dif = 1;
      while(dif>acc){
        prew = w;
        temp = dnorm(t,tg[i],h)*(1-p);
        pt = sum(temp*x)/(sum(temp))/N;
        temp1 = dnorm(t,tg[i],h)%*%p;
        w = sum(temp1)/(sum(temp1)+sum(temp));
        dif = abs(w-prew);
        p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt^x[ind1]*(1-pt)^(N-x[ind1]));
      }
      PT[i] = pt; W[i] = 1-w;
      out=list(pt=PT,w=W)
    }
  }

  out$h=h;
  return(out)
}



#mixbinosemi(tg,t,x,N)
#fits the semiparametric mixture of binomial distribution with one degenerate component
#w*B(N,0)+(1-w)*B(N,p(t))
#INPUT
#ditto
#h: the bandwidth
#OUTPUT
#w contains the estimate of component proportion
#pt is the estimate of component probability function



#' Semiparametric Mixture of Binomial Regression with a Degenerate Component
#'
#'Fits the semiparametric mixture of binomial distribution with one degenerate component
#'\cr \eqn{w*B(N,0)+(1-w)*B(N,p(t))}
#'
#' @param tg the grid points on which we want to evaluate the function w(t) and p(t).
#' @param t the time variable associated with x.
#' @param x the integer observation.
#' @param N the number of experiments for Binomial distribution.
#' @param h the bandwidth. If NULL, will be estimated by local constant regression suggested by Bowman and Azzalini (1997) p.31.
#'
#' @return
#' w contains the estimate of first component proportion function
#' pt is the estimate of component probability function for the first component
#' h the bandwidth
#' @export
#'
#' @examples
#' n=100;tg=seq(from=0,to=1,length.out=50);t=seq(from=0,to=1,length.out=n);
#' pt=0.5*(1-cos(2*pi*t));
#' b=rbinom(n,1,0.2);
#' y=numeric();
#' for(i in 1:n){
#'  if(b[i]==1){
#'   y[i]=0;
#'   }else{
#'       y[i]=rbinom(1,7,pt[i]);
#'         }
#'         }
#'  ft=mixbinosemi(tg,t,y,7);
mixbinosemi<-function(tg,t,x,N,h = NULL){
  if(length(x)!=length(t)){
    errorCondition("x and t are in different sizes")
  }

  x = as.matrix(x); t = as.matrix(t);
  if(dim(x)[1]==1){
    x = x; t = t;
  }else{
    x = t(x); t=t(t);
  }

  lx = length(x);
  temp = 1:lx;
  ind1 = temp[x==0]; ind2 = temp[x!=0];
  w = length(ind1)/lx;
  pt = mean(x[ind2])/N;
  p = numeric();
  p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt^x[ind1]*(1-pt)^(N-x[ind1]));
  p[ind2] = 0;
  acc = 0.001;

  #Find the optimal bandwidth h for estimating the pt
  #optimal bandwidth for local constant regression suggested by Bowman and Azzalini (1997) p.31
  if(is.null(h)){
    subx = x[ind2]/N;
    subt = t[ind2];#subx and subt must from the second component
    lsx = length(subx);
    hx = median(abs(subx-median(subx)))/0.6745*(4/3/lsx)^0.2;
    hy = median(abs(subt-median(subt)))/0.6745*(4/3/lsx)^0.2;
    h = sqrt(hy*hx);
    acc1 = min(0.01,0.01*h);
    difh = 1;

    while(difh>acc1){
      preh = h; dif = 1;
      while(dif>acc){
        prew = w;
        #M step to find pt and w
        pt = numeric();
        for(i in 1:lx){
          temp = dnorm(t,t[i],h)*(1-p);
          pt[i] = sum(temp*x)/(sum(temp))/N;
        }
        w = sum(p)/lx;
        dif = abs(w-prew);
        p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt[ind1]^x[ind1]*(1-pt[ind1])^(N-x[ind1]));
      }
      lsx = lx-sum(p);
      mx = (1-p)%*%t(x)/lsx; mx = as.numeric(mx);
      mt = (1-p)%*%t(t)/lsx; mt = as.numeric(mt);
      hx = sqrt((1-p)%*%t(x-mx)^2/lsx)/N*(4/3/lsx)^0.2;
      hy = sqrt((1-p)%*%t(t-mt)^2/lsx)/N*(4/3/lsx)^0.2;
      h = sqrt(hy*hx); difh = abs(h-preh);
    }
  }else{
    dif = 1;
    while(dif>acc){
      prew = w;
      #M step to find pt and w
      pt = numeric();
      for(i in 1:lx){
        temp = dnorm(t,t[i],h)*(1-p);
        pt[i] = sum(temp*x)/(sum(temp))/N;
      }
      w = sum(p)/lx;
      dif = abs(w-prew);
      p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt[ind1]^x[ind1]*(1-pt[ind1])^(N-x[ind1]));
    }
  }

  pt = approx(t,pt,tg)$y;
  out=list(pt=pt,w=w,h=h)
  return(out)
}



#' One-Step Estimation of Semiparametric Mixture of Binomial Regression with a Degenerate Component
#'
#'Fits the semiparametric mixture of binomial distribution with one degenerate component using one-step
#'backfitting procedure
#'
#' @param tg the grid points on which we want to evaluate the function w(t) and p(t).
#' @param t the time variable associated with x.
#' @param x the integer observation.
#' @param N the number of experiments for Binomial distribution.
#' @param tune the bandwidth will be h*tune. Default will be 1.
#'
#' @return
#' w contains the estimate of first component proportion function
#' pt is the estimate of component probability function for the first component
#' h the bandwidth
#' @export
#'
#' @examples
#'
#' nobs=50;
#' tobs=seq(from=0,to=1,length.out=nobs);
#' pi1Tru=0.4*matrix(rep(1,nobs),nrow=nobs);
#' ptTru=0.3*(1.5+cos(2*pi*tobs));
#' nfine=nobs;
#' tfine=seq(from=0,to=1,length.out=nfine);
#' yobs=numeric();
#' for(i in 1:nobs){
#'   b=rbinom(1,1,pi1Tru[i]);
#'   if(b==1){
#'       yobs[i]=0;
#'    }else{
#'       yobs[i]=rbinom(1,7,ptTru[i]);
#'       }
#'  }
#'  ftCon=mixbinosemi(tfine,tobs,yobs,7);
#'  pi1Con=ftCon$w;
#'  ptCon=ftCon$pt;
#'  ftonestep=mixbinosemionestep(tfine,tobs,yobs,7);
#'  pi1onestep=ftonestep$w;
#'  ptonestepMat=ftonestep$pt;
mixbinosemionestep<-function(tg,t,x,N,tune = NULL){
  if(length(x)!=length(t)){
    errorCondition("x and t are in different sizes")
  }
  if(is.null(tune)){
    tune = 1;
  }

  #make x and t row vector
  x = as.matrix(x); t = as.matrix(t);
  if(dim(x)[1]==1){
    x = x; t = t;
  }else{
    x = t(x); t=t(t);
  }

  lx = length(x);
  temp = 1:lx;
  ind1 = temp[x==0]; ind2 = temp[x!=0];
  out = mixbino(tg,t,x,N,tune,1);
  pt = out$pt;
  w = mean(out$w);
  p = numeric();
  p[ind2] = 0;
  acc = 0.001; dif = 1;

  while(dif>acc){
    prew = w;
    #M step to find pt and w
    p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt[ind1]^x[ind1]*(1-pt[ind1])^(N-x[ind1]));
    w = sum(p)/lx;
    dif = abs(w-prew);
  }
  lsx = lx-sum(p);
  mx = (1-p)%*%t(x)/lsx; mx = as.numeric(mx);
  mt = (1-p)%*%t(t)/lsx; mt = as.numeric(mt);
  hx = sqrt((1-p)%*%t(x-mx)^2/lsx)/N*(4/3/lsx)^0.2;
  hy = sqrt((1-p)%*%t(t-mt)^2/lsx)/N*(4/3/lsx)^0.2;
  h = sqrt(hy*hx); run = 0; dif = 1;

  while(dif>acc  && run<1000){
    prep = pt; run = run+1;
    #M step to find pt and w
    pt = numeric();
    for(i in 1:lx){
      temp = dnorm(t,t[i],h)*(1-p);
      pt[i] = sum(temp*x)/(sum(temp))/N;
    }
    dif = max(abs(pt-prep));
    p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt[ind1]^x[ind1]*(1-pt[ind1])^(N-x[ind1]));
  }

  pt = approx(t,pt,tg)$y;
  out=list(pt=pt,w=w,h=h)
  return(out)
}



#'Estimation of Semiparametric Mixture of Binomial Regression with a Degenerate Component Using Full Iterative
#'Backfitting Procedure
#'
#'Fits the semiparametric mixture of binomial distribution with one degenerate component using full iterative
#'backfitting procedure
#' @param tg the grid points on which we want to evaluate the function w(t) and p(t).
#' @param t the time variable associated with x.
#' @param x the integer observation.
#' @param N the number of experiments for Binomial distribution.
#' @param tune the bandwidth will be h*tune. Default will be 1.
#'
#' @return
#' w contains the estimate of first component proportion function
#' pt is the estimate of component probability function for the first component
#' h the bandwidth
#' @export
#'
#' @examples
#'
#' nobs=50;
#' tobs=seq(from=0,to=1,length.out=nobs);
#' pi1Tru=0.4*matrix(rep(1,nobs),nrow=nobs);
#' ptTru=0.3*(1.5+cos(2*pi*tobs));
#' nfine=nobs;
#' tfine=seq(from=0,to=1,length.out=nfine);
#' yobs=numeric();
#' for(i in 1:nobs){
#'   b=rbinom(1,1,pi1Tru[i]);
#'   if(b==1){
#'       yobs[i]=0;
#'    }else{
#'       yobs[i]=rbinom(1,7,ptTru[i]);
#'       }
#'  }
#'  ftfull=mixbinosemifull(tfine,tobs,yobs,7);
#'  pi1full=ftfull$w;
#'  ptfull=ftfull$pt;
mixbinosemifull<-function(tg,t,x,N,tune = NULL){
  if(length(x)!=length(t)){
    errorCondition("x and t are in different sizes")
  }
  if(is.null(tune)){
    tune = 1;
  }

  #make x and t row vector
  x = as.matrix(x); t = as.matrix(t);
  if(dim(x)[1]==1){
    x = x; t = t;
  }else{
    x = t(x); t=t(t);
  }

  lx = length(x);
  temp = 1:lx;
  ind1 = temp[x==0]; ind2 = temp[x!=0];
  out = mixbino(tg,t,x,N,tune=tune);
  pt = out$pt;
  w = mean(out$w);
  p = numeric();
  p[ind2] = 0;
  acc = 0.001; dif1 = 1; numiter = 0;
  while(dif1>acc*0.01 && numiter<500){
    prew1 = w; dif = 1; numiter = numiter+1;
    while(dif>acc){
      prew = w;
      #M step to find pt and w
      p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt[ind1]^x[ind1]*(1-pt[ind1])^(N-x[ind1]));
      w = sum(p)/lx;
      dif = abs(w-prew);
    }
    lsx = lx-sum(p);
    mx = (1-p)%*%t(x)/lsx; mx = as.numeric(mx);
    mt = (1-p)%*%t(t)/lsx; mt = as.numeric(mt);
    hx = sqrt((1-p)%*%t(x-mx)^2/lsx)/N*(4/3/lsx)^0.2;
    hy = sqrt((1-p)%*%t(t-mt)^2/lsx)/N*(4/3/lsx)^0.2;
    h = sqrt(hy*hx); dif = 1; run = 0;
    while(dif>acc  && run<1000){
      prep = pt; run = run+1;
      #M step to find pt and w
      pt = numeric();
      for(i in 1:lx){
        temp = dnorm(t,t[i],h)*(1-p);
        pt[i] = sum(temp*x)/(sum(temp))/N;
      }
      dif = max(abs(pt-prep));
      p[ind1] = w/(w+(1-w)*gamma(N+1)/gamma(x[ind1]+1)/gamma(N+1-x[ind1])*pt[ind1]^x[ind1]*(1-pt[ind1])^(N-x[ind1]));
    }
    dif1 = abs(prew1-w);
  }

  pt = approx(t,pt,tg)$y;
  out=list(pt=pt,w=w,h=h)
  return(out)
}



