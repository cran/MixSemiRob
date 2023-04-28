




#' Minimum profile Hellinger distance estimation for a semiparametric mixture model
#'
#'A effective estimator for a class of semiparametric mixture models where one component has known distribution
#'with possibly unknown parameters while the other component density and the mixing proportion are unknown.
#'The proposed estimator is based on the minimum profile Hellinger distance (MPHD).
#'
#' @param x vector of observations
#' @param sigma sd of the known component; Can be set as NULL if the parameter is unknown.
#' @param ini initial values of parameters of components.  Format: list(mu,pi,sigma)
#' If ini is not specified, the function will derive initial values from x with equal
#' unknown variance normal mixture assumption.
#' @param true if the sigma is known: true = list(p,mu); if the sigma is unknown true =list(p,sigma,mu)
#'
#' @return list of estimations:  numiter: number of iteration
#' initialtrue: If the true value specified by the user fits the simulated data better return 1, 0 o.w.;
#' For below parameters, the returned value is based on initialtrue. If initialtrue = 1 return true parameters specified by the user,
#' if initialtrue=0 return estimated parameters:
#' sigma:estimated sigma when sigma is unknown;pi: estimated proportion; mu: estimated mean;
#' @export
#'
#'@importFrom stats nlminb
#'
#' @examples
#'
#' #Model: X-0.3N(0,1)+0.7N(3,1)
#' set.seed(4)
#' n=100;sigma2=1;p=0.3;mu=1.5;sigma=1;k=2;n1=rbinom(1,n,p);
#' x1=rnorm(sum(n1),0,sigma);x2=rnorm(n-sum(n1),mu,sigma2);x=c(x1,x2);
#'
#' true=c(p,mu);
#' temp=mixonekn(x,sigma);
#' mphdest=mphd(x,sigma,ini = temp,true = true)
#'
#'
#'
mphd<-function(x,sigma = NULL,ini = NULL,true = NULL){
  if(is.null(sigma)){#true = c(p,sigma,mu);
    out=mhdem1_x(x,ini,true)
  }else{#true = c(p,mu);
    out=mhde2m1_x(x,sigma,ini,true)
  }

  return(out)
}





#Reliable and extremely fast kernel density estimator for one-dimensional data
#Gaussian kernel is assumed and the bandwidth is chosen automatically
#dct1d;fixed_point
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

#using the interp1 to interpolate the response curse at xgrid.
interpcut<-function(x,y,xgrid){
  a = sort(x,index.return = T)$x; b = sort(x,index.return = T)$ix
  x = a; y = y[b];
  c = which(duplicated(x)==FALSE)
  x = x[c]; y = y[c];
  my = 0;#if xgrid is out of bound, we set them equal to 0.
  i = which(xgrid>max(x));
  j = which(xgrid<min(x));
  out = approx(x,y,xgrid)$y;
  out[i] = my; out[j] = my;
  out = mmax(0,out)
}



#estimate the two component mixture model assuming the first component has mean 0.
#If sig is given, then the first component has standard deviation "sig"
#mixonecompkn

#' Two normal component mixture estimation
#'
#' Estimate the two component mixture model assuming the first component has mean 0.
#' This function is used in the example of mphd to give initial estimation.
#' @param x observations of two mixture normal
#' @param sig standard deviation of first component. By default sig is set to be NULL.
#'
#' @return mu: estimated mean; sigma: estimated sigma; pi:estimated proportion; lh: estimated likelihood for observations.
#' @export
#'
#' @examples
#' #see example in mphd
mixonekn<-function(x,sig = NULL){
  n = length(x);
  s2 = var(x)*(n-1)/n;
  sx = sort(x);
  numini = 10;
  mugrid = seq(from=mean(x),to=sx[round(n*0.95)],length.out=10);

  l = numeric();
  mu = matrix(rep(0,2*numini),nrow=numini); sigma = mu; prop = mu;

  #%run mixonecompkn based on ten intial value
  if(is.null(sig)){
    for(i in 1:numini){
      muini = c(0,mugrid[1]);
      ind1 = which(x^2<(x-muini[2])^2);
      ind2 = which(x^2>(x-muini[2])^2);
      sigmaini = c(sd(x[ind1]),sd(x[ind2]));
      propini = c(length(ind1),length(ind2))/n;
      temp = mixonecompkn(x,muini,sigmaini,propini);
      l[i] = temp$lh; mu[i,] = temp$mu;
      sigma[i,] = temp$sigma; prop[i,] = temp$pi;
    }
  }else{
    for(i in 1:numini){
      muini = c(0,mugrid[1]);
      ind1 = which(x^2<(x-muini[2])^2);
      ind2 = which(x^2>(x-muini[2])^2);
      sigmaini = c(sig,sd(x[ind2]));
      propini = c(length(ind1),length(ind2))/n;
      temp = mixonecompkn(x,muini,sigmaini,propini,sig);
      l[i] = temp$lh; mu[i,] = temp$mu;
      sigma[i,] = temp$sigma; prop[i,] = temp$pi;
    }
  }

  #choose the converged value which has the largest likelihood
  #y = sort((-l),index.return = T)$x;
  id = sort((-l),index.return = T)$ix; id = id[1];
  mu = mu[id,]; prop = prop[id,];sigma = sigma[id,];
  #a = sort((mu),index.return = T)$x;
  b = sort((mu),index.return = T)$ix;

  out=list(mu=mu[b],sigma=sigma[b],pi=prop[b],lh=l)
  return(out)
}



#mmax
mixonecompkn<-function(x,mu,sigma,prop,sig = NULL){
  trimprop = 0.05;#the proportion of samples being trimmed.
  n = length(x); acc = 10^(-5)*n; oldx = x; n1 = ceiling(n*(1-trimprop));
  mu = c(0,3);sigma = c(1,1);prop = c(0.5,0.5);

  if(is.null(sig)){
    stop=0; run=0; dif=10; difml=10; ml=10^8;

    #use em algorithm to calculate the mle
    while(stop == 0){

      denm = t(cbind(as.matrix(dnorm(oldx,mu[1],sigma[1])),as.matrix(dnorm(oldx,mu[2],sigma[2]))));
      f = apply(prop*denm,2,sum);
      #a = sort((-f),index.return = T)$x;
      b = sort((-f),index.return = T)$ix;
      x = oldx[b[1:n1]]; x = matrix(x,nrow=1);
      run = run+1;
      denm = denm[,b[1:n1]];
      f = f[b[1:n1]];
      # lh[run] = sum(log(f));#use lh to record the sequence of lh
      lh = c(lh,sum(log(f))) #updated oct23 2022 by Xin Shen to make lh a vector instead of ts to fix a bug

      if(run > 2){
        denom = lh[run-1]-lh[run-2];
        if(denom == 0){
          stop = 1;
        }else{
          c = (lh[run]-lh[run-1])/denom;#convergence rate of Em algorithm
          if(c>0 && c<1){
            preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
            dif = ml-lh[run]; difml = abs(preml-ml);
          }
        }
      }

      #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
      if(dif<0.001 && difml<0.0001 && run>20 || run>200 || is.nan(min(mu)) == TRUE){
        stop = 1;
      }else{
        p = matrix(rep(prop,n1),ncol=n1)*denm/rbind(f,f);#E-step, classification probability of x from each component
        sp1 = sum(p[1,]);
        mu = c(0,x%*%p[2,]/sum(p[2,]));#M step for mu
        prop = apply(p,1,sum)/n1;#M step for pi
        sigma = sqrt(apply((rbind(x,x)-matrix(rep(mu,n1),nrow=2))^2*p,1,sum)/apply(p,1,sum));# M step for sigma
        sigma = mmax(0.1*max(sigma),sigma);
      }
    }

    #a = sort((mu),index.return = T)$x;
    b = sort((mu),index.return = T)$ix;
    out=list(mu=mu[b],sigma=sigma[b],pi=prop[b],lh=lh[run],p=p[b,])

  }else{
    stop=0; run=0; dif=10; difml=10; ml=10^8;

    while(stop == 0){

      denm = t(cbind(as.matrix(dnorm(oldx,mu[1],sigma[1])),as.matrix(dnorm(oldx,mu[2],sigma[2]))));
      f = apply(prop*denm,2,sum);
      #a = sort((-f),index.return = T)$x;
      b = sort((-f),index.return = T)$ix;
      x = oldx[b[1:n1]]; x = matrix(x,nrow=1);
      run = run+1;
      denm = denm[,b[1:n1]];
      f = f[b[1:n1]];
      # lh[run] = sum(log(f));#use lh to record the sequence of lh
      lh = c(lh,sum(log(f))) #updated oct23 2022 by Xin Shen to make lh a vector instead of ts to fix a bug

      if(run > 2){
        denom = lh[run-1]-lh[run-2];
        if(denom == 0){
          stop = 1;
        }else{
          c = (lh[run]-lh[run-1])/denom;#convergence rate of Em algorithm
          if(c>0 && c<1){
            preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
            dif = ml-lh[run]; difml = abs(preml-ml);
          }
        }
      }

      if(dif<0.001 && difml<0.0001 && run>20 || run>200 || is.nan(min(mu)) == TRUE){
        stop = 1;
      }else{
        p = matrix(rep(prop,n1),ncol=n1)*denm/rbind(f,f);#E-step, classification probability of x from each component
        #sp1 = sum(p[1,]);
        mu = c(0,x%*%p[2,]/sum(p[2,]));#M step for mu
        prop = apply(p,1,sum)/n1;#M step for pi
        sigma = c(sig,mmax(0.1*sig,sqrt((x-mu[2])^2%*%p[2,]/sum(p[2,]))));# M step for sigma
      }
    }

    #a = sort((mu),index.return = T)$x;
    b = sort((mu),index.return = T)$ix;
    out=list(mu=mu[b],sigma=sigma[b],pi=prop[b],lh=lh[run],p=p[b,])
  }

  return(out)
}



#uses Em algorithm to estimate k components norm mixture model with equal unknown variance
#we can also get p the classification probability matrix
#rsample;mmax
mixnveq<-function(x,k = NULL,ini = NULL){
  if(is.null(k)){
    k = 2;
  }

  x = as.vector(x);
  n = length(x); acc = 10^(-5)*n; a = ceiling(n/2);
  s2 = var(x)*(n-1)/n; sigma=matrix(rep(NaN,10),nrow=5);
  sx = sort(x); lowvar = ((range(x)[2]-range(x)[1])/38)^2;

  if(k<2){
    warning("The number of components (second input) is less than 2")
  }else if(k==2){
    #choose five intial value
    prop = 1/2*matrix(rep(1,10),nrow=5);#set intial value of pi 1/2 for each componetn
    mu = rbind(rsample(x,2),rsample(x,2),cbind(min(x),max(x)),cbind(mean(sx[1:a]),mean(sx[(a+1):n])),cbind(mean(x)-0.5*sd(x),mean(x)+0.5*sd(x)));
    sigma[(1:2),] = matrix(rep(sqrt(runif(2,0,1)*s2),2),ncol=2);
    sigma[(3:5),] = matrix(rep(sqrt(mmax(lowvar,s2-diag(var(t(mu[(3:5),])))/k)),2),ncol=2);

    if(is.null(ini) == FALSE){
      mu = rbind(mu,ini$mu);
      sigma = rbind(sigma,ini$sigma);
      prop = rbind(prop,ini$pi);
    }

    numini = dim(mu)[1];
    #run em algorithm for a while for the five intial values
    l = rep(NaN,numini);
    for(i in 1:numini){
      stop = 0; run=0; l2 = -10^8;
      while(stop==0){
        run = run+1;
        denm = rbind(dnorm(x,mu[i,1],sigma[i,1]),dnorm(x,mu[i,2],sigma[i,2]));
        f = prop[i,]%*%denm;
        l1 = l2; l2 = sum(log(f)); dif = l2-l1;

        if(dif<10^(-5)*n && run>30 || dif<10^(-8)*n){
          stop = 1;
        }else{
          p = matrix(rep(prop[i,],n),ncol=n)*denm/rbind(f,f);#E step, classification probability of x from each component
          mu[i,] = x%*%t(p)/apply(p,1,sum);#M step for mu
          prop[i,] = apply(p,1,sum)/n;#M step for prop
          sigma[i,] = sqrt(sum(sum(p*(rbind(x,x)-matrix(rep(mu[i,],n),ncol=n))^2))/n)*rep(1,2);#M step for sigma
        }
      }
      l[i] = l2;
    }

    #choose the initial value of the 4, which has the largest likelihood
    I = sort(l,index.return = T)$ix; #y = sort(l,index.return = T)$x;
    id = I[numini];
    mu = mu[id,]; prop = prop[id,]; sigma = sigma[id,];
    stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8;

    #use em algorithm to calculate the mle
    while(stop == 0){

      denm = rbind(dnorm(x,mu[1],sigma[1]),dnorm(x,mu[2],sigma[2]));
      f = prop%*%denm;
      dd = max(apply((denm-rbind(f,f))/rbind(f,f),1,sum));
      run = run+1;
      lh[run] = sum(log(f));#use lh to record the sequence of lh

      if(run > 2){
        denom = lh[run-1]-lh[run-2];
        if(denom == 0){
          stop = 1;
        }else{
          c = (lh[run]-lh[run-1])/denom;#convergence rate of Em algorithm
          if(c>0 && c<1){
            preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
            dif = ml-lh[run]; difml = abs(preml-ml);
          }
        }
      }

      #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
      if(dif<0.001 && difml<0.0001 && dd<acc){
        stop = 1;
      }else{
        p = matrix(rep(prop,n),ncol=n)*denm/rbind(f,f);#E-step, classification probability of x from each component
        #sp1 = sum(p[1,]);
        mu = x%*%t(p)/apply(p,1,sum);#M step for mu
        prop = apply(p,1,sum)/n;#M step for pi
        sigma = sqrt(sum((rbind(x,x)-matrix(rep(mu,n),nrow=2))^2*p)/n)*rep(1,2);# M step for sigma
      }
    }

    out=list(mu=mu,sigma=sigma,pi=prop,lh=lh[run],p=p)

  }else{#the situation when k>2
    mu = matrix(rep(NaN,4*k),nrow=4);
    sigma = matrix(rep(NaN,4*k),nrow=4);
    prop = 1/k*matrix(rep(1,4*k),nrow=4);
    mu[1,] = rsample(x,k); mu[2,] = rsample(x,k);
    mu[(3:4),1] = rbind(sx[1],mean(sx[1:ceiling(n/k)]));
    for(i in 1:(k-1)){
      mu[(3:4),(i+1)] = rbind(sx[ceiling(i*n/(k-1))],mean(sx[(ceiling(i*n/k)+1):ceiling((i+1)*n/k)]));
    }
    sigma[(1:2),] = matrix(rep(sqrt(runif(2,0,1)*s2),k),ncol=k);
    sigma[3,] = sqrt(mmax(lowvar,s2-var(mu[3,])*(k-1)/k))*rep(1,k);
    sigma[4,] = sqrt(mmax(lowvar,s2-var(mu[4,])*(k-1)/k))*rep(1,k);

    if(is.null(ini) == FALSE){
      mu = rbind(mu,ini$mu);
      sigma = rbind(sigma,ini$sigma);
      prop = rbind(prop,ini$pi);
    }

    numini = dim(mu)[1];
    #run em algorithm for a while for the five intial values
    l = rep(NaN,numini);
    for(i in 1:numini){
      stop = 0; run=0; l2 = -10^8;
      while(stop==0){
        run = run+1;
        denm = matrix(rep(NaN,n*k),nrow=k);
        for(j in 1:k){
          denm[j,] = dnorm(x,mu[i,j],sigma[i,j]);
        }
        f = prop[i,]%*%denm;
        l1 = l2; l2 = sum(log(f)); dif = l2-l1;

        if(dif<10^(-5)*n && run>30 || dif<10^(-8)*n){
          stop = 1;
        }else{
          p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
          mu[i,] = x%*%t(p)/apply(p,1,sum);
          prop[i,] = apply(p,1,sum)/n;
          sigma[i,] = sqrt(sum(sum(p*(t(matrix(rep(x,k),ncol=k))-matrix(rep(mu[i,],n),ncol=n))^2))/n)*rep(1,k);
        }
      }
      l[i] = l2;
    }

    #choose the initial value of the five which has the largest likelihood
    I = sort(l,index.return = T)$ix; #y = sort(l,index.return = T)$x;
    id = I[numini];
    mu = mu[id,]; prop = prop[id,]; sigma = sigma[id,];
    stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8;

    #use em algorithm to calculate the mle
    while(stop == 0){

      for(j in 1:k){
        denm[j,] = dnorm(x,mu[j],sigma[j]);
      }
      f = prop%*%denm;
      dd = max(apply((denm-t(matrix(rep(f,k),ncol=k)))/t(matrix(rep(f,k),ncol=k)),1,sum));
      run = run+1;
      lh[run] = sum(log(f));#use lh to record the sequence of lh

      if(run > 2){
        c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);#convergence rate of Em algorithm
        if(c>0 && c<1){
          preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
          dif = ml-lh[run]; difml = abs(preml-ml);
        }
      }

      #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
      if(dif<0.001 && difml<0.0001 && dd<acc){
        stop = 1;
      }else{
        p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
        mu = x%*%t(p)/apply(p,1,sum);
        prop = apply(p,1,sum)/n;
        sigma = sqrt(sum(sum(p*(t(matrix(rep(x,k),ncol=k))-matrix(rep(mu,n),ncol=n))^2))/n)*rep(1,k);
      }
    }

    out=list(mu=mu,sigma=sigma,pi=prop,lh=lh[run],p=p)
  }

  return(out)
}



#Sigma known
#ini: the initial values for mu and prop
#h:bandwidth. h=1.06*n^(-1/5) by default
#kdebw;mixnveq;mmax;interpcut
mhde2m1<-function(x,p,sigma,mu,ini = NULL){
  stopiter = 30; k = 2; n = length(x); x = matrix(x,nrow=n);
  h = kdebw(x,2^14);
  true = c(p,mu);

  if(is.null(ini)){
    ini = mixnveq(x,k);
  }

  if(ini$mu[2]>ini$mu[1]){
    prop = ini$pi[1]; mu = ini$mu[2];
  }else{
    prop = ini$pi[2]; mu = ini$mu[1];
  }
  est = c(prop,mu);

  xgridmin = min(x)-5*h; xgridmax = max(x)+5*h; lxgrid = 100;
  xgrid = seq(from=xgridmin,to=xgridmax,length.out=lxgrid);
  hspace = (xgridmax-xgridmin)/lxgrid;
  acc = 10^(-5)/hspace;

  #nonparametric estimator
  deng<-function(t){
    y = apply(exp(-(matrix(rep(x,length(t)),nrow=n)-t(matrix(rep(t,n),ncol=n)))^2/2/h^2),2,mean)/h/sqrt(2*pi);
  }
  dengx = deng(xgrid)^(1/2);#(hn_hat)^(1/2)

#calculate the MHDE using temp
  dif = acc+1; numiter = 0; fval = 10^10;

  while(dif>acc && numiter<stopiter){
    numiter = numiter+1; pfval = fval;

    denf1<-function(t){
      y = dnorm(t,0,sigma)
    }

    #Find alpha and M by iteration
    difa = 1; step = 0; a = 1;
    while(difa>10^(-3) && step<20){
      prea = a; step = step+1;

      mfun<-function(t){
        y = a*deng(t)>prop*denf1(t);
      }
      temp<-function(t){
        y = denf1(t)*mfun(t);
      }
      temp1<-function(t){
        y = deng(t)*mfun(t);
      }
      #integrate积分
      a = min((prop*integrate(temp,xgridmin,xgridmax)$value+1-prop)/max(integrate(temp1,xgridmin,xgridmax)$value,1-prop),1);
      difa = abs(prea-a);
    }
    if(a>0.99){
      a = 1;
    }

    #Given theta update f
    denfmu = (mmax(0,a*deng(xgrid)-prop*denf1(xgrid))+mmax(0,a*deng(2*mu-xgrid)-prop*denf1(2*mu-xgrid)))/2/(1-prop);
    preest=est;
    #Given f, update theta
    denf<-function(t){
      y = interpcut(c(xgrid-mu,mu-xgrid),c(denfmu,denfmu),t);
    }
    obj<-function(t){#step2
      y = sum(((min(0.95,max(t[1],0.05))*dnorm(xgrid,0,sigma)+(1-min(0.95,max(t[1],0.05)))*denf(xgrid-min(max(x),max(0,t[2]))))^(1/2)-dengx)^2);
    }

    fmin = nlminb(start=preest,objective=obj)#初始值+优化函数->非线性函数极小值
    est = fmin$par; fval = fmin$objective;
    est = c(min(est[1],0.95),min(est[2],max(x)));
    est = c(max(est[1],0.05),max(est[2],0));
    dif = pfval-fval;
    if(dif<0){
      est = preest; fval = pfval;
    }

    prop = est[1]; mu = est[2];
  }

  res=list(fval=fval,pi=prop,mu=mu,numiter=numiter)

#calculate the MHDE using true
  dif = acc+1; numiter = 0; fval = 10^10;
  est = true; prop = true[1]; mu = true[2];

  while(dif>acc && numiter<stopiter){
    numiter = numiter+1; pfval = fval;

    denf1<-function(t){
      y = dnorm(t,0,sigma)
    }

    #Find alpha and M by iteration
    difa = 1; step = 0; a = 1;
    while(difa>10^(-3) && step<20){
      prea = a; step = step+1;

      mfun<-function(t){
        y = a*deng(t)>prop*denf1(t);
      }
      temp<-function(t){
        y = denf1(t)*mfun(t);
      }
      temp1<-function(t){
        y = deng(t)*mfun(t);
      }

      a = min((prop*integrate(temp,xgridmin,xgridmax)$value+1-prop)/max(integrate(temp1,xgridmin,xgridmax)$value,1-prop),1);
      difa = abs(prea-a);
    }
    if(a>0.99){
      a = 1;
    }

    #Given theta update f
    denfmu = (mmax(0,a*deng(xgrid)-prop*denf1(xgrid))+mmax(0,a*deng(2*mu-xgrid)-prop*denf1(2*mu-xgrid)))/2/(1-prop);
    preest = est;
    #Given f, update theta
    denf<-function(t){
      y = interpcut(c(xgrid-mu,mu-xgrid),c(denfmu,denfmu),t);
    }
    obj<-function(t){
      y = sum(((min(0.95,max(t[1],0.05))*dnorm(xgrid,0,sigma)+(1-min(0.95,max(t[1],0.05)))*denf(xgrid-min(max(x),max(0,t[2]))))^(1/2)-dengx)^2);
    }

    fmin = nlminb(start=preest,objective=obj)
    est = fmin$par; fval = fmin$objective;
    est = c(min(est[1],0.95),min(est[2],max(x)));
    est = c(max(est[1],0.05),max(est[2],0));
    dif = pfval-fval;
    if(dif<0){
      est = preest; fval = pfval;
    }

    prop = est[1]; mu = est[2];
  }

  if(res$fval<fval){
    out=list(pi=res$pi,mu=res$mu,numiter=res$numiter,initialtrue=0)
  }else{
    out=list(pi=prop,mu=mu,numiter=numiter,initialtrue=1)
  }

  return(out)
}



#MHDE sigma unknown
#ini: the initial values for mu and prop.
#h:bandwidth. h=1.06*n^(-1/5) by default.
#kdebw;mixnveq;mmax;interpcut
mhdem1<-function(x,p,sigma,mu,ini = NULL){
  stopiter = 30; k = 2; n = length(x); x = matrix(x,nrow=n);
  h = kdebw(x,2^14);
  true = c(p,sigma,mu);

  if(is.null(ini)){
    ini = mixnveq(x,k);
  }

  if(ini$mu[2]>ini$mu[1]){
    prop = ini$pi[1]; mu = ini$mu[2]; sigma = ini$sigma[1];
  }else{
    prop = ini$pi[2]; mu = ini$mu[1]; sigma = ini$sigma[2];
  }
  est = c(prop,sigma,mu);

  xgridmin = min(x)-5*h; xgridmax = max(x)+5*h; lxgrid = 100;
  xgrid = seq(from=xgridmin,to=xgridmax,length.out=lxgrid);
  hspace = (xgridmax-xgridmin)/lxgrid;
  acc = 10^(-5)/hspace;

  #nonparametric estimator
  deng<-function(t){
    y = apply(exp(-(matrix(rep(x,length(t)),nrow=n)-t(matrix(rep(t,n),ncol=n)))^2/2/h^2),2,mean)/h/sqrt(2*pi);
  }
  dengx = deng(xgrid)^(1/2);

  #calculate the MHDE using temp
  dif = acc+1; numiter = 0; fval = 10^10;

  while(dif>acc && numiter<stopiter){
    numiter = numiter+1; pfval = fval;

    denf1<-function(t){
      y = dnorm(t,0,sigma)
    }

    #Find alpha and M by iteration
    difa = 1; step = 0; a = 1;
    while(difa>10^(-3) && step<20){
      prea = a; step = step+1;

      mfun<-function(t){
        y = a*deng(t)>prop*denf1(t);
      }
      temp<-function(t){
        y = denf1(t)*mfun(t);
      }
      temp1<-function(t){
        y = deng(t)*mfun(t);
      }

      a = min((prop*integrate(temp,xgridmin,xgridmax)$value+1-prop)/max(integrate(temp1,xgridmin,xgridmax)$value,1-prop),1);
      difa = abs(prea-a);
    }
    if(a>0.99){
      a = 1;
    }

    #Given theta update f
    denfmu = (mmax(0,a*deng(xgrid)-prop*denf1(xgrid))+mmax(0,a*deng(2*mu-xgrid)-prop*denf1(2*mu-xgrid)))/2/(1-prop);
    preest = est;
    #Given f, update theta
    denf<-function(t){
      y = interpcut(c(xgrid-mu,mu-xgrid),c(denfmu,denfmu),t);
    }
    obj<-function(t){
      y = sum(((min(0.95,max(t[1],0.05))*dnorm(xgrid,0,min(sd(x),max(0.1*sd(x),t[2])))+(1-min(0.95,max(t[1],0.05)))*denf(xgrid-min(max(x),max(0,t[3]))))^(1/2)-dengx)^2);
    }

    fmin = nlminb(start=preest,objective=obj)
    est = fmin$par; fval = fmin$objective;
    est = c(min(est[1],0.95),min(est[2],sd(x)),min(est[3],max(x)));
    est = c(max(est[1],0.05),max(est[2],sd(x)*0.1),max(est[3],0));
    dif = pfval-fval;
    if(dif<0){
      est = preest; fval = pfval;
    }

    prop = est[1]; sigma = est[2]; mu = est[3];
  }

  res=list(fval=fval,pi=prop,sigma=sigma,mu=mu,numiter=numiter)

  #calculate the MHDE using true
  dif = acc+1; numiter = 0; fval = 10^10;
  est = true; prop = true[1]; sigma=true[2]; mu = true[3];

  while(dif>acc && numiter<stopiter){
    numiter = numiter+1; pfval = fval;

    denf1<-function(t){
      y = dnorm(t,0,sigma)
    }

    #Find alpha and M by iteration
    difa = 1; step = 0; a = 1;
    while(difa>10^(-3) && step<20){
      prea = a; step = step+1;

      mfun<-function(t){
        y = a*deng(t)>prop*denf1(t);
      }
      temp<-function(t){
        y = denf1(t)*mfun(t);
      }
      temp1<-function(t){
        y = deng(t)*mfun(t);
      }

      a = min((prop*integrate(temp,xgridmin,xgridmax)$value+1-prop)/max(integrate(temp1,xgridmin,xgridmax)$value,1-prop),1);
      difa = abs(prea-a);
    }
    if(a>0.99){
      a = 1;
    }

    #Given theta update f
    denfmu = (mmax(0,a*deng(xgrid)-prop*denf1(xgrid))+mmax(0,a*deng(2*mu-xgrid)-prop*denf1(2*mu-xgrid)))/2/(1-prop);
    preest = est;
    #Given f, update theta
    denf<-function(t){
      y = interpcut(c(xgrid-mu,mu-xgrid),c(denfmu,denfmu),t);
    }
    obj<-function(t){
      y = sum(((min(0.95,max(t[1],0.05))*dnorm(xgrid,0,min(sd(x),max(0.1*sd(x),t[2])))+(1-min(0.95,max(t[1],0.05)))*denf(xgrid-min(max(x),max(0,t[3]))))^(1/2)-dengx)^2);
    }

    fmin = nlminb(start=preest,objective=obj)
    est = fmin$par; fval = fmin$objective;
    est = c(min(est[1],0.95),min(est[2],sd(x)),min(est[3],max(x)));
    est = c(max(est[1],0.05),max(est[2],sd(x)*0.1),max(est[3],0));
    dif = pfval-fval;
    if(dif<0){
      est = preest; fval = pfval;
    }

    prop = est[1]; sigma = est[2]; mu = est[3];
  }

  if(res$fval<fval){
    out=list(pi=res$pi,sigma=res$sigma,mu=res$mu,numiter=res$numiter,initialtrue=0)
  }else{
    out=list(pi=prop,sigma=sigma,mu=mu,numiter=numiter,initialtrue=1)
  }

  return(out)
}




#true = c(p,mu);
mhde2m1_x<-function(x,sigma,ini = NULL,true = NULL){
  stopiter = 30; k = 2; n = length(x); x = matrix(x,nrow=n);
  h = kdebw(x,2^14);

  if(is.null(ini)){
    ini = mixnveq(x,k);
  }

  if(ini$mu[2]>ini$mu[1]){
    prop = ini$pi[1]; mu = ini$mu[2];
  }else{
    prop = ini$pi[2]; mu = ini$mu[1];
  }
  est = c(prop,mu);

  xgridmin = min(x)-5*h; xgridmax = max(x)+5*h; lxgrid = 100;
  xgrid = seq(from=xgridmin,to=xgridmax,length.out=lxgrid);
  hspace = (xgridmax-xgridmin)/lxgrid;
  acc = 10^(-5)/hspace;

  #nonparametric estimator
  deng<-function(t){
    y = apply(exp(-(matrix(rep(x,length(t)),nrow=n)-t(matrix(rep(t,n),ncol=n)))^2/2/h^2),2,mean)/h/sqrt(2*pi);
  }
  dengx = deng(xgrid)^(1/2);

  #calculate the MHDE using temp
  dif = acc+1; numiter = 0; fval = 10^10;

  while(dif>acc && numiter<stopiter){
    numiter = numiter+1; pfval = fval;

    denf1<-function(t){
      y = dnorm(t,0,sigma)
    }

    #Find alpha and M by iteration
    difa = 1; step = 0; a = 1;
    while(difa>10^(-3) && step<20){
      prea = a; step = step+1;

      mfun<-function(t){
        y = a*deng(t)>prop*denf1(t);
      }
      temp<-function(t){
        y = denf1(t)*mfun(t);
      }
      temp1<-function(t){
        y = deng(t)*mfun(t);
      }

      a = min((prop*integrate(temp,xgridmin,xgridmax)$value+1-prop)/max(integrate(temp1,xgridmin,xgridmax)$value,1-prop),1);
      difa = abs(prea-a);
    }
    if(a>0.99){
      a = 1;
    }

    #Given theta update f
    denfmu = (mmax(0,a*deng(xgrid)-prop*denf1(xgrid))+mmax(0,a*deng(2*mu-xgrid)-prop*denf1(2*mu-xgrid)))/2/(1-prop);
    preest=est;
    #Given f, update theta
    denf<-function(t){
      y = interpcut(c(xgrid-mu,mu-xgrid),c(denfmu,denfmu),t);
    }
    obj<-function(t){
      y = sum(((min(0.95,max(t[1],0.05))*dnorm(xgrid,0,sigma)+(1-min(0.95,max(t[1],0.05)))*denf(xgrid-min(max(x),max(0,t[2]))))^(1/2)-dengx)^2);
    }

    fmin = nlminb(start=preest,objective=obj)
    est = fmin$par; fval = fmin$objective;
    est = c(min(est[1],0.95),min(est[2],max(x)));
    est = c(max(est[1],0.05),max(est[2],0));
    dif = pfval-fval;
    if(dif<0){
      est = preest; fval = pfval;
    }

    prop = est[1]; mu = est[2];
  }

  res=list(fval=fval,pi=prop,mu=mu,numiter=numiter)
  out=list(fval=fval,pi=prop,mu=mu,numiter=numiter)

  if(is.null(true)){
    out = out
  }else{#calculate the MHDE using true
    dif = acc+1; numiter = 0; fval = 10^10;
    est = true; prop = true[1]; mu = true[2];

    while(dif>acc && numiter<stopiter){
      numiter = numiter+1; pfval = fval;

      denf1<-function(t){
        y = dnorm(t,0,sigma)
      }

      #Find alpha and M by iteration
      difa = 1; step = 0; a = 1;
      while(difa>10^(-3) && step<20){
        prea = a; step = step+1;

        mfun<-function(t){
          y = a*deng(t)>prop*denf1(t);
        }
        temp<-function(t){
          y = denf1(t)*mfun(t);
        }
        temp1<-function(t){
          y = deng(t)*mfun(t);
        }

        a = min((prop*integrate(temp,xgridmin,xgridmax)$value+1-prop)/max(integrate(temp1,xgridmin,xgridmax)$value,1-prop),1);
        difa = abs(prea-a);
      }
      if(a>0.99){
        a = 1;
      }

      #Given theta update f
      denfmu = (mmax(0,a*deng(xgrid)-prop*denf1(xgrid))+mmax(0,a*deng(2*mu-xgrid)-prop*denf1(2*mu-xgrid)))/2/(1-prop);
      preest = est;
      #Given f, update theta
      denf<-function(t){
        y = interpcut(c(xgrid-mu,mu-xgrid),c(denfmu,denfmu),t);
      }
      obj<-function(t){
        y = sum(((min(0.95,max(t[1],0.05))*dnorm(xgrid,0,sigma)+(1-min(0.95,max(t[1],0.05)))*denf(xgrid-min(max(x),max(0,t[2]))))^(1/2)-dengx)^2);
      }

      fmin = nlminb(start=preest,objective=obj)
      est = fmin$par; fval = fmin$objective;
      est = c(min(est[1],0.95),min(est[2],max(x)));
      est = c(max(est[1],0.05),max(est[2],0));
      dif = pfval-fval;
      if(dif<0){
        est = preest; fval = pfval;
      }

      prop = est[1]; mu = est[2];
    }

    if(res$fval<fval){
      out=list(pi=res$pi,mu=res$mu,numiter=res$numiter,initialtrue=0)
    }else{
      out=list(pi=prop,mu=mu,numiter=numiter,initialtrue=1)
    }
  }

  return(out)
}

#true = c(p,sigma,mu);
mhdem1_x<-function(x,ini = NULL,true = NULL){
  stopiter = 30; k = 2; n = length(x); x = matrix(x,nrow=n);
  h = kdebw(x,2^14);

  if(is.null(ini)){
    ini = mixnveq(x,k);
  }

  if(ini$mu[2]>ini$mu[1]){
    prop = ini$pi[1]; mu = ini$mu[2]; sigma = ini$sigma[1];
  }else{
    prop = ini$pi[2]; mu = ini$mu[1]; sigma = ini$sigma[2];
  }
  est = c(prop,sigma,mu);

  xgridmin = min(x)-5*h; xgridmax = max(x)+5*h; lxgrid = 100;
  xgrid = seq(from=xgridmin,to=xgridmax,length.out=lxgrid);
  hspace = (xgridmax-xgridmin)/lxgrid;
  acc = 10^(-5)/hspace;

  #nonparametric estimator
  deng<-function(t){
    y = apply(exp(-(matrix(rep(x,length(t)),nrow=n)-t(matrix(rep(t,n),ncol=n)))^2/2/h^2),2,mean)/h/sqrt(2*pi);
  }
  dengx = deng(xgrid)^(1/2);

  dif = acc+1; numiter = 0; fval = 10^10;

  while(dif>acc && numiter<stopiter){
    numiter = numiter+1; pfval = fval;

    denf1<-function(t){
      y = dnorm(t,0,sigma)
    }

    #Find alpha and M by iteration
    difa = 1; step = 0; a = 1;
    while(difa>10^(-3) && step<20){
      prea = a; step = step+1;

      mfun<-function(t){
        y = a*deng(t)>prop*denf1(t);
      }
      temp<-function(t){
        y = denf1(t)*mfun(t);
      }
      temp1<-function(t){
        y = deng(t)*mfun(t);
      }

      a = min((prop*integrate(temp,xgridmin,xgridmax)$value+1-prop)/max(integrate(temp1,xgridmin,xgridmax)$value,1-prop),1);
      difa = abs(prea-a);
    }
    if(a>0.99){
      a = 1;
    }

    #Given theta update f
    denfmu = (mmax(0,a*deng(xgrid)-prop*denf1(xgrid))+mmax(0,a*deng(2*mu-xgrid)-prop*denf1(2*mu-xgrid)))/2/(1-prop);
    preest = est;
    #Given f, update theta
    denf<-function(t){
      y = interpcut(c(xgrid-mu,mu-xgrid),c(denfmu,denfmu),t);
    }
    obj<-function(t){
      y = sum(((min(0.95,max(t[1],0.05))*dnorm(xgrid,0,min(sd(x),max(0.1*sd(x),t[2])))+(1-min(0.95,max(t[1],0.05)))*denf(xgrid-min(max(x),max(0,t[3]))))^(1/2)-dengx)^2);
    }

    fmin = nlminb(start=preest,objective=obj)
    est = fmin$par; fval = fmin$objective;
    est = c(min(est[1],0.95),min(est[2],sd(x)),min(est[3],max(x)));
    est = c(max(est[1],0.05),max(est[2],sd(x)*0.1),max(est[3],0));
    dif = pfval-fval;
    if(dif<0){
      est = preest; fval = pfval;
    }

    prop = est[1]; sigma = est[2]; mu = est[3];
  }

  res=list(fval=fval,pi=prop,sigma=sigma,mu=mu,numiter=numiter)
  out=list(fval=fval,pi=prop,sigma=sigma,mu=mu,numiter=numiter)

  if(is.null(true)){
    out = out
  }else{
    dif = acc+1; numiter = 0; fval = 10^10;
    est = true; prop = true[1]; sigma=true[2]; mu = true[3];

    while(dif>acc && numiter<stopiter){
      numiter = numiter+1; pfval = fval;

      denf1<-function(t){
        y = dnorm(t,0,sigma)
      }

      #Find alpha and M by iteration
      difa = 1; step = 0; a = 1;
      while(difa>10^(-3) && step<20){
        prea = a; step = step+1;

        mfun<-function(t){
          y = a*deng(t)>prop*denf1(t);
        }
        temp<-function(t){
          y = denf1(t)*mfun(t);
        }
        temp1<-function(t){
          y = deng(t)*mfun(t);
        }

        a = min((prop*integrate(temp,xgridmin,xgridmax)$value+1-prop)/max(integrate(temp1,xgridmin,xgridmax)$value,1-prop),1);
        difa = abs(prea-a);
      }
      if(a>0.99){
        a = 1;
      }

      #Given theta update f
      denfmu = (mmax(0,a*deng(xgrid)-prop*denf1(xgrid))+mmax(0,a*deng(2*mu-xgrid)-prop*denf1(2*mu-xgrid)))/2/(1-prop);
      preest = est;
      #Given f, update theta
      denf<-function(t){
        y = interpcut(c(xgrid-mu,mu-xgrid),c(denfmu,denfmu),t);
      }
      obj<-function(t){
        y = sum(((min(0.95,max(t[1],0.05))*dnorm(xgrid,0,min(sd(x),max(0.1*sd(x),t[2])))+(1-min(0.95,max(t[1],0.05)))*denf(xgrid-min(max(x),max(0,t[3]))))^(1/2)-dengx)^2);
      }

      fmin = nlminb(start=preest,objective=obj)
      est = fmin$par; fval = fmin$objective;
      est = c(min(est[1],0.95),min(est[2],sd(x)),min(est[3],max(x)));
      est = c(max(est[1],0.05),max(est[2],sd(x)*0.1),max(est[3],0));
      dif = pfval-fval;
      if(dif<0){
        est = preest; fval = pfval;
      }

      prop = est[1]; sigma = est[2]; mu = est[3];
    }

    if(res$fval<fval){
      out=list(pi=res$pi,sigma=res$sigma,mu=res$mu,numiter=res$numiter,initialtrue=0)
    }else{
      out=list(pi=prop,sigma=sigma,mu=mu,numiter=numiter,initialtrue=1)
    }
  }

  return(out)
}
