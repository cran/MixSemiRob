#' Complete Likelihood Frequency Method for Label Switching
#'
#' Solving the label switching problem by maximizing the complete likelihood.
#'
#' @param est Records the MLE and the estimated classification probabilities, can
#' be derived by the function: est=mixnorm(x,2,ini).
#' @param lat The latent component labels for all observations, a m by n 0-1 matrix
#' if (i,j)th cell is 1, then jth observation is from ith component.
#'
#' @return
#' The estimation by complete posterior distribution
#' the result vector consist of mu of first dimension of first component, second component...,
#' mu of second dimension of first component, second component ...
#' pi of first component , second component ...
#'
#' @export
#'@importFrom mvtnorm rmvnorm dmvnorm
#' @examples
#'
#' #Example1 : Mixture of Two Univariate Normal
#'
#' #Simulate the data
#' n=200;prop=0.3;mudif=1.5;
#' n1=rbinom(1,n,prop);
#' x1=rnorm(n1,0,1);x2=rnorm(n-n1,mudif,1);
#' x=c(x1,x2);
#' pm=c(2,1,3,5,4);
#'
#' #Use mixnorm to get MLE and the estimated classification probabilities
#' out=mixnorm(x,2);
#'
#' #Prepare latent component label
#' lat=rbind(cbind(matrix(rep(1,n1),nrow=1),matrix(rep(0,n-n1),nrow=1)),
#' cbind(matrix(rep(0,n1),nrow=1),matrix(rep(1,n-n1),nrow=1)));
#'
#' #Fit the complhfrequency/ distlatfrequency function
#' clhest=complhfrequency(out,lat);
#' ditlatest=distlatfrequency(out,lat);
#'
#' #Order constraint labeling based on mu
#' est=c(out$mu,out$sigma[1],out$pi);
#' if(est[1]<est[2]){
#'  ordmuest=est;
#' }else{
#'  ordmuest=est[pm];
#'  }
#'
#' #Order constraint labeling based on pi
#' if(est[4]<est[5]){
#'  ordpiest=est;
#'   }else{
#'  ordpiest=est[pm];
#'  }
#'
#'#Example 2: Mixture of Two Multivariate Normal
#'
#' #Simulate the data
#'  mu1=0.5;mu2=0.5;
#'  prop=0.3;
#'  n=400;n1=rbinom(1,n,prop);
#'  pm=c(2,1,4,3,6,5);
#'  sigma=diag(c(1,1));mu=matrix(c(0,mu1,0,mu2),ncol=2);pi=c(prop,1-prop);
#'  ini=list(sigma=sigma,mu=mu,pi=pi)
#'  x1=mvtnorm::rmvnorm(n1,c(0,0),ini$sigma);x2=mvtnorm::rmvnorm(n-n1,c(mu1,mu2),ini$sigma);
#'  x=rbind(x1,x2);
#'
#' #Use mixnorm to get MLE and the estimated classification probabilities
#'  out=mixmnorm(x,2,1,ini);
#'
#' #Prepare latent component label
#'  lat=rbind(cbind(matrix(rep(1,n1),nrow=1),matrix(rep(0,n-n1),nrow=1)),
#'  cbind(matrix(rep(0,n1),nrow=1),matrix(rep(1,n-n1),nrow=1)));
#'
#' #Fit the complhfrequency/ distlatfrequency function
#'  clhest=complhfrequency(out,lat);
#'  distlatest=distlatfrequency(out,lat);
#'
#' #Order constraint labeling based on mu1 or mu2
#'  est=c(out$mu[,1],out$mu[,2],out$pi,out$sigma[1,1],
#'  out$sigma[2,1],out$sigma[2,2]);
#'  if(est[1]<est[2]){
#'  ordmu1est=est[1:6];
#'  }else{
#'  ordmu1est=est[pm];
#'  }
#'  if(est[3]<est[4]){
#'   ordmu2est=est[1:6];
#'  }else{
#'   ordmu2est=est[pm];
#'  }
#'
#' #Order constraint labeling based on pi
#'  if(est[5]<est[6]){
#'   ordpiest=est[1:6];
#'  }else{
#'   ordpiest=est[pm];
#'   }
#'
#'
#' #Example 3 Mixture of Three Multivariate Normal
#'
#' #Simulate the data
#' mu1=1;mu2=1;
#' prop=c(0.2,0.3,0.5);n=100;m=3;
#' sigma=diag(c(1,1));mu=matrix(c(0,mu1,2*mu1,0,mu2,2*mu2),ncol=2);pi=prop;
#' ini=list(sigma=sigma,mu=mu,pi=pi)
#' u=runif(n,0,1);n1=sum(u<prop[1]);n3=sum(u>(prop[1]+prop[2]));n2=n-n1-n3;
#' x1=mvtnorm::rmvnorm(n1,c(0,0),ini$sigma);
#' x2=mvtnorm::rmvnorm(n2,c(mu1,mu2),ini$sigma);
#' x3=mvtnorm::rmvnorm(n3,c(2*mu1,2*mu2),ini$sigma);
#' x=rbind(x1,x2,x3);
#'
#' #Use mixnorm to get MLE and the estimated classification probabilities
#' out=mixmnorm(x,m,1,ini);
#'
#' #Prepare latent component label
#' lat=rbind(cbind(matrix(rep(1,n1),nrow=1),matrix(rep(0,n-n1),nrow=1)),
#' cbind(matrix(rep(0,n1),nrow=1),matrix(rep(1,n2),nrow=1),matrix(rep(0,n3),nrow=1)),
#' cbind(matrix(rep(0,n-n3),nrow=1),matrix(rep(1,n3),nrow=1)));
#'
#' #Fit the complhfrequency/ distlatfrequency function
#' clhest=complhfrequency(out,lat);
#' distlatest=distlatfrequency(out,lat);
#'
#' #Order constraint labeling based on mu1 or mu2
#' est=c(out$mu[,1],out$mu[,2],out$pi);
#' b=sort(est[1:3],index.return = TRUE)$ix; pm=c(b,b+3,b+6);
#' ordmu1est=est[pm];
#' b=sort(est[4:6],index.return = TRUE)$ix; pm=c(b,b+3,b+6);
#' ordmu2est=est[pm];
#'
#' #Order constraint labeling based on pi
#' b=sort(est[7:9],index.return = TRUE)$ix; pm=c(b,b+3,b+6);
#' ordpiest=est[pm];
#'
complhfrequency<-function(est,lat){
  s = dim(est$mu);
  if(s[1]==1){ #one dimensional data
    m = length(est$mu); #the number of components
    pm = permfull(1:m); numpm = max(dim(pm));
    p = est$p;

    if(est$sigma[1]!=est$sigma[2]){ #the variance is unequal
      est = c(est$mu,est$sigma,est$pi);
      pm = cbind(pm,pm+m,pm+2*m);
      dist = numeric();
      for(j in 1:numpm){
        dist[j] = sum(log(p[pm[j,1:m],])*lat);
      }
      ind = sort(dist,index.return = T)$ix;
      if(ind[numpm]!=1){
        est = est[pm[ind[numpm],]];
      }
    }else{ #the variance is equal
      est = c(est$mu,est$sigma[1],est$pi);
      pm = cbind(pm,matrix(rep(m+1,numpm),ncol=1),pm+m+1);
      dist = numeric();
      for(j in 1:numpm){
        dist[j] = sum(log(p[pm[j,1:m],])*lat);
      }
      ind = sort(dist,index.return = T)$ix;
      if(ind[numpm]!=1){
        est = est[pm[ind[numpm],]];
      }
    }
  }else{ #high dimensional data equal covariance
    m = length(est$pi);
    p = est$p;
    pm = permfull(1:m); numpm = max(dim(pm));
    s=s[2];
    temp=t(est$mu[,1])
    for(i in 2:s){
      temp = cbind(temp,t(est$mu[,i]));
    }
    temp = c(temp,est$pi);
    est = temp;
    temp=pm;
    for(i in 2:(s+1)){
      temp = cbind(temp,pm+(i-1)*m);
    }
    pm = temp;
    dist = numeric();
    for(j in 1:numpm){
      dist[j] = sum(log(p[pm[j,1:m],])*lat);
    }
    ind = sort(dist,index.return = T)$ix;
    if(ind[numpm]!=1){
      est = est[pm[ind[numpm],]];
    }
  }

  return(est)
}


#' Euclidean Distance Based Labeling Method for Label Switching
#'
#' Solving the label switching problem  by minimizing the distance between the
#' classification probabilities and the true latent labels.
#'
#' @param est Records the MLE and the estimated classification probabilities, can
#' be derived by the function: est=mixnorm(x,2,ini).
#' @param lat The latent component labels for all observations, a m by n 0-1 matrix
#' lat if (i,j)th cell is 1, then jth observation is from ith component.
#'
#' @return
#' The estimation by complete posterior distribution
#' the result vector consist of mu of first dimension of first component, second component...,
#' mu of second dimension of first component, second component ...
#' pi of first component , second component ...
#' @export
#'
#' @examples
#' #See examples in complhfrequency function.
distlatfrequency<-function(est,lat){
  s = dim(est$mu);
  if(s[1]==1){ #one dimensional data
    m = length(est$mu); #the number of components
    pm = permfull(1:m); numpm = max(dim(pm));
    p = est$p;

    if(est$sigma[1]!=est$sigma[2]){ #the variance is unequal
      est = c(est$mu,est$sigma,est$pi);
      pm = cbind(pm,pm+m,pm+2*m);
      dist = numeric();
      for(j in 1:numpm){
        dist[j] = sum(log(p[pm[j,1:m],])*lat);
      }
      ind = sort(dist,index.return = T)$ix;
      if(ind[numpm]!=1){
        est = est[pm[ind[numpm],]];
      }
    }else{ #the variance is equal
      est = c(est$mu,est$sigma[1],est$pi);
      pm = cbind(pm,matrix(rep(m+1,numpm),ncol=1),pm+m+1);
      dist = numeric();
      for(j in 1:numpm){
        dist[j] = sum(log(p[pm[j,1:m],])*lat);
      }
      ind = sort(dist,index.return = T)$ix;
      if(ind[numpm]!=1){
        est = est[pm[ind[numpm],]];
      }
    }
  }else{ #high dimensional data equal covariance
    m = length(est$pi);
    p = est$p;
    pm = permfull(1:m); numpm = max(dim(pm));
    s=s[2];
    temp=t(est$mu[,1])
    for(i in 2:s){
      temp = cbind(temp,t(est$mu[,i]));
    }
    temp = c(temp,est$pi);
    est = temp;
    temp=pm;
    for(i in 2:(s+1)){
      temp = cbind(temp,pm+(i-1)*m);
    }
    pm = temp;
    dist = numeric();
    for(j in 1:numpm){
      dist[j] = sum(log(p[pm[j,1:m],])*lat);
    }
    ind = sort(dist,index.return = T)$ix;
    if(ind[numpm]!=1){
      est = est[pm[ind[numpm],]];
    }
  }

  out=list(est=est,obj=dist[ind[numpm]],listobj=dist)
  return(out)
}



#' Normal Likelihood Based Method for Label Switching
#'
#'
#' Solving the label switching problem by using Em algorithm to estimate normal
#' mixture model with equal variance. The function support one dimensional and high
#' dimensional observations.
#'
#' @param x The observations. See details below.
#' @param k k components (2 by default); k can also be a vector. See details below.
#' @param acc Stopping criteria. When x is one dimensional, the default value of acc is
#'  0.001; When x is multidimensional the stopping criteria is 10^-5.
#' @param ini List of initial value of mu, pi, sigma (when k>2). \cr Example: ini = list(mu=c(0,1),pi=c(0.5,0.5), sigma=c(1,1)),
#'  If left as NULL, initial values will be calculated from random sample from observations x by default.
#'
#'@details
#'If the observation is one dimensional
#'
#'When k (k=2 by default) is a scalar:
#'Input x: the data vector, k: k components (2 by default);
#'Output: 3 by k matrix, each column corresponding to one component and 1st line:mu,2nd:pi,3rd:sigma.
#'
#'When k is a vector:
#'Input k:the standard deviation vector of k=length(k) components.
#'Output: 2 by k matrix, each column corresponding to one component and 1st line:mu, 2nd:pi.
#'
#'If the observation is high dimensional
#'
#'When k (k=2 by default) is a scalar:
#'Input x: the data matrix with each row corresponding to a observation, the dimension of data s is column dimension of x
#'Output: (k+s+k) by s matrix, first k rows are component means, next s rows are sigma matrix, the last k rows are pi estimation
#'
#'When k is a vector:
#'Input k:the standard deviation matrix of all the components with number of components k=k/s
#'Output: (k+k) by s matrix,first k rows are component means, next k rows are pi estimation
#'
#'@importFrom mvtnorm dmvnorm rmvnorm
#'
#' @return List of estimations:
#' mu: estimated means of each group;
#' pi: estimated proportion of each group;
#' sigma: estimated sd of each group;
#' p: the classification probability matrix;
#' lh: the likelihood value of the estimations;
#' id: which method we use for initial value;
#' inilh: likelihood of five random starting point. The program choose the point with largest
#' likelihood to continue run the estimation process until converge;
#' iniid: Which one of the initial point has the largest likelihood;
#'
#' @export
#'
#' @examples
#' #See example in complhfrequency
mixnorm<-function(x,k = NULL,ini = NULL,acc = NULL){
  x = as.matrix(x);
  s = dim(x);

  if(min(s)>1){
    if(is.null(k)){
      k = 2; equal = 1;
    }else{
      k = as.matrix(k);
      temp = dim(k);
      if(max(temp)==1){
        equal = 1;
      }else{
        equal = 2;
      }
    }
    main = mixmnorm(x,k,equal,ini);
  }else{
    if(is.null(acc)){
      acc = 0.001;
    }
    if(is.null(k)){
      k = 2;
    }
    if(s[1]>1){
      x=t(x);
    }
    if(length(k)>1){ #The situation when the variance is known
      sigma = k; k = length(sigma); n = length(x); acc=10^(-5)*n;
      sx=sort(x);

      if(k==2){
        prop = matrix(rep(1/2,10*2),nrow=10);
        stdmix = sqrt(max(10^(-4)*var(x),var(x)-mean(sigma^2)));
        mu = rbind(sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),cbind(mean(x)-stdmix,mean(x)+stdmix),cbind(min(x),max(x)),cbind(mean(sx[1:round(n/2+0.0001)]),mean(sx[(round(n/2)+1+0.0001):n])),cbind(mean(x)-0.5*sd(x),mean(x)+0.5*sd(x)));
        if(is.null(ini)==FALSE){
          mu = rbind(mu,ini$mu);
          prop = rbind(prop,ini$pi);
        }

        numini = dim(mu)[1]; l = rep(NaN,numini);
        for(i in 1:numini){
          stop = 0; run=0; l2 = -10^8;
          while(stop==0){
            run = run+1;
            denm = rbind(dnorm(x,mu[i,1],sigma[1]),dnorm(x,mu[i,2],sigma[2]));
            f = prop[i,]%*%denm;
            l1 = l2; l2 = sum(log(f)); dif = l2-l1;

            if(dif<10^(-4)*n && run>15 || dif<10^(-8)*n || run>1000){
              stop = 1;
            }else{
              p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,2),ncol=2));
              mu[i,] = x%*%t(p)/apply(p,1,sum);
              prop[i,] = apply(p,1,sum)/n;
            }
          }
          l[i] = l2;
        }

        I = sort(l,index.return = T)$ix;
        id = I[numini];
        mu = mu[id,]; prop = prop[id,];
        stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8;

        while(stop == 0){
          denm = rbind(dnorm(x,mu[1],sigma[1]),dnorm(x,mu[2],sigma[2]));
          f = prop%*%denm;
          run = run+1;
          lh[run] = sum(log(f));

          if(run > 2){


            if(lh[run-1]-lh[run-2] == 0){
              stop = 1;
            }else{
              c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);
              if(c>0 && c<1){
                preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);
                dif = ml-lh[run]; difml = abs(preml-ml);
              }
            }
          }

          if(dif<acc && difml<10^(-4) || run>50000){
            stop = 1;
          }else{
            p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,2),ncol=2));
            mu[i,] = x%*%t(p)/apply(p,1,sum);
            prop[i,] = apply(p,1,sum)/n;
          }
        }

        main=list(mu=mu,pi=prop,inilh=l,lh=lh[run],p=p,iniid=id)

      }else{ #the situation when k>2
        prop = matrix(rep(1/k,10*k),nrow=10);
        mu = rbind(sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)));
        mmu = rbind(sx[1],mean(sx[1:round(n/k+0.0001)]));
        for(i in 1:(k-1)){
          mmu = cbind(mmu,rbind(sx[round(i*n/(k-1)+0.0001)],mean(sx[(round(i*n/k)+1+0.0001):round((i+1)*n/k+0.0001)])));
        }
        mu = rbind(mu,mmu);
        if(is.null(ini)==FALSE){
          mu = rbind(mu,ini$mu);
          prop = rbind(prop,ini$pi);
        }

        numini = dim(mu)[1]; l = rep(NaN,numini);
        for(i in 1:numini){
          stop = 0; run=0; l2 = -10^8;
          while(stop==0){
            run = run+1;
            denm = matrix(rep(0,n*k),nrow=k);
            for(j in 1:k){
              denm[j,] = dnorm(x,mu[i,j],sigma[j]);
            }
            f = prop[i,]%*%denm;
            l1 = l2; l2 = sum(log(f)); dif = l2-l1;

            if(dif<10^(-4)*n && run>15 || dif<10^(-8)*n || run>1000){
              stop = 1;
            }else{
              p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
              mu[i,] = x%*%t(p)/apply(p,1,sum);
              prop[i,] = apply(p,1,sum)/n;
            }
          }
          l[i] = l2;
        }

        I = sort(l,index.return = T)$ix;
        id = I[numini];
        mu = mu[id,]; prop = prop[id,];
        stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8;

        while(stop == 0){
          for(j in 1:k){
            denm[j,] = dnorm(x,mu[j],sigma[j]);
          }
          f = prop%*%denm;
          run = run+1;
          lh[run] = sum(log(f));

          if(run > 2){


            if(lh[run-1]-lh[run-2] == 0){
              stop = 1;
            }else{
              c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);
              if(c>0 && c<1){
                preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);
                dif = ml-lh[run]; difml = abs(preml-ml);
              }
            }
          }

          if(dif<acc && difml<10^(-4) || run>50000){
            stop = 1;
          }else{
            p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
            mu = x%*%t(p)/apply(p,1,sum);
            prop = apply(p,1,sum)/n;
          }
        }

        main=list(mu=mu,pi=prop,inilh=l,lh=lh[run],p=p,iniid=id)
      }

    }else{ #when the equal variance is unknown
      k = k; n = length(x); acc = acc*n; a = round(n/2+0.0001);#Xin Shen updated acc from 10^(-5) to acc
      s2 = var(as.numeric(x))*(n-1)/n;
      sx = sort(x);
      lowvar = ((range(x)[2]-range(x)[1])/35)^2;

      if(k==2){
        prop = matrix(rep(1/2,10*2),nrow=10);
        mu = rbind(sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),sort(rsample(x,2)),cbind(min(x),max(x)),cbind(mean(sx[1:a]),mean(sx[(a+1+0.0001):n])),cbind(mean(x)-0.5*sd(x),mean(x)+0.5*sd(x)));
        sigma = matrix(rep(sqrt(runif(7,0,1)*s2),2),nrow=7);
        sigma = rbind(sigma,matrix(rep(sqrt(mmax(lowvar,s2-diag(var(t(mu[(8:10),])))/k)),2),ncol=2));

        if(is.null(ini)==FALSE){
          mu = rbind(mu,ini$mu);
          prop = rbind(prop,ini$pi);
          sigma = rbind(sigma,ini$sigma);
        }

        numini = dim(mu)[1]; l = rep(NaN,numini);
        for(i in 1:numini){
          stop = 0; run=0; l2 = -10^8;
          while(stop==0){
            run = run+1;
            denm = rbind(dnorm(x,mu[i,1],sigma[i,1]),dnorm(x,mu[i,2],sigma[i,2]));
            f = prop[i,]%*%denm;
            l1 = l2; l2 = sum(log(f)); dif = l2-l1;

            if(dif<10^(-4)*n && run>30 || dif<10^(-8)*n || run>1000){
              stop = 1;
            }else{
              p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,2),ncol=2));
              mu[i,] = x%*%t(p)/apply(p,1,sum);
              prop[i,] = apply(p,1,sum)/n;
              sigma[i,] = sqrt(sum(sum(p*(rbind(x,x)-matrix(rep(mu[i,],n),ncol=n))^2))/n)*rep(1,2);
            }
          }
          l[i] = l2;
        }

        I = sort(l,index.return = T)$ix;
        id = I[numini];
        mu = mu[id,]; prop = prop[id,]; sigma = sigma[id,];
        stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8;

        while(stop == 0){
          denm = rbind(dnorm(x,mu[1],sigma[1]),dnorm(x,mu[2],sigma[2]));
          f = prop%*%denm;
          run = run+1;
          lh[run] = sum(log(f));

          if(run > 2){

            if(lh[run-1]-lh[run-2] == 0){
              stop = 1;
            }else{
              c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);
              if(c>0 && c<1){
                preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);
                dif = ml-lh[run]; difml = abs(preml-ml);
              }
            }
          }

          if(dif<acc && difml<10^(-4) || run>50000){
            stop = 1;
          }else{
            p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,2),ncol=2));
            mu = x%*%t(p)/apply(p,1,sum);
            prop = apply(p,1,sum)/n;
            sigma = sqrt(sum((rbind(x,x)-matrix(rep(mu,n),nrow=2))^2*p)/n)*rep(1,2);
          }
        }

        main=list(mu=mu,pi=prop,sigma=sigma,inilh=l,lh=lh[run],p=p,iniid=id)

      }else{ #the situation when k>2
        prop = matrix(rep(1/k,10*k),nrow=10);
        mu = rbind(sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)),sort(rsample(x,k)));
        mmu = rbind(sx[1],mean(sx[1:round(n/k+0.0001)]));
        for(i in 1:(k-1)){
          mmu = cbind(mmu,rbind(sx[round(i*n/(k-1)+0.0001)],mean(sx[(round(i*n/k)+1+0.0001):round((i+1)*n/k+0.0001)])));
        }
        mu = rbind(mu,mmu);
        sigma = matrix(rep(sqrt(runif(8,0,1)*s2),k),nrow=8);
        sigma = rbind(sigma,matrix(rep(sqrt(mmax(lowvar,s2-var(mu[9,])*(k-1)/k)),k),ncol=k));
        sigma = rbind(sigma,matrix(rep(sqrt(mmax(lowvar,s2-var(mu[10,])*(k-1)/k)),k),ncol=k));

        if(is.null(ini)==FALSE){
          mu = rbind(mu,ini$mu);
          prop = rbind(prop,ini$pi);
          sigma = rbind(sigma,ini$sigma);
        }

        numini = dim(mu)[1]; l = rep(NaN,numini);
        for(i in 1:numini){
          stop = 0; run=0; l2 = -10^8;
          while(stop==0){
            run = run+1;
            denm = matrix(rep(0,n*k),nrow=k);
            for(j in 1:k){
              denm[j,] = dnorm(x,mu[i,j],sigma[i,j]);
            }
            f = prop[i,]%*%denm;
            l1 = l2; l2 = sum(log(f)); dif = l2-l1;

            if(dif<10^(-4)*n && run>15 || dif<10^(-8)*n || run>1000){
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

        I = sort(l,index.return = T)$ix;
        id = I[numini];
        mu = mu[id,]; prop = prop[id,]; sigma = sigma[id,];
        stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8;lh=numeric()

        while(stop == 0){
          for(j in 1:k){
            denm[j,] = dnorm(x,mu[j],sigma[j]);
          }
          f = prop%*%denm;
          run = run+1;
          lh[run] = sum(log(f));

          if(run > 2){


            if(lh[run-1]-lh[run-2] == 0){
              stop = 1;
            }else{
              c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);
              if(c>0 && c<1){
                preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);
                dif = ml-lh[run]; difml = abs(preml-ml);
              }
            }
          }

          if(dif<acc && difml<10^(-4) || run>50000){
            stop = 1;
          }else{
            p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
            mu = x%*%t(p)/apply(p,1,sum);
            prop = apply(p,1,sum)/n;
            sigma = sqrt(sum(sum(p*(t(matrix(rep(x,k),ncol=k))-matrix(rep(mu,n),ncol=n))^2))/n)*rep(1,k);
          }
        }

        main=list(mu=mu,pi=prop,sigma=sigma,inilh=l,lh=lh[run],p=p,iniid=id)

      }
    }
  }

  return(main)
}



#' Parameter Estimation for Multivariate Normal Mixture
#'
#' Estimate parameters of multivariate normal mixtures using EM-algorithm.
#'
#' @param x observations.
#' @param sigma_k Number of component. See details.
#' @param equal If the covariate matrix of each component is equal. See details.
#' @param ini initial value of parameters
#'
#'@details
#'equal:possible value: 0,1,.... When equal = 0, fit k component normal mixture with unknown covariance matrix; when equal = 1,
#'fit k component normal mixture with unknown equal covariance matrix; when equal = other value, fit k component normal mixture with covariance matrix is known.
#'
#'
#'sigma_k: Default is NULL. When both of sigma_k and equal is NULL then estimate the data as two component
#' multivariate normal with equal unknown covariates.
#'
#'
#' @return
#' mu: estimated mean
#' sigma: estimated sd
#' pi:estimated proportion
#' lh: likelihood of estimated parameters
#' p: density of observations
#' @export
#'
mixmnorm<-function(x,sigma_k = NULL,equal = NULL,ini = NULL){

  if(is.null(sigma_k) && is.null(equal)){
    k = 2; equal = 1;
  }
  if(is.null(sigma_k)==FALSE){
    temp = dim(as.matrix(sigma_k));
    if(temp[1]<temp[2]){
      sigma_k = t(sigma_k); temp = dim(sigma_k);
    }
    if(temp[1]>1){
      if(temp[1]==temp[2]){
        sigma = matrix(rep(sigma_k,2),nrow=2); k = 2;
        if(is.null(equal)==FALSE){
          k = equal;
        }
      }else{
        sigma = sigma_k; k = temp[1]/temp[2];
      }
      equal = 2;
    }else{
      k = sigma_k;
      if(is.null(equal)){
        equal = 1;
      }
    }
  }

  #if(dim(x)[1]<dim(x)[2]){
  if(NCOL(x)>NROW(x)) {
    #n = dim(x)[2]; s = dim(x)[1]; x = t(x);
    n = NCOL(x); s = NROW(x); x = t(x);
  }else{
    #n = dim(x)[1]; s = dim(x)[2];
    n = NROW(x); s = NCOL(x);
  }

  if(equal==0){
    #fit k component normal mixture with unknown covariance matrix
    mu = rbind(rsample(x,k),rsample(x,k),rsample(x,k),rsample(x,k),rsample(x,k));
    prop = matrix(rep(1/k,5*k),nrow=5);
    sigma = matrix(rep(0,k*s*5*s),ncol=s);
    for(i in 1:5){
      sigma[(k*s*(i-1)+1):(i*s*k),] = t(matrix(rep(diag(mmax(0.5,runif(s,0,1))*diag(var(x))),k),nrow=s));
    }

    l = rep(NaN,5);
    for(i in 1:5){
      stop = 0; run=0; l2 = -10^8;
      while(stop==0){
        run = run+1;
        denm = matrix(rep(0,n*k),nrow=k);
        for(j in 1:k){
          denm[j,] = dmvnorm(x,mu[(k*(i-1)+j),],sigma[(((i-1)*k+j-1)*s+1):(((i-1)*k+j)*s),]);
        }
        f = prop[i,]%*%denm;
        l1 = l2; l2 = sum(log(f)); dif = l2-l1;

        if(dif<10^(-5)*n && run>20 || dif<10^(-8)*n || run>10000){
          stop = 1;
        }else{
          p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
          mu[(k*(i-1)+1):(k*i),] = p%*%x/apply(p,1,sum);
          prop[i,] = apply(p,1,sum)/n;
          for(j in 1:k){
            sigma[(((i-1)*k+j-1)*s+1):(((i-1)*k+j)*s),] = t(x-t(matrix(rep(mu[k*(i-1)+j,],n),nrow=s)))%*%(matrix(rep(p[j,],s),ncol=s)*(x-t(matrix(rep(mu[k*(i-1)+j,],n),nrow=s))))/sum(p[j,]);
          }
        }
        if(min(prop[i,])<4/n){
          mu[(k*i-1+1):(k*i),] = rsample(x,k);  prop[i,] = rep(1/k,k); sigma[((i-1)*k*s+1):(i*k*s),] = t(matrix(rep(diag(mmax(0.05,runif(s,0,1))*diag(var(x))),k),nrow=s));
          stop = 0; run = 0;
        }
      }
      l[i] = l2;
    }
    I = sort(l,index.return = T)$ix; #y = sort(l,index.return = T)$x;
    id = I[5];
    mu = mu[(k*(id-1)+1):(k*id),]; prop = prop[id,]; sigma = sigma[(k*s*(id-1)+1):(k*s*id),];
    stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8; lh = numeric();

    while(stop == 0){

      denm = matrix(rep(0,n*k),nrow=k);
      for(j in 1:k){
        denm[j,] = dmvnorm(x,mu[j,],sigma[((j-1)*s+1):(j*s),]);
      }
      f = prop%*%denm;
      run = run+1;
      lh[run] = sum(log(f));#use lh to record the sequence of lh

      if(run > 3){

        if((lh[run-1]-lh[run-2])<10^(-8)){
          warning("Likelihood in EM algorithm is not increasing.")
        }
        if(lh[run-1]-lh[run-2] == 0){
          stop = 1;
        }else{
          c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);#convergence rate of Em algorithm
          if(c>0 && c<1){
            preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
            dif = ml-lh[run]; difml = abs(preml-ml);
          }
        }
      }

      #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
      if(dif<0.001 && difml<10^(-6) || run>100){
        stop = 1;
      }else{
        p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
        mu = p%*%x/apply(p,1,sum);
        prop = apply(p,1,sum)/n;
        for(j in 1:k){
          sigma[((j-1)*s+1):(j*s),] = t(x-t(matrix(rep(mu[j,],n),nrow=2)))%*%(matrix(rep(p[j,],s),ncol=s)*(x-t(matrix(rep(mu[j,],n),nrow=2))))/sum(p[j,]);
        }
      }
      if(min(prop)<2/n){
        est=mixmnorm(x,sigma_k,equal);
      }else{
        est=list(mu=mu,sigma=sigma,pi=prop,lh=lh[run],p=p)
      }
    }

    est=est

  }else if(equal==1){
    #unknown equal covariance matrix
    if(k==2){
      prop = 1/2*matrix(rep(1,4*2),nrow=4);
      b = sort(x[,1],index=T)$ix; c = sort(x[,2],index=T)$ix;
      sigma = matrix(rep(0,2*s*4),ncol=2);
      mu = rbind(rsample(x,k),rsample(x,k),x[b[1],],x[b[n],],x[c[1],],x[c[n],]);
      sigma[1:(2*s),] = rbind(diag(mmax(0.1,runif(s,0,1))*diag(var(x))),diag(mmax(0.1,runif(s,0,1))*diag(var(x))));
      sigma[(2*s+1):(4*s),] = rbind(diag(mmax(0.1,diag(var(x))*(n-1)/n-diag(var(mu[5:6,]))/2)),diag(mmax(0.1,diag(var(x))*(n-1)/n-diag(var(mu[7:8,]))/2)));

      l = rep(NaN,4);
      for(i in 1:4){
        stop = 0; run=0; l2 = -10^8;
        while(stop==0){
          run = run+1;
          denm = rbind(dmvnorm(x,mu[(2*i-1),],sigma[((i-1)*s+1):(i*s),]),dmvnorm(x,mu[(2*i),],sigma[((i-1)*s+1):(i*s),]));
          f = prop[i,]%*%denm;
          l1 = l2; l2 = sum(log(f)); dif = l2-l1;

          if(dif<10^(-5)*n && run>20 || dif<10^(-8)*n || run>1000){
            stop = 1;
          }else{
            p = matrix(rep(prop[i,],n),ncol=n)*denm/rbind(f,f);
            mu[(2*i-1):(2*i),] = p%*%x/apply(p,1,sum);
            prop[i,] = apply(p,1,sum)/n;
            sigma[((i-1)*s+1):(i*s),] = (t(x-t(matrix(rep(mu[2*i-1,],n),nrow=2)))%*%(matrix(rep(p[1,],s),ncol=s)*(x-t(matrix(rep(mu[2*i-1,],n),nrow=2))))+t(x-t(matrix(rep(mu[2*i,],n),nrow=2)))%*%(matrix(rep(p[2,],s),ncol=s)*(x-t(matrix(rep(mu[2*i,],n),nrow=2)))))/n;
          }
        }
        l[i] = l2;
      }

      I = sort(l,index.return = T)$ix; #y = sort(l,index.return = T)$x;
      id = I[4];
      mu = rbind(mu[(2*id-1),],mu[(2*id),]); prop = prop[id,]; sigma = sigma[(s*(id-1)+1):(s*id),];
      stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8; lh = numeric();

      while(stop == 0){

        denm = rbind(dmvnorm(x,mu[1,],sigma),dmvnorm(x,mu[2,],sigma));
        f = prop%*%denm;
        run = run+1;
        lh[run] = sum(log(f));#use lh to record the sequence of lh

        if(run > 2){

          if(lh[run-1]-lh[run-2] == 0){
            stop = 1;
          }else{
            c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);#convergence rate of Em algorithm
            if(c>0 && c<1){
              preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
              dif = ml-lh[run]; difml = abs(preml-ml);
            }
          }
        }

        #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
        if(dif<0.0001 && difml<10^(-7) || run>50000){
          stop = 1;
        }else{
          p = matrix(rep(prop,n),ncol=n)*denm/rbind(f,f);
          mu = p%*%x/rep(apply(p,1,sum),s);
          prop = apply(p,1,sum)/n;
          sigma = (t(x-t(matrix(rep(mu[1,],n),nrow=2)))%*%(matrix(rep(p[1,],s),ncol=s)*(x-t(matrix(rep(mu[1,],n),nrow=2))))+t(x-t(matrix(rep(mu[2,],n),nrow=2)))%*%(matrix(rep(p[2,],s),ncol=s)*(x-t(matrix(rep(mu[2,],n),nrow=2)))))/n;
        }
      }
      p = matrix(rep(prop,n),ncol=n)*denm/rbind(f,f);

      est=list(mu=mu,sigma=sigma,pi=prop,likelihood=lh[run],p=p)
    }else{
      mu = rbind(rsample(x,k),rsample(x,k),rsample(x,k),rsample(x,k),rsample(x,k));
      prop = matrix(rep(1/k,5*k),nrow=5);
      sigma = matrix(rep(0,s*5*2),ncol=2);
      for(i in 1:5){
        sigma[(s*(i-1)+1):(i*s),] = diag(mmax(0.1,runif(s,0,1))*diag(var(x)));
      }

      l = rep(NaN,5);
      for(i in 1:5){
        stop = 0; run=0; l2 = -10^8;
        while(stop==0){
          run = run+1;
          denm = matrix(rep(0,n*k),nrow=k);
          for(j in 1:k){
            denm[j,] = dmvnorm(x,mu[(k*(i-1)+j),],sigma[((i-1)*s+1):(i*s),]);
          }
          f = prop[i,]%*%denm;
          l1 = l2; l2 = sum(log(f)); dif = l2-l1;

          if(dif<10^(-5)*n && run>20 || dif<10^(-8)*n || run>10000){
            stop = 1;
          }else{
            p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
            mu[(k*(i-1)+1):(k*i),] = p%*%x/apply(p,1,sum);
            prop[i,] = apply(p,1,sum)/n;
            temp = 0;
            for(j in 1:k){
              temp = temp+t(x-t(matrix(rep(mu[k*(i-1)+j,],n),nrow=2)))%*%(matrix(rep(p[j,],s),ncol=s)*(x-t(matrix(rep(mu[k*(i-1)+j,],n),nrow=2))));
            }
            sigma[((i-1)*s+1):(i*s),] = temp/n;
          }
        }
        l[i] = l2;
      }
      I = sort(l,index.return = T)$ix; #y = sort(l,index.return = T)$x;
      id = I[5];
      mu = mu[(k*(id-1)+1):(k*id),]; prop = prop[id,]; sigma = sigma[(s*(id-1)+1):(s*id),];

      stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8; lh = numeric();
      while(stop == 0){

        denm = matrix(rep(0,n*k),nrow=k);
        for(j in 1:k){
          denm[j,] = dmvnorm(x,mu[j,],sigma);
        }
        f = prop%*%denm;
        run = run+1;
        lh[run] = sum(log(f));#use lh to record the sequence of lh

        if(run > 2){

          if(lh[run-1]-lh[run-2] == 0){
            stop = 1;
          }else{
            c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);#convergence rate of Em algorithm
            if(c>0 && c<1){
              preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
              dif = ml-lh[run]; difml = abs(preml-ml);
            }
          }
        }

        #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
        if(dif<0.0001 && difml<10^(-7) || run>50000){
          stop = 1;
        }else{
          p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
          mu = p%*%x/apply(p,1,sum);
          prop = apply(p,1,sum)/n;
          temp = 0;
          for(j in 1:k){
            temp = temp+t(x-t(matrix(rep(mu[j,],n),nrow=2)))%*%(matrix(rep(p[j,],s),ncol=s)*(x-t(matrix(rep(mu[j,],n),nrow=2))));
          }
          sigma = temp/n;
        }
      }
      p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
      est=list(mu=mu,sigma=sigma,pi=prop,lh=lh[run],p=p,inilh=l,iniid=id)
    }

  }else{
    #covariance matrix is known
    mu = rbind(rsample(x,k),rsample(x,k),rsample(x,k),rsample(x,k),rsample(x,k));
    prop = matrix(rep(1/k,5*k),nrow=5);

    l = rep(NaN,5);
    for(i in 1:5){
      stop = 0; run=0; l2 = -10^8;
      while(stop==0){
        run = run+1;
        denm = matrix(rep(0,n*k),nrow=k);
        for(j in 1:k){
          denm[j,] = dmvnorm(x,mu[(k*(i-1)+j),],sigma[((j-1)*s+1):(j*s),]);
        }
        f = prop[i,]%*%denm;
        l1 = l2; l2 = sum(log(f)); dif = l2-l1;

        if(dif<10^(-5)*n && run>20 || dif<10^(-8)*n || run>10000){
          stop = 1;
        }else{
          p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
          mu[(k*(i-1)+1):(k*i),] = p%*%x/rep(apply(p,1,sum),s);
          prop[i,] = apply(p,1,sum)/n;
        }
      }
      l[i] = l2;
    }
    I = sort(l,index.return = T)$ix; #y = sort(l,index.return = T)$x;
    id = I[5];
    mu = mu[(k*(id-1)+1):(k*id),]; prop = prop[id,];

    stop = 0; run = 0; dif = 10; difml = 10; ml = 10^8; lh = numeric();
    while(stop == 0){

      denm = matrix(rep(0,n*k),nrow=k);
      for(j in 1:k){
        denm[j,] = dmvnorm(x,mu[j,],sigma[((j-1)*s+1):(j*s),]);
      }
      f = prop%*%denm;
      run = run+1;
      lh[run] = sum(log(f));#use lh to record the sequence of lh

      if(run > 2){

        if(lh[run-1]-lh[run-2] == 0){
          stop = 1;
        }else{
          c = (lh[run]-lh[run-1])/(lh[run-1]-lh[run-2]);#convergence rate of Em algorithm
          if(c>0 && c<1){
            preml = ml; ml=lh[run-1]+(lh[run]-lh[run-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
            dif = ml-lh[run]; difml = abs(preml-ml);
          }
        }
      }

      #algorithm stopped when the direction derivative was smaller than acc and the difference of the predicted maximum likelihood value and the current likelihood value is small than 0.005
      if(dif<0.0001 && difml<10^(-7) || run>50000){
        stop = 1;
      }else{
        p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
        mu = p%*%x/rep(apply(p,1,sum),s);
        prop = apply(p,1,sum)/n;
      }
    }
    p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,k),ncol=k));
    est=list(mu=mu,pi=prop,lh=lh[run],p=p)
  }

  return(est)
}



permfull<-function(x){
  lx = length(x);
  fact = 1:lx;
  if(lx==1){
    out = x;
  }else{
    n = gamma(lx);
    out = matrix(rep(0,lx*n*lx),ncol=lx)
    for(i in 1:lx){
      temp = permfull(x[fact!=i]); temp = as.matrix(temp);
      out[((i-1)*n+1):(i*n),] = cbind(x[i]*matrix(rep(1,n),nrow=n),temp);
    }
  }
  return(out)
}



mmax<-function(a,b){
  n = length(b); c=numeric();
  for(i in 1:n){
    if(a>b[i]){
      c[i] = a
    }else{
      c[i] = b[i]
    }
  }
  return(c)
}



rsample<-function(x,n = NULL,replace = NULL){
  x = as.matrix(x);
  s = dim(x);
  if(s[1]==1){
    x = t(x);
  }
  s = dim(x);
  dim1 = s[1]; dim2 = s[2];

  if(is.null(replace)){
    replace = 0;
  }
  if(is.null(n)){
    n = dim1;
  }

  if(replace==0){
    if(n>dim1){
      stop("??????? ERROR:the number of samples to be chosen is bigger than that of x in the situation without replacement")
    }else{
      b = sample(dim1);
      a = x[b[1:n],];
    }
  }else{
    a=x[sample(1:dim1,size=n),];
  }

  if(s[2]==1){
    a = t(a);
  }

  return(a)
}






