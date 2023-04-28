
#' Profile Likelihood Method for Normal Mixture with Unequal Variance
#'
#'
#'pfmix finds the profile likelihood of the m-component normal mixture model using the EM algorithm
#'when assuming that the ratio of the smallest variance to the largest variance is k.
#'The function also provide a visualized way to find MLE of k by finding the maximum interior mode for likelihood.
#'To achieve this, one extra step of plot out likelihood vs k plot for different k is needed.
#'
#' @param x Observations. Can be vector or matrix.
#' @param k Ratio of the smallest variance to the largest variance.
#' @param m The number of components, equals to 2 by default.
#' @param numini The number of initial values used to find the constrained--The default value is 20
#'
#' @return a list of lh:likelihood; mu:estimated means of components;sigma: estimated sd of
#' components; pi: estimated proportion of components.
#' @export
#'
#' @examples
#'
#' #Example 1
#'
#'
#' n=100;
#' u=runif(n,0,1);
#' x1=(u<=0.3)*rnorm(n,0,0.5)+(u>0.3)*rnorm(n,1,1);
#' k1=0.4378;
#' est1=pfmix(x1,k1)
#'
#'# Example 2
#'
#'
#' x2=(u<=0.3)*rnorm(n,0,0.5)+(u>0.3)*rnorm(n,1.5,1);
#' k2=0.4879;
#' est2=pfmix(x2,k2)
pfmix<-function(x,k,m = NULL,numini = NULL){
  x = as.matrix(x);
  if(dim(x)[1]>dim(x)[2]){
    x = t(x);
  }
  if(is.null(m)){
    m = 2;
  }
  if(is.null(numini)){
    numini = 20;
  }

  n = length(x); sx = sort(x);
  vk = k; acc = 10^(-6)*n;
  pm = permfull(1:m); lpm = max(dim(pm)[1],dim(pm)[2]);
  likelihood = numeric();

  for(t in 1:length(vk)){
    k = vk[t];
    prop = 1/m*matrix(rep(1,lpm*numini*m),ncol=m);
    sig = matrix(rep(0,lpm*numini*m),ncol=m);
    mu = matrix(rep(0,numini*m),ncol=m);
    for(i in 1:(numini-1)){
      mu[i,] = rsample(x,m);
    }
    mu[numini,] = sx[round((1:m)/m*n+0.0001)];
    for(i in 1:(lpm-1)){
      mu = rbind(mu,mu[(1:numini),pm[i,]]);
    }
    sig[,(1:(m-1))] = matrix(rep(sqrt(0.5*var(as.numeric(x))),lpm*numini*(m-1)),ncol=m-1);
    sig[,m] = sig[,1]*k;

    lh = numeric();
    for(i in 1:(lpm*numini)){
      stop = 0; run = 0; l2 = -10^8;
      while(stop==0){
        run = run+1;
        denm = matrix(rep(0,n*m),nrow=m);
        for(j in 1:m){
          denm[j,] = dnorm(x,mu[i,j],sig[i,j]);
        }
        f = prop[i,]%*%denm;
        l1 = l2; l2 = sum(log(f)); dif = l2-l1;

        if(dif<10^(-6)*n && run>30 || dif<10^(-8)*n || run>1000){
          stop = 1;
        }else{
          p = matrix(rep(prop[i,],n),ncol=n)*denm/t(matrix(rep(f,m),ncol=m));
          nc = apply(p,1,sum)
          mu[i,] = x%*%t(p)/nc;
          prop[i,] = nc/n;
          s = numeric();
          for(j in 1:m){
            s[j] = (x-mu[i,j])^2%*%p[j,];
          }

          b = sort(s/nc,index.return = T)$ix;
          ind = c(1,1); q = 10^10;
          if(s[b[1]]/nc[b[1]]<k^2*s[b[m]]/nc[b[m]]){
            for(l in 1:(m-1)){
              minn = m; j = 1;
              sigtemp = numeric();

              while(j<min(minn,m-l+1)){
                temp = (sum(s[b[1:l]])+k^2*sum(s[b[(m-j+1):m]]))/(sum(nc[b[1:l]])+sum(nc[b[(m-j+1):m]]));
                if(temp<s[b[l+1]]/nc[b[l+1]] && temp>k^2*s[b[m-j]]/nc[b[m-j]]){
                  minn = min(minn,j);
                  sigtemp[b[1:l]] = matrix(rep(sqrt(temp),l),nrow=1);
                  if((l+1)<=(m-j)){
                    sigtemp[b[(l+1):(m-j)]] = sqrt(s[b[(l+1):(m-j)]]/nc[b[(l+1):(m-j)]]);
                  }
                  sigtemp[b[(m-j+1):m]] = matrix(rep(sqrt(temp)/k,j),nrow=1);
                  oldq = q; sigtemp = matrix(sigtemp,nrow=1);
                  q = log(sigtemp)%*%nc+s%*%t(1/(sigtemp^2))/2;
                  if(q<oldq){
                    ind=c(1,j);
                    break
                  }else{
                    j = j+1;
                  }
                }
              }
            }

            l = ind[1]; j = ind[2];
            temp=(sum(s[b[1:l]])+k^2*sum(s[b[(m-j+1):m]]))/(sum(nc[b[1:l]])+sum(nc[b[(m-j+1):m]]));
            sig[i,b[1:l]] = matrix(rep(sqrt(temp),l),nrow=1);
            if((l+1)<=(m-j)){
              sig[i,b[(l+1):(m-j)]] = sqrt(s[b[(l+1):(m-j)]]/nc[b[(l+1):(m-j)]]);
            }
            sig[i,b[(m-j+1):m]] = matrix(rep(sqrt(temp)/k,j),nrow=1);
          }else{
            ind = c(1,2); q = 10^10;
            for(l in 1:(m-1)){
              sigtemp = numeric();
              for(j in (l+1):m){
                temp = (s[b[l]]+k^2*s[b[j]])/(nc[b[l]]+nc[b[j]]);
                if(temp<s[b[1]]/nc[b[1]] && temp>k^2*s[b[m]]/nc[b[m]]){
                  sigtemp[b[1:m]] = sqrt(s[b[1:m]]/nc[b[1:m]]);
                  sigtemp[b[l]] = sqrt(temp);
                  sigtemp[b[j]] = sqrt(temp)/k;
                  oldq = q; sigtemp = matrix(sigtemp,nrow=1);
                  q = log(sigtemp)%*%nc+s%*%t(1/(sigtemp^2))/2;
                  if(q<oldq){
                    ind=c(l,j);
                  }
                }
              }
            }

            l = ind[1]; j = ind[2];
            temp = (s[b[l]]+k^2*s[b[j]])/(nc[b[l]]+nc[b[j]]);
            sig[i,b[1:m]] = sqrt(s[b[1:m]]/nc[b[1:m]]);
            sig[i,b[l]] = sqrt(temp);
            sig[i,b[j]] = sqrt(temp)/k;
          }
        }
      }

      for(j in 1:m){
        denm[j,] = dnorm(x,mu[i,j],sig[i,j]);
      }
      f = prop[i,]%*%denm;
      lh[i] = sum(log(f));
    }

    I = sort(lh,index.return = T)$ix; #y = sort(l,index.return = T)$x;
    id = I[lpm*numini];
    mu = mu[id,]; prop = prop[id,]; sig = sig[id,];
    stop = 0; numiter = 0; dif = 10; difml = 10; ml = 10^8; lh = numeric();

    while(stop==0){
      for(j in 1:m){
        denm[j,] = dnorm(x,mu[j],sig[j]);
      }
      f = prop%*%denm;
      numiter = numiter+1;
      lh[numiter] = sum(log(f));
      if(numiter>2){
        temp = lh[numiter-1]-lh[numiter-2];
        if(temp==0){
          dif = 0;
          break
        }else{
          c=(lh[numiter]-lh[numiter-1])/temp;
          if(c>0 && c<1){
            preml = ml; ml=lh[numiter-1]+(lh[numiter]-lh[numiter-1])/(1-c);#use Aitken acceleration to predict the maximum likelihood value
            dif = ml-lh[numiter]; difml = abs(preml-ml);
          }
        }
      }

      if(dif<0.001 && difml<10^(-8) || numiter>50000){
        stop = 1;
      }else{
        p = matrix(rep(prop,n),ncol=n)*denm/t(matrix(rep(f,m),ncol=m));
        nc = apply(p,1,sum);
        mu = x%*%t(p)/nc;
        prop = nc/n;
        s = numeric();
        for(j in 1:m){
          s[j] = (x-mu[j])^2%*%matrix(p[j,]);
        }

        b = sort(s/nc,index.return = T)$ix;
        ind = c(1,1); q = 10^10;
        if(s[1]/nc[1]<k^2*s[m]/nc[m]){
          for(l in 1:(m-1)){
            minn = m; j = 1;
            sigtemp = numeric();

            while(j<min(minn,m-l+1)){
              temp = (sum(s[b[1:l]])+k^2*sum(s[b[(m-j+1):m]]))/(sum(nc[b[1:l]])+sum(nc[b[(m-j+1):m]]));
              if(temp<s[b[l+1]]/nc[b[l+1]] && temp>k^2*s[b[m-j]]/nc[b[m-j]]){
                minn = min(minn,j);
                sigtemp[b[1:l]] = matrix(rep(sqrt(temp),l),nrow=1);
                if((l+1)<=(m-j)){
                  sigtemp[b[(l+1):(m-j)]] = sqrt(s[b[(l+1):(m-j)]]/nc[b[(l+1):(m-j)]]);
                }
                sigtemp[b[(m-j+1):m]] = matrix(rep(sqrt(temp)/k,j),nrow=1);
                oldq = q; sigtemp = matrix(sigtemp,nrow=1);
                q = log(sigtemp)%*%nc+s%*%t(1/(sigtemp^2))/2;
                if(q<oldq){
                  ind=c(l,j);
                  break
                }else{
                  j = j+1;
                }
              }
            }
          }

          l = ind[1]; j = ind[2];
          temp=(sum(s[b[1:l]])+k^2*sum(s[b[(m-j+1):m]]))/(sum(nc[b[1:l]])+sum(nc[b[(m-j+1):m]]));
          sig[b[1:l]] = matrix(rep(sqrt(temp),l),nrow=1);
          if((l+1)<=(m-j)){
            sig[b[(l+1):(m-j)]] = sqrt(s[b[(l+1):(m-j)]]/nc[b[(l+1):(m-j)]]);
          }
          sig[b[(m-j+1):m]] = matrix(rep(sqrt(temp)/k,j),nrow=1);
        }else{

          ind = c(1,2); q = 10^10;
          for(l in 1:(m-1)){
            sigtemp = numeric();
            for(j in (l+1):m){
              temp = (s[b[l]]+k^2*s[b[j]])/(nc[b[l]]+nc[b[j]]);
              if(temp<s[b[1]]/nc[b[1]] && temp>k^2*s[b[m]]/nc[b[m]]){
                sigtemp[b[1:m]] = sqrt(s[b[1:m]]/nc[b[1:m]]);
                sigtemp[b[l]] = sqrt(temp);
                sigtemp[b[j]] = sqrt(temp)/k;
                oldq = q; sigtemp = matrix(sigtemp,nrow=1);
                q = log(sigtemp)%*%nc+s%*%t(1/(sigtemp^2))/2;
                if(q<oldq){
                  ind=c(l,j);
                }
              }
            }
          }

          l = ind[1]; j = ind[2];
          temp = (s[b[l]]+k^2*s[b[j]])/(nc[b[l]]+nc[b[j]]);
          sig[b[1:m]] = sqrt(s[b[1:m]]/nc[b[1:m]]);
          sig[b[l]] = sqrt(temp);
          sig[b[j]] = sqrt(temp)/k;
        }
      }
    }

    for(j in 1:m){
      denm[j,] = dnorm(x,mu[j],sig[j]);
    }
    f = prop%*%denm;
    likelihood[t] = sum(log(f));
  }

  if(length(vk)==1){
    out=list(lh=likelihood,mu=mu,sigma=sig,pi=prop)
  }else{
    out=list(lh=likelihood)
  }

  return(out)
}



#finds all the permuted values of the vector x.
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

