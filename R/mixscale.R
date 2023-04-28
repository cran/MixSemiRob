#' Continuous scale mixture approach for parametric mixture model
#'
#'Estimation procedure for a class of semiparametric
#'mixture models that is a mixture of unknown location-shifted symmetric distributions. The
#'proposed method assumes that the nonparametric symmetric distribution falls in a rich
#'class of continuous normal scale mixture distributions.
#'
#' @param x univariate data vector
#' @param ini  ini:ini$mu: initial mu vector; ini$sigma: initial common variance;
#' ini$pi: initial proportion vector vector
#' @param maxiter the maximum iteration the algorithm run, default is 50
#'
#' @return
#'out$mu: mu vector; out$pi: proportion vector; out$bb: support of Q; out$weight: weight of Q corresponding to ini.bb;
#'out$lik: log-likelihood (Since the algorithm converges to different values you may compare this)
#' @export
#'
#'
#'@importFrom pracma lsqlincon
#'@importFrom quadprog solve.QP
#' @examples
#' n<- 10
#' mu <- c(-2.5,0); sd <- c(0.8,0.6); pi <- c(0.3,0.7)
#' set.seed(2023)
#' n1=rbinom(n,1,pi[1])
#' x=c(rnorm(sum(n1),mu[1],sd[1]),rnorm(n-sum(n1),mu[2],sd[2]))
#' ini=list(pi=pi,mu=mu,sigma=sd);
#' out=mixscale(x,ini, maxiter=1);
mixscale<-function(x,ini, maxiter = 50){
  #This program only works for two components
  #x: univariate data
  #ini:
  #ini$mu: initial mu vector
  #ini$sigma: initial common variance
  #ini$pi: initial proportion vector vector
  #This program will return "out" object containing
  #out$mu: mu vector
  #out$pi: proportion vector
  #out$bb: support of Q
  #out$weight: weight of Q corresponding to ini.bb
  #out$lik: log-likelihood (Since the algorithm converges to different values you may compare this)

  x=as.matrix(x);
  data=x;
  n=max(dim(x));
  stop=n/100;
  cons=1/n;

  a=ini$mu;
  bb=ini$sigma;
  bb=c(bb[1]/2,bb[1],2*bb[1],(max(x)-mean(x))^2);
  c=ini$pi;
  weight=rep(1/length(bb),length(bb));

  for(itr in 1:maxiter){
    zz=numeric()
    for(i in 1:n){
      zz[i]=c[1]*normal_mixture(data[i],a[1]*rep(1,length(bb)),bb+cons,weight)/(c[1]*normal_mixture(data[i],a[1]*rep(1,length(bb)),bb+cons,weight)+c[2]*normal_mixture(data[i],a[2]*rep(1,length(bb)),bb+cons,weight));
    }

    f1=find_mu(data,bb,weight,a[1],zz,cons);
    a[1]=f1;
    f2=find_mu(data,bb,weight,a[2],1-zz,cons);
    a[2]=f2;
    c=c(mean(zz),1-mean(zz));
    # data,beta=rbind(a,c),weight,variance=bb,cons
    f3=find_sigma(data,rbind(a,c),weight,bb,cons);
    new_sigma=f3$newx;
    ttt=f3$ttt;
    bb=c(bb,new_sigma);
    weight=c(weight,0*new_sigma);

    f4=update_weight(data,bb,weight,rbind(a,c),cons);
    bb=f4$bb;
    weight=f4$weight;

    lik=loglik_s(data,a,bb+cons,c,weight);

    if(max(ttt)<stop&&itr>1){break}

  }

  out=list(mu=a,pi=c,bb=bb,weight=weight,lik=lik,itr=itr,ttt=max(ttt))
  out
}



######
#univariate noraml mixture
normal_mixture<-function(x,loc,variance,weight){
  answer=0;
  for(i in 1:length(loc)){
    answer=answer+weight[i]*normal_pdf(x,loc[i],variance[i]);
  }
  answer
}


######
normal_mixture2<-function(x,variance,weight,beta){
  answer=0;
  mu=beta[1,];
  p=beta[2,];
  d=length(mu);
  for(i in 1:length(weight)){
    answer=answer+weight[i]*normal_mixture(x,mu,variance[i]*rep(1,d),p);
  }
  answer
}



######
#univariate normal density with mean mu and standard error sigma
normal_pdf<-function(x,mu,sigma){
  y=x-mu;
  pdf = dnorm(y,sd=sigma);
  pdf
}



######
# find_mu(data,bb,weight,mu=a[1],zz,cons)
find_mu<-function(data,bb,weight,mu,zz,cons){
  nn=max(dim(data));

  myfun2<-function(mu){
    g=0;
    for(i in 1:nn){
      gg=0;
      for(j in 1:length(weight)){
        n_pdf = normal_pdf(data[i],mu,bb[j]+cons)

        if (identical(n_pdf, numeric(0))) {
          gg=gg

        }else{

          gg=gg+weight[j]*n_pdf;
        }
      }
      g=g-zz[i]*log(gg);
    }
    g

  }

  #fmincon=optimize(f=myfun2,lower=-20,upper=20);
  fmincon=optimize(f=myfun2,lower=min(data),upper=max(data));
  out=fmincon$minimum;
  out
}


# find_sigma(data,rbind(a,c),weight,bb,cons);
# beta=rbind(a,c)
# variance = bb
######
find_sigma<-function(data,beta,weight,variance,cons){
  mu=beta[1,];
  p=beta[2,];
  d=length(mu);
  n=max(dim(data));

  grad<-function(pointer){
    answer=0;
    tmp=numeric();
    for(i in 1:n){
      tmp[i]=0;
      for(j in 1:length(weight)){
        tmp[i]=tmp[i]+weight[j]*normal_mixture(data[i],mu,(variance[j]+cons)*rep(1,d),p);
      }
      answer=answer+normal_mixture(data[i],mu,(pointer+cons)*rep(1,d),p)/tmp[i];
    }
    answer=n-answer;
    answer
  }

  tmp1=data;
  max_range=max((mean(tmp1)-max(tmp1))^2,(mean(tmp1)-min(tmp1))^2);
  min_range=max(cons,(min(abs(tmp1))/2)^2);
  # x=seq(from=min_range,to=max_range,by=max_range/100);
  x=seq(from=min_range,to=max_range,by=(max_range-min_range)/100)
  x1=numeric()
  for(i in 1:15){
    x1[i]=0.001*2^(i-6);
  }
  x=sort(c(x1,x,variance/2));x=matrix(x,nrow=1);

  newx=numeric()
  ran=numeric()

  br=numeric()
  for(ii in 1:length(x)){
    br[ii]=(gradient_sigma(data,x[ii]+1.0e-10,beta,weight,variance,cons)-gradient_sigma(data,x[ii]-1.0e-10,beta,weight,variance,cons))/2.0e-10;
    if(ii>1&&br[ii-1]>0&&br[ii]<0){
      newx=cbind(newx,(x[ii-1]+x[ii])/2);
      ran=rbind(ran,cbind(x[ii-1],x[ii]));
    }
  }

  if(br[length(x)]>0){
    newx=cbind(newx,max_range+1);
    ran=rbind(ran,cbind(x[dim(x)[1],2],500));
  }

  final_can=numeric()
  for(t in 1:dim(newx)[1]){
    # fmincon=donlp2(newx[t],grad,ran[t,2],ran[t,1],as.matrix(1),500,-Inf);
    #answer=fmincon$par;
    fmincon=optimize(f=grad,lower=ran[t,1],upper=min(ran[t,2],500));
    # fmincon=optimize(f=grad,c(ran[t,1],min(ran[t,2],500)),pointer=newx[t]);
    answer=fmincon$minimum;
    fval=fmincon$objective;
    if(fval<-0.0000001){
      final_can=cbind(final_can,answer);
    }
  }

  if(!is.null(final_can)){
    newx=final_can;
  }else{
    ttt=0;
    return(warning('final_can is null'))
  }
  if(is.null(newx)){
    ttt=0;
    return(warning('newx is null'))
  }

  newx=unique(newx*1000000000000)/1000000000000;


  newxx=numeric()
  ttt=numeric()
  for(jj in 1:length(newx)){
    g_sigma=gradient_sigma(data,newx[jj],beta,weight,variance,cons);
    ttt=cbind(ttt,g_sigma);
    if(g_sigma>0.0000001){
      newxx=cbind(newxx,newx[jj]);
    }
  }

  newx=newxx;
  out=list(newx=newx,ttt=ttt)
  out
}



#####
update_weight<-function(data,bb,weight,beta,cons){
  mu=beta[1,];
  p=beta[2,];
  d=length(mu);
  n=max(dim(data));

  m=length(bb);
  old_weight=weight;
  S=matrix(numeric(n*m),nrow=n);
  z=S;

  for(itr in 1:50){
    for(i in 1:n){
      for(j in 1:m){
        S[i,j]=normal_mixture(data[i],mu,(bb[j]+cons)*rep(1,d),p)/normal_mixture2(data[i],bb+cons,weight,beta);
      }
    }
    eps=1.0e-14;
    weight=lsqlincon(S,2*rep(1,n),-diag(rep(1,m)),rep(0,m),rep(1,m),1,eps*rep(1,m),(1-eps)*rep(1,m));
    #weight=t(weight);

    if(max(svd(weight-old_weight)$d)<0.0000000000001){break}
    old_weight=weight;
  }

  for(i in 1:n){
    for(j in 1:m){
      S[i,j]=normal_mixture(data[i],mu,(bb[j]+cons)*rep(1,d),p)/normal_mixture2(data[i],bb+cons,weight,beta);
      z[i,j]=weight[j]*S[i,j]/bb[j];
    }
  }

  W=diag(apply(z,1,sum));
  k=which(weight<=1.0e-14);
  if (!identical(k,integer(0))){
  bb=bb[-k];
  weight=weight[-k];}

  out=list(bb=bb,weight=weight,W=W)
  out

}


######
gradient_sigma<-function(data,phi,beta,weight,variance,cons){
  answer=0;
  mu=beta[1,];
  p=beta[2,];
  d=length(mu);
  n=max(dim(data));

  tmp=numeric()
  for(i in 1:n){
    tmp[i]=0;
    for(j in 1:length(weight)){
      tmp[i]=tmp[i]+weight[j]*normal_mixture(data[i],mu,(variance[j]+cons)*rep(1,d),p);
    }
    answer=answer+normal_mixture(data[i],mu,(phi+cons)*rep(1,d),p)/tmp[i];
  }

  answer=answer-n;
  answer
}



######
loglik_s<-function(data,phi,bb,c,weight){
  lik=0;
  for(i in 1:length(data)){
    tmp=0;
    for(j in 1:length(bb)){
      tmp=tmp+weight[j]*(c[1]*normal_pdf(data[i],phi[1],bb[j])+c[2]*normal_pdf(data[i],phi[2],bb[j]));
    }
    lik=lik+log(tmp);
  }
  lik
}











