
#########################
##mixlaplace based on 20 initial values
##############################################
#' Robust Mixture Regression Based on Laplace Distribution
#'
#' @param y response variable vector
#' @param x explanatory variables matrix with rows for each observation
#' @param k number of components, default is 2
#' @param numini number of initial values, default is 20
#'
#' @return
#' theta: estimated parameters matrix, the columns are beta(intercept, slopes), sigma, proportion for each component.
#' likelihood: likelihood of the estimated parameters
#' run: number of iterations to converge
#' diflh: difference of likelihood of last iteration and second to last iteration
#' lh: likelihood for all initial values
#' @export
#'
#' @examples
#' library(fpc)
#' data('tonedata')
#' y=tonedata$tuned
#' x=tonedata$stretchratio
#' k=160;x[151:k]=0;y[151:k]=5
#'\donttest{est_lap=mixreg_Lap(x,y,2,numini=1);}
#'
mixreg_Lap <-function(x,y,k=2,numini=20)
{ n=length(y); x=matrix(x,nrow=n)
a=dim(x);p=a[2]+1; n1=2*p; lh= rep(0,numini);
bet= matrix(rep(0,k*p),nrow=k);sig=0;
for(j in seq(k))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/k,k);sig=sig/n1/k;sig[sig<0.1] <- 0.5;sig[sig>3] <- 1;
est=mixlap_one(x=x,y=y,bet,sig,pr,k);lh[1]=est$likelihood;
for(i in seq(numini-1))
{bet= matrix(rep(0,k*p),nrow=k);sig=0;
for(j in seq(k))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);sig[sig<0.1] <- 0.5;sig[sig>3] <- 1;}
pr=rep(1/k,k);sig=sig/n1/k;sig[sig<0.1] <- 0.5;sig[sig>3] <- 1;
temp=mixlap_one(x,y,bet,sig,pr,k);lh[i+1]=temp$likelihood;
if(lh[i+1]>est$likelihood& min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;}}
est$inilikelihoodseq=lh;
est}


####################################################
#Mixture regression fitting by Laplace distribution
####################################################
#mixlap_one estimates the mixture regression parameters based on ONE initial value
mixlap_one<-function(x,y,bet,sig,pr,m=2)
{##function left{
  run=0
  n=length(y)
  x=matrix(x,nrow=n)
  X=cbind(rep(1,n),x)
  if(length(sig)>1 ){ #the case when the variance is unequal
    r=matrix(rep(0,m*n),nrow=n)

    lh=0
    for(j in seq(m))
    {
      r[,j]=y-X%*%bet[j,]
      r[abs(r)<0.0001] <- 0.0001
      lh=lh+(pr[j]/(sqrt(2)*sig[j]))*exp(-(sqrt(2)*abs(r[,j]))/sig[j])
    }
    lh=sum(log(lh))
    #E-steps
    repeat
    { ##repeat left{
      prest=c(bet,sig,pr)
      run=run+1
      plh=lh
      pk=matrix(rep(0,m*n),nrow=n)
      delta=matrix(rep(0,m*n),nrow=n)
      for(j in seq(m))
      {
        pk[,j]=pr[j]*pmax(10^(-300),(exp(-(sqrt(2)*abs(r[,j]))/sig[j]))/sig[j])
        delta[,j]=sig[j]/(sqrt(2)*abs(r[,j]))
      }
      pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n)
      #M-step
      np=apply(pk,2,sum)
      pr=np/n
      lh=0
      for(j in seq(m))
      {
        w1=diag(pk[,j])
        w2=diag(delta[,j])
        w=w1%*%w2
        w_star=matrix(rep(0,n),nrow=n)
        for(i in seq(n))
        {
          w_star[i,]=w[i,i]
        }
        bet[j,]=ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y
        r[,j]= y-X%*%bet[j,]
        r[abs(r)<0.0001] <- 0.0001
        sig[j]=sqrt(2*t(w_star)%*%(r[,j]^2)/np[j])
        sig[sig<0.1] <- 0.5        ##truncate the extremely small sigma
        sig[sig>3] <- 1
        lh=lh+(pr[j]/(sqrt(2)*sig[j]))*exp(-(sqrt(2)*abs(r[,j]))/sig[j])
      }
      lh=sum(log(lh))
      dif=lh-plh
      if(dif<10^(-5)|run>500){break}}}
  else{  #the case when the variance is equal

    r=matrix(rep(0,m*n),nrow=n)
    lh=0
    for(j in seq(m))
    {
      r[,j]=y-X%*%bet[j,]
      r[abs(r)<0.0001] <- 0.0001
      lh=lh+(pr[j]/(sqrt(2)*sig))*exp(-(sqrt(2)*abs(r[,j]))/sig)
    }
    lh=sum(log(lh))
    #E-steps
    repeat
    { ##repeat left{
      prest=c(bet,sig,pr)
      run=run+1
      plh=lh
      pk=matrix(rep(0,m*n),nrow=n)
      delta=matrix(rep(0,m*n),nrow=n)
      for(j in seq(m))
      {
        pk[,j]=pr[j]* pmax(10^(-300),(exp(-(sqrt(2)*abs(r[,j]))/sig))/sig)
        delta[,j]=sig/(sqrt(2)*abs(r[,j]))
      }
      pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n)
      #M-step
      np=apply(pk,2,sum)
      pr=np/n
      w_star=matrix(rep(0,n*m),nrow=n)
      for(j in seq(m))
      {
        w1=diag(pk[,j])
        w2=diag(delta[,j])
        w=w1%*%w2
        for(i in seq(n))
        {
          w_star[i,j]=w[i,i]
        }
        bet[j,]=ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y
        r[,j]= y-X%*%bet[j,]
        r[abs(r)<0.0001] <- 0.0001
      }
      sig=sqrt(2*sum(w_star*(r^2))/n)
      sig[sig<0.1] <- 0.5        ##truncate the extremely small sigma
      sig[sig>3] <- 1             ##truncate the extremely large sigma
      lh=0
      for(j in seq(m))
      {
        lh=lh+(pr[j]/(sqrt(2)*sig))*exp(-(sqrt(2)*abs(r[,j]))/sig)
      }
      lh=sum(log(lh))
      dif=lh-plh
      if(dif<10^(-5)|run>100){break}}
    sig=sig*rep(1,m)}#the case when the variance is equal
  est=list(theta= matrix(c(bet,sig,pr),nrow=m),likelihood=lh,run=run,diflh=dif)
  est}

