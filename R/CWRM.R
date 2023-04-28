
####################################################
#Mixture regression fitting by CWRM
####################################################
#CWRM_one estimates the mixture regression parameters based on ONE initial value
CWRM_one<-function(x,y,alpha=0.1,bet,sig,pr,m=2)
{##function left{
           run=0
           n=length(y)
           n1=floor(n*alpha)
           x=matrix(x,nrow=n)
           X=cbind(rep(1,n),x)
           a=dim(x)[2]
           b=bet[-1,]
           b0=matrix(bet[1,],nrow=1)
     if(sig[1,]!=sig[2,1]){ #the case when the variance is unequal
             r=matrix(rep(0,m*n),nrow=n)
              lh=0
             for(j in seq(m))
             {
             r[,j]=y-x%*%b[,j]-b0[,j]
               r[abs(r)<0.0001] <- 0.0001
              lh=lh+pr[j,]*pmax(10^(-300),dnorm(r[,j],0,sig[1,]))
              }
              ind=rank(-lh)
              obj=sum(log(lh[ind<(n-n1)]))

      #E-steps
        repeat
       { ##repeat left{
              prest=c(bet,sig,pr)
              run=run+1
              pobj=obj
              r=matrix(rep(0,m*n),nrow=n)
              pk=matrix(rep(0,m*n),nrow=n)
              for(j in seq(m))
            {
              r[,j]=y-x%*%b[,j]-b0[,j]
              r[abs(r)<0.0001] <- 0.0001
              pk[,j]=pr[j,]*pmax(10^(-300),dnorm(r[,j],mean=0,sd=sig[1,]))
             }
              D=rank(apply(pk,1,sum))
              for (i in seq(n))
              {
               if (D[i]< n1){pk[i,]=0
               }else{
                 pk[i,]=pk[i,]/matrix(apply(pk,1,sum),ncol=1)[i,]}
              }
           tao_ij=pk

      #M-step
             np=apply(tao_ij,2,sum)
             #pr=np/(n-n1)
             pr=as.matrix(np/(n-n1),nrow=m)
             sig=matrix(rep(0,m),nrow=m)
             mu=matrix(rep(0,m*a),nrow=m)
             T=array(dim=c(a,a,m),NA)
             J=matrix(rep(1,n),nrow=n)
             lh=0
            for (j in seq(m))
             {
              w=diag(tao_ij[,j])
              mu[j,]=(t(x)%*%w%*%J)/(apply(tao_ij,2,sum)[j])
              T[,,j]=(t(x-mu[j,])%*%w%*%(x-mu[j,]))/(apply(tao_ij,2,sum)[j])
             b[,j]=((t(x)%*%w%*%y)/(apply(tao_ij,2,sum)[j])-(sum(w%*%y)/(apply(tao_ij,2,sum)[j]))*((t(x)%*%w%*%J)/(apply(tao_ij,2,sum)[j])))/(t(t(tao_ij[,j])%*%(x^2))/(apply(tao_ij,2,sum)[j])-((t(x)%*%w%*%J)/(apply(tao_ij,2,sum)[j]))^2)
             b0[,j]=(sum(w%*%y)/(apply(tao_ij,2,sum)[j]))-t(b[,j])%*%((t(x)%*%w%*%J)/(apply(tao_ij,2,sum)[j]))
             r[,j]= y-x%*%b[,j]-b0[,j]
             r[abs(r)<0.0001] <- 0.0001
           sig[j,]=sqrt(sum(w%*%(r[,j]^2))/(apply(tao_ij,2,sum)[j]))
            sig[sig<0.1] <- 0.5        ##truncate the extremely small sigma
            sig[sig>3] <- 1
             lh=lh+pr[j,]*pmax(10^(-300),dnorm(r[,j],0,sig[j,]))
           }
            ind=rank(-lh)
              obj=sum(log(lh[ind<(n-n1)]))

           #lh=sum(log(lh))
           dif=obj-pobj
           if(dif<10^(-5)|run>500){break}}}

    else{  #the case when the variance is equal
            r=matrix(rep(0,m*n),nrow=n)
              lh=0
             for(j in seq(m))
             {
             r[,j]=y-x%*%b[,j]-b0[,j]
               r[abs(r)<0.0001] <- 0.0001
              lh=lh+pr[j,]*pmax(10^(-300),dnorm(r[,j],0,sig[1,]))
              }
              ind=rank(-lh)
              obj=sum(log(lh[ind<(n-n1)]))

      #E-steps
        repeat
       { ##repeat left{
              prest=c(bet,sig,pr)
              run=run+1
              pobj=obj
              r=matrix(rep(0,m*n),nrow=n)
              pk=matrix(rep(0,m*n),nrow=n)
              for(j in seq(m))
            {
              r[,j]=y-x%*%b[,j]-b0[,j]
              r[abs(r)<0.0001] <- 0.0001
              pk[,j]=pr[j,]*pmax(10^(-300),dnorm(r[,j],mean=0,sd=sig[1,]))
             }
              D=rank(apply(pk,1,sum))
              for (i in seq(n))
              {
               if (D[i]< n1){pk[i,]=0
               }else{
                 pk[i,]=pk[i,]/matrix(apply(pk,1,sum),ncol=1)[i,]}
              }
           tao_ij=pk
      #M-step
             np=apply(tao_ij,2,sum)
             pr=as.matrix(np/(n-n1),nrow=m)
             r=matrix(rep(0,m*n),nrow=n)
              #sig=matrix(rep(0,m),nrow=m)
             mu=matrix(rep(0,m*a),nrow=m)
             T=array(dim=c(a,a,m),NA)
             J=matrix(rep(1,n),nrow=n)
              for (j in seq(m))
             {
              w=diag(tao_ij[,j])
              mu[j,]=(t(x)%*%w%*%J)/(apply(tao_ij,2,sum)[j])
              T[,,j]=(t(x-mu[j,])%*%w%*%(x-mu[j,]))/(apply(tao_ij,2,sum)[j])
             b[,j]=((t(x)%*%w%*%y)/(apply(tao_ij,2,sum)[j])-(sum(w%*%y)/(apply(tao_ij,2,sum)[j]))*((t(x)%*%w%*%J)/(apply(tao_ij,2,sum)[j])))/(t(t(tao_ij[,j])%*%(x^2))/(apply(tao_ij,2,sum)[j])-((t(x)%*%w%*%J)/(apply(tao_ij,2,sum)[j]))^2)
             b0[,j]=(sum(w%*%y)/(apply(tao_ij,2,sum)[j]))-t(b[,j])%*%((t(x)%*%w%*%J)/(apply(tao_ij,2,sum)[j]))
             r[,j]= y-x%*%b[,j]-b0[,j]
             r[abs(r)<0.0001] <- 0.0001
             }
            sig=sqrt(sum(w%*%(r^2))/(n-n1))  #the case when the variance is equal
            sig[sig<0.1] <- 0.5        ##truncate the extremely small sigma
            sig[sig>3] <- 1
            sig=matrix(sig*rep(1,m),nrow=m)
            lh=0
            for (j in seq(m))
             {
             lh=lh+pr[j,]*pmax(10^(-300),dnorm(r[,j],0,sig[1,]))
             }
            ind=rank(-lh)
              obj=sum(log(lh[ind<(n-n1)]))

           #lh=sum(log(lh))
        dif=obj-pobj
        if(dif<10^(-5)|run>100){break}}}
        beta=cbind(t(b0),t(b))
     est=list(theta= matrix(c(beta,sig,pr),nrow=m),likelihood=lh,run=run,diflh=dif)
     est}


  #est=CWRM_one(x=x,y=y,bet,sig,pr,m=2,alpha=0.1)

#########################
## CWRM based on 20 initial values
##############################################
#' Trimmed Cluster Weighted Restricted Model
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
#'\donttest{est_CWRM=mixreg_CWRM(x,y,2,numini=1);}
#'
mixreg_CWRM <-function(x,y,k=2,numini=20)
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






