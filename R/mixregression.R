


#######EM algorithm###################
#' Robust Mixture Regression with Thresholding-Embedded EM Algorithm for Penalized Estimation
#'
#' @param y response variable vector
#' @param x explanatory variables matrix with rows for each observation
#' @param ininum number of random initial values used if initial value is not specified, default is 20
#' @param pi0 initial value vector for proportion
#' @param beta0 initial value matrix for beta, each row for each component.
#' @param sig0 initial value vector for sigma
#' @param gamma0 initial value for gamma the mean shift vector
#' @param method which threshold to use. Options are "HARD" and "SOFT", default is "HARD"
#' @param cont cutoff for minimum and maximum component sigma ratio. default is 0.001
#' @param lambda tuning parameter in penalty term
#' @param m number of components of mixture model, default is 2
#' @param runnum1 maximum number of iteration for E step, default is 50
#' @param runnum2 maximum number of iteration for M step, default is 10
#' @param outer.conv converging criteria for E step, default is 10^-2
#' @param inner.conv converging criteria for M step, default is 10^-2
#'
#' @return
#' pi: estimated proportion
#' beta: estimated regression coefficient beta
#' sigma: estimated sigma for components
#' gamma: estimated mean shift vectors
#' pij: posterior probability
#' run: total iterations
#' @importFrom stats optim
#'
#' @export
#'
#' @examples
#' library(fpc)
#' data('tonedata')
#' y=tonedata$tuned
#' x=tonedata$stretchratio
#' k=160;x[151:k]=0;y[151:k]=5
#' \donttest{est_RM2=RM2_mixreg(y,x,lambda=1, ininum=1, runnum1=1, runnum2=1)}
RM2_mixreg<-function(y,x,ininum = 20,pi0=NULL,beta0=NULL,sig0=NULL,gamma0=NULL,
method= c("HARD","SOFT")[1],
cont=0.001,lambda,m=2,runnum1=50,runnum2=10,
outer.conv=1e-2,inner.conv=1e-2){
# RM2_mixreg<-function(y,X,pi0=pi_0,beta0=beta_0,sig0=sig_0,gamma0=gamma_0,
#                      method= c("HARD","SOFT")[1],
#                      cont=0.1,lambda=lambda.k,m=2,runnum1=50,runnum2=10,
#                      outer.conv=1e-2,inner.conv=1e-2){


  #Add intercept variable Xin Shen
  # if (!intercept){

      X=cbind(1,x)

  #   } else{
  #
  #     X=x
  # }

  n=length(y)
  p=ncol(X)
  y = matrix(y)
  #Add default value for pi, beta, sigma, gamma Xin Shen
  #If the initial value is not assigned, estimate by mix regression with Gaussian assumption.
  #If the output of mixlin is updated to be more user friendly, be sure to update here as well
  if (is.null(pi0)|is.null(beta0)|is.null(sig0)|is.null(gamma0)){
    ini_est = mixlin(x,y,m,20)
  }

  if (is.null(pi0)){

    pi0 = matrix(ini_est$theta[,p+2])

  }

  if (is.null(beta0)){

    beta0 = t(ini_est$theta[,1:p])

  }


  if (is.null(sig0)){

    sig0 = matrix(ini_est$theta[,p+1])

  }

  if (is.null(gamma0)){

    gamma0 = matrix(0,nrow=n)

  }

  f_plh <- switch(method,
                  "HARD" = hard_plh,
                  "SOFT" = soft_plh)
  f_ipod <- switch(method,
                   "HARD" = hard_ipod,
                   "SOFT" = soft_ipod)


  pi1=as.matrix(pi0)
  beta1=as.matrix(beta0)
  sig1=as.matrix(sig0)
  gamma1=gamma0


  out_lh=f_plh(y=y,X=X,pi=pi1,beta=beta1,sigma=sig1,gamma=gamma1,lambda=lambda)

  dif=1
  obj<-vector()
  run=0
  restart = 0
  ##EM loop
  while (dif > outer.conv & run<runnum1 & restart<10)
  {
    out_lh0=out_lh
    run=run+1
    oldgamma=gamma1
    oldpi=pi1
    oldbeta=beta1
    oldsig=sig1

    pk=matrix(rep(0,m*n),nrow=n)
    for(j in seq(m))
    {
      pk[,j]=oldpi[j,]*pmax(10^(-300),dnorm(y-oldgamma*oldsig[j,]-X%*%oldbeta[,j],mean=0,sd=oldsig[j,]))
    }

    p_ij=pk/matrix(rep(apply(pk,1,sum),m),nrow=n)
    n_p=apply(p_ij,2,sum)
    ##M steps loop
    pi2=oldpi
    beta2=oldbeta
    sig2=oldsig
    gamma2=oldgamma


    inner_lh=f_plh(y=y,X=X,pi=oldpi,beta=oldbeta,sigma=oldsig,gamma=oldgamma,lambda=lambda)

    run1=0
    dif1=1
    while (dif1 > inner.conv & run1<runnum2)
    {
      inner_lh0=inner_lh
      run1=run1+1
      prepi=pi2
      prebeta=beta2
      presig=sig2
      pregamma=gamma2
      pi=matrix(rep(0,m),nrow=m)
      beta=matrix(rep(0,p*m),nrow=p)
      sig=matrix(rep(0,m),nrow=m)

      ##update gamma first
      gamma=matrix(c(rep(0,n)),nrow=n)
      for(i in seq(n))
      {
        gamma[i]=f_ipod(theta=((p_ij[i,])%*%((y[i,]-(X%*%beta2)[i,])/(sig2))),lambda=lambda)
      }
      gamma2=as.matrix(gamma,nrow=n)


      for(k in seq(m))
      {
        pi[k]=n_p[k]/n
        wt=diag(p_ij[,k])
        beta[,k]=ginv(t(X)%*%wt%*%X)%*%t(X)%*%wt%*%(y-gamma2*presig[k,])
      }
      sig=optim(presig[1,],y=y,X=X,beta=beta,p_ij=p_ij,gamma=gamma2,fn=sigmaf,method="BFGS",hessian=FALSE,m=m)$par
      pi2=as.matrix(pi,nrow=m)
      beta2=as.matrix(beta,nrow=p)
      sig2=matrix(c(rep(sig,m)),nrow=m)

      #check if there is extreme value
      if (min(sig2)/max(sig2)<cont){

        warning('One of the sigma is going to 0, start again with another initial value.')

        #generate one new initial value from mixline
        ini_est = mixlin(x,y,m,1)


        gamma_new = matrix(0,nrow=n)
        pi_new=matrix(ini_est$theta[,p+2])
        beta_new=t(ini_est$theta[,1:p])
        sig_new=matrix(ini_est$theta[,p+1])
        out1 = list(pi=pi_new,beta=beta_new,sig=sig_new,gamma=gamma_new,run=0)
        out_lh0=0
        obj<-vector()
        restart = restart+1
        break
      }
      # sig2[sig2<cont] <- 0.5        ##truncate the extremely small sigma
      # sig2[sig2>3] <- 1

      #updated loglikelihood in the inner-loop
      inner_lh=f_plh(y=y,X=X,pi=pi2,beta=beta2,sigma=sig2,gamma=gamma2,lambda=lambda)
      dif1=inner_lh-inner_lh0
      #if(dif1 < -0.1 & sigma.flag==FALSE) message(paste("Warning: Inner loop decreases.  DIFF = ", round(dif1,4),"  Sig.truncated = ", sigma.flag, sep=" "))
      ####
      if (dif1> (-1e-10) & dif1<0) {dif1=0}
      out1=list(pi=pi2,beta=beta2,sig=sig2,gamma=gamma2,run=run1)
    }
    pi1=out1$pi
    beta1=out1$beta
    sig1=out1$sig
    gamma1=out1$gamma


    out_lh=f_plh(y=y,X=X,pi=pi1,beta=beta1,sigma=sig1,gamma=gamma1,lambda=lambda)

    obj[run] <- out_lh
    dif=out_lh-out_lh0
    if(dif < 0 & sum((beta1-beta2)^2) > 1e-2) warning(paste("Outer loop decreases.  DIFF=  ", dif, sep=" "))
    if (dif> (-1e-10) & dif<0) {dif=0}
  }

  #check if the algorithm converged
  cov=T
  if(run>=runnum1){
    warning('The algorithm didn\'t converge. stopped due to reaching the maximum iteration.')
    cov=F
  }
  if (restart>=10){
    warning('The algorithm didn\'t converge. stopped due to too many restarts.')
    cov=F
    }


  out=list(pi=pi1,beta=beta1,sigma=sig1,gamma=gamma1,pij=p_ij,run=run,cov=cov)#,diff=dif1,dif=dif) We seems do not need these value to return
}













#library(gregmisc)
# library(robustbase)

huberpsi<-function(t,k=1.345){ out=pmax(-k,pmin(k,t));out}
bisquare<-function(t,k=4.685){out=t*pmax(0,(1-(t/k)^2))^2;out}
biscalew<-function(t){ t[which(t==0)]=min(t[which(t!=0)])/10; out=pmin(1-(1-t^2/1.56^2)^3,1)/t^2;out}

################
#the EM algorithm to fit the mixture of linear regression
#########################
#mixlinone estimates the mixture regression parameters by MLE based on ONE initial value
mixlinone<-function(x,y,bet,sig,pr,m=2){
  run=0; n=length(y);
  X=cbind(rep(1,n),x);
  if(length(sig)>1 ){ #the case when the variance is unequal
    r=matrix(rep(0,m*n),nrow=n);pk=r;lh=0;
    for(j in seq(m))
    {r[,j]=y-X%*%bet[j,];lh=lh+pr[j]*dnorm(r[,j],0,sig[j]);}
    lh=sum(log(lh));
    #E-steps
    repeat
    {  prest=c(bet,sig,pr);run=run+1;plh=lh;
    for(j in seq(m))
    { pk[,j]=pr[j]*pmax(10^(-300),dnorm(r[,j],0,sig[j]))}
    pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
    #M-step
    np=apply(pk,2,sum);pr=np/n;lh=0;
    for(j in seq(m))
    {w=diag(pk[,j]);
    bet[j,]=ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y;
    r[,j]= y-X%*%bet[j,]; sig[j]=sqrt(t(pk[,j])%*%(r[,j]^2)/np[j]);
    lh=lh+pr[j]*dnorm(r[,j],0,sig[j]);}
    lh=sum(log(lh));dif=lh-plh;
    if(dif<10^(-5)|run>500){break}}}
  else{  #the case when the variance is equal
    r=matrix(rep(0,m*n),nrow=n);pk=r; lh=0
    for(j in seq(m))
    {r[,j]=y-X%*%bet[j,];lh=lh+pr[j]*dnorm(r[,j],0,sig);}
    lh=sum(log(lh));
    #E-steps
    repeat
    {  prest=c(bet,sig,pr);run=run+1;plh=lh;
    for(j in seq(m))
    {      pk[,j]=pr[j]* pmax(10^(-300),dnorm(r[,j],0,sig)) }
    pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
    #M-step
    np=apply(pk,2,sum);pr=np/n;
    for(j in seq(m))
    {   w=diag(pk[,j]);
    bet[j,]=ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y;
    r[,j]= y-X%*%bet[j,]; }
    sig=sqrt(sum(pk*(r^2))/n);lh=0;
    for(j in seq(m))
    {lh=lh+pr[j]*dnorm(r[,j],0,sig);}
    lh=sum(log(lh));
    dif=lh-plh;
    if(dif<10^(-5)|run>500){break}}
    sig=sig*rep(1,m)}
  est=list(theta= matrix(c(bet,sig,pr),nrow=m),likelihood=lh,run=run,diflh=dif)
  est}

##############################################
##mixlin based on multiple initial values
##############################################
#----------------------------------------------------------------------------------------------
# Programming Note:
# This function is used to estimate the initial values of parameters in RM2_mixreg function
# If any update is made about the output of this function, make sure to update corresponding part in RM2_mixreg.
#----------------------------------------------------------------------------------------------
#' Traditional MLE based mixture regression assuming the normal error density
#'
#' @param y response variable vector
#' @param x explanatory variables matrix with rows for each observation
#' @param k number of components, default is 2
#' @param numini number of initial values, default is 20.
#'
#' @return
#' theta: estimated parameters matrix, the columns are beta(intercept, slopes), sigma, proportion for each component.
#' likelihood: likelihood of the estimated parameters
#' run: number of iterations to converge
#' diflh: difference of likelihood of last iteration and second to last iteration
#'
#' @export
#'
#' @examples
#' library(fpc)
#' data('tonedata')
#' y=tonedata$tuned
#' x=tonedata$stretchratio
#' k=160;x[151:k]=0;y[151:k]=5
#'est=mixlin(x,y,2,numini = 1);
#'mle=est$theta
#'
mixlin <-function(x,y,k=2,numini=20)
{ n=length(y); x=matrix(x,nrow=n)
a=dim(x);p=a[2]+1; n1=2*p; lh= rep(0,numini);
bet= matrix(rep(0,k*p),nrow=k);sig=0;
for(j in seq(k))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/k,k);sig=sig/n1/k;
est=mixlinone(x,y,bet,sig,pr,k);lh[1]=est$likelihood;
for(i in seq(numini-1))
{bet= matrix(rep(0,k*p),nrow=k);sig=0;
for(j in seq(k))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/k,k);sig=sig/n1/k;
temp=mixlinone(x,y,bet,sig,pr,k);lh[i+1]=temp$likelihood;
if(lh[i+1]>est$likelihood& min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;}}
est$inilikelihoodseq=lh;
est}

################
#the robust EM algorithm to fit the mixture of linear regression based on bisquare function
#########################

#' Robust EM algorithm to fit the mixture of linear regression based on bisquare function with one initial value
#'
#' mixlinrb_bione estimates the mixture regression parameters robustly using bisquare function based on one initial value
#' @param x explanatory variables matrix with rows for each observation
#' @param y response variable vector
#' @param bet initial value of beta
#' @param sig initial value of sigma
#' @param pr initial value of proportion
#' @param m number of components, default is 2
#'
#' @return
#' theta: estimated parameters matrix, the columns are beta(intercept, slopes), sigma, proportion for each component.
#' run:number of iterations to converge
#'
#'
mixlinrb_bione<-function(x,y,bet,sig,pr,m=2){
  run=0;acc=10^(-4)*max(abs(c(bet,sig,pr))); n=length(y);
  X=cbind(rep(1,n),x);p=dim(X)[2];
  if(length(sig)>1)
  {
    r=matrix(rep(0,m*n),nrow=n);pk=r;
    for(j in seq(m))
      r[,j]=(y-X%*%bet[j,])/sig[j];
    #E-steps
    repeat
    {  prest=c(sig,bet,pr);run=run+1;
    for(j in seq(m))
    {
      pk[,j]=pr[j]*pmax(10^(-300),dnorm(r[,j],0,1))/sig[j]
    }
    pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
    #M-step
    np=apply(pk,2,sum);pr=np/n;
    r[which(r==0)]=min(r[which(r!=0)])/10;
    for(j in seq(m))
    {
      w=diag(pk[,j]*bisquare(r[,j])/r[,j]);
      bet[j,]= solve(t(X)%*%w%*%X+10^(-10)*diag(rep(1,p)))%*%t(X)%*%w%*%y;
      r[,j]= (y-X%*%bet[j,])/sig[j];
      sig[j]=sqrt(sum(r[,j]^2*sig[j]^2*pk[,j]*biscalew(r[,j]))/np[j]/0.5);
    }
    dif=max(abs(c(sig,bet,pr)-prest))
    if(dif<acc|run>500){break}
    }
  }
  else{ r=matrix(rep(0,m*n),nrow=n);pk=r;
  for(j in seq(m))
    r[,j]=( y-X%*%bet[j,])/sig;
  #E-steps
  repeat
  {  prest=c(sig,bet,pr);run=run+1;
  for(j in seq(m))
  {
    pk[,j]=pr[j]* pmax(10^(-300),dnorm(r[,j],0,1))/sig
  }
  pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
  #M-step
  np=apply(pk,2,sum);pr=np/n; r[which(r==0)]=min(r[which(r!=0)])/10;
  for(j in seq(m))
  { w=diag(pk[,j]*bisquare(r[,j])/r[,j]);
  bet[j,]=solve(t(X)%*%w%*%X+10^(-10)*diag(rep(1,p)))%*%t(X)%*%w%*%y;
  r[,j]=( y-X%*%bet[j,])/sig;
  }
  sig=sqrt(sum(pk*(r^2*sig[1]^2)*biscalew(r))/n/0.5)
  dif=max(abs(c(sig,bet,pr)-prest))
  if(dif<acc|run>500){break}
  }
  sig=rep(sig,m);
  }
  theta=matrix(c(bet,sig,pr),nrow=m);
  est=list(theta=theta,difpar=dif,run=run)
  # est=list(theta=theta,run=run)
  est
}

###

#' Robust EM algorithm to fit the mixture of linear regression based on bisquare function
#'
#'mixlinrb_bi estimates the mixture regression parameters robustly using bisquare function based on multiple initial values. The solution is found by the modal solution
#' @param x explanatory variables matrix with rows for each observation
#' @param y response variable vector
#' @param m number of components, default is 2
#' @param numini number of initial values, default is 20.
#'
#' @return
#'theta: estimated parameters matrix for the best estimator (the estimated value with most initial values converge to),
#' the columns are beta(intercept, slopes), sigma, proportion for each component
#'estall: all estimated parameters
#'uniqueest: unique estimated parameters
#'countuniqueest: number of unique estimated parameters
#'uniqueestindex: which initial values gives the unique estimators
#'bestindex: index of initial value gives the best estimated parameter
#'estindex: matrix of each initial value converge to which set of estimated parameters.
#' @export
#'
#'@examples
#' library(fpc)
#' data('tonedata')
#' y=tonedata$tuned
#' x=tonedata$stretchratio
#' k=160;x[151:k]=0;y[151:k]=5
#'est_bi=mixlinrb_bi(x,y,2,numini=2);
#'bi=est_bi$theta
mixlinrb_bi<-function(x,y,m=2, numini=20)
{ n=length(y); x=matrix(x,nrow=n);a=dim(x);p=a[2]+1; n1=2*p;
perm=permutations(m,m);sig=0;ind1=c(); bet=matrix(rep(0,2*p),nrow=2);
for(j in seq(m))
{ ind1= sample(1:n,n1); X=cbind(rep(1,n1),x[ind1,]);bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind1];
sig=sig+sum((y[ind1] -X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=max(sig/n1/m);
est=mixlinrb_bione(x,y,bet,sig,pr,m); lenpar=length(c(est$theta));
theta= matrix(rep(0,lenpar*(numini)),ncol=lenpar);
theta[1,]=c(est$theta); minsig=10^(-1)*sig;
trimbet=matrix(theta[1,1:(p*m)],nrow=m);
trimbet=matrix(rep(matrix(t(trimbet),ncol=p*m,byrow=T),gamma(m+1)),ncol=p*m,byrow=T);
ind=matrix(rep(0,numini),nrow=1);ind[1]=1;numsol=1;solindex=1; sol=matrix(theta[1,],nrow=1);

for(i in 2:numini)
{sig=0;ind1=c();
for(j in seq(m))
{ ind1= sample(1:n,n1); X=cbind(rep(1,n1),x[ind1,]);bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind1];
sig=sig+sum((y[ind1] -X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=max(sig/n1/m,minsig);
est=mixlinrb_bione(x,y,bet,sig,pr,m); theta[i,]=est$theta;
temp= matrix(theta[i,1:(p*m)],nrow=m);temp=matrix(t(temp[t(perm),]),ncol=p*m,byrow=T);
dif=apply((trimbet-temp )^2,1,sum);temp1=which(dif==min(dif));
theta[i,]=c(c(matrix(temp[temp1[1],],nrow=m,byrow=T)),theta[i,p*m+perm[temp1[1],]],theta[i,p*m+m+perm[temp1[1],]]);
dif= apply((matrix(rep(theta[i,1:(p*m)],numsol),nrow=numsol,byrow=T)-sol[,1:(p*m)])^2,1,sum);
if(min(dif)>0.1){sol=rbind(sol,theta[i,]); numsol=numsol+1; solindex=c(solindex,i);ind=rbind(ind,rep(0,numini));ind[numsol,i]=1}else{ind1=which(dif==min(dif));ind[ind1,i]=1;} }
num=apply(ind,1,sum); ind1=order(-num); bestindex=ind1;
for(j in seq(numsol)){ if(min(sol[ind1[j], (p*m+m+1):(p*m+2*m)])>0.05){index=1; est=matrix(sol[ind1[j],],nrow=m);
for(l in seq(m-1)){ temp=matrix(rep(est[l,1:p],m-l),nrow=m-l,byrow=T)-est[(l+1):m,1:p];
temp=matrix(temp,nrow=m-l); dif=apply(temp^2,1,sum);if(min(dif)<0.1){index=0;break} }
if(index==1){bestindex=ind1[j];break}}}
est= sol[bestindex[1],];
out=list(theta=matrix(est,nrow=m),estall=theta, uniqueest=sol, countuniqueest=num,uniqueestindex=solindex, bestindex= solindex[bestindex],estindex=ind);out }

##############################
#trimmed likelihood estimator
###################################
#trimmixone uses the trimmed likelihood estimator based on ONE initial value.
trimmixone<-function(x,y,k=2,alpha=0.95,bet,sig,pr){
  n=length(y);n1=round(n*alpha); x=matrix(x,nrow=n);a=dim(x);p=a[2]+1;
  X=cbind(rep(1,n),x);  if(dim(bet)[2]==k)  bet=t(bet);
  lh=0
  for (i in seq(k)) {
    lh=lh+pr[i]*dnorm(y-X%*%bet[i,],0,sig[1])}
  ind=order(-lh);run=0; acc=10^(-4);
  obj=sum(log(lh[ind[1:n1]]));
  repeat
  {pobj=obj;run=run+1;
  x1=x[ind[1:n1],];y1=y[ind[1:n1]];
  fit=mixlinone(x1,y1,bet,sig,pr,k);fit=fit$theta;
  bet=matrix(fit[1:(p*k)],nrow=k);sig=fit[p*k+1];pr=fit[(p*k+k+1):(p*k+2*k)];lh=0;
  for(i in seq(k))
  {lh=lh+pr[i]*dnorm(y-X%*%bet[i,],0,sig[1]);}
  ind=order(-lh);obj=sum(log(lh[ind[1:n1]]));dif=obj-pobj;
  if(dif<acc|run>50){break}
  }
  if(length(sig)<2)  sig=rep(sig,k);
  # theta=matrix(c(bet,sig,pr),nrow=m);est=list(theta=theta,likelihood=obj,diflikelihood=dif,run=run) #m is not working, updated with k by Xin Shen
  theta=matrix(c(bet,sig,pr),nrow=k);est=list(theta=theta,likelihood=obj,diflikelihood=dif,run=run)
  est}

#' Trimmed Likelihood Estimator
#'
#' @param x explanatory variables matrix with rows for each observation
#' @param y response variable vector
#' @param k number of components, default is 2
#' @param alpha proportion of data to be keep according to likelihood range from 0-1, default is 0.95
#' @param numini number of initial values, default is 20.
#'
#' @return
#' theta: estimated best parameters according to likelihood: beta for intercept ,beta for variables, sigma, proportion.
#' finallikelihood: likelihood of the best estimated parameters
#' likelihoodseq: all likelihood for all parameters estimated from initial values
#' @export
#'
#' @examples
#' library(fpc)
#' data('tonedata')
#' y=tonedata$tuned
#' x=tonedata$stretchratio
#' k=160;x[151:k]=0;y[151:k]=5
#' est_TLE=trimmix(x,y,2,0.95,numini=1);
#' TLE=est_TLE$theta
trimmix<-function(x,y,k=2,alpha=0.95,numini=20)
{ n=length(y); x=matrix(x,nrow=n)
a=dim(x);p=a[2]+1; n1=2*p;
bet= matrix(rep(0,k*p),nrow=k);sig=0;
for(j in seq(k))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);
}
pr=rep(1/k,k);sig=sig/n1/k; lh=rep(0,numini);
est=trimmixone(x,y,k,alpha,bet,sig,pr);lh[1]=est$likelihood;
for(i in seq(numini-1))
{ sig=0;
for(j in seq(k))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);
}
pr=rep(1/k,k);sig=sig/n1/k;
temp=trimmixone(x,y,k,alpha,bet,sig,pr);lh[i+1]=temp$likelihood;
if(lh[i+1]>est$likelihood){est=temp;}
}
est=list(theta=est$theta,finallikelihood=est$likelihood,likelihoodseq=lh)
est}

################
#the robust EM algorithm to fit the mixture of linear regression using t-distribution
#########################
#Definition of t density
dent<-function(y,mu,sig,v){
  est=gamma((v+1)/2)*sig^(-1)/((pi*v)^(1/2)*gamma(v/2)*(1+(y-mu)^2/(sig^2*v))^(0.5*(v+1)));
  est}

# mixlint estimates the mixture regression parameters robustly assuming the error distribution #is t-distribution
mixlintonev<-function(x,y,bet,sig,pr,v=3,m=2,acc=10^(-5)){
  run=0; n=length(y); if(length(v)==1){v=rep(v,m);}
  X=cbind(rep(1,n),x); lh=-10^10; a=dim(x);p=a[2]+1;
  r=matrix(rep(0,m*n),nrow=n);pk=r;u=r;logu=r;
  if(length(sig)>1)  #the component variance are different
  {
    for(j in seq(m))
      r[,j]=(y-X%*%bet[j,])/sig[j];
    #E-steps
    repeat
    { # prest=c(sig,bet,pr);
      run=run+1;prelh=lh;
      for(j in seq(m))
      {
        pk[,j]=pr[j]*pmax(10^(-300),dent(r[,j]*sig[j],0,sig[j],v[j]));
        u[,j]=(v[j]+1)/(v[j]+r[,j]^2)
        #logu[,j]=log(u[,j])+(digamma((v[j]+1)/2)-log((v[j]+1)/2))
      }
      lh= sum(log(apply(pk,1,sum)));pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
      dif=lh-prelh; if(dif<acc|run>500){break}
      #M-step
      np=apply(pk,2,sum);pr=np/n;
      for(j in seq(m))
      {   w=diag(pk[,j]*u[,j]);
      bet[j,]= ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y;
      sig[j]=sqrt(sum((y-X%*%bet[j,])^2*pk[,j]*u[,j])/sum(pk[,j]));
      r[,j]= (y-X%*%bet[j,])/sig[j];  }
    }
  }
  else{
    for(j in seq(m))
      r[,j]=( y-X%*%bet[j,])/sig;
    #E-steps
    repeat
    {  run=run+1; prelh=lh;
    for(j in seq(m))
    {
      pk[,j]=pr[j]*dent(y, X%*%bet[j,],sig,v[j]);
      u[,j]=(v[j]+1)/(v[j]+r[,j]^2)
      #logu[,j]=log(u[,j])+(digamma((v[j]+1)/2)-log((v[j]+1)/2))
    }
    lh= sum(log(apply(pk,1,sum)));
    pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
    dif=lh-prelh;
    if(dif<acc|run>500){break}
    #M-step
    np=apply(pk,2,sum);pr=np/n;
    for(j in seq(m))
    {  w=diag(pk[,j]*u[,j]);
    bet[j,]= ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y;
    r[,j]=y-X%*%bet[j,];
    }
    sig=sqrt(sum(pk*r^2*u)/sum(pk)) ;r=r/sig;
    }
    sig=sig*rep(1,m);
  }
  theta= matrix(c(bet,sig,pr),nrow=m);
  est=list(theta=theta,likelihood =lh,vdegree=v,dif=dif,run=run);est
}

#
#' Robust Mixture Regression Based on T-distribution
#' mixlintv adaptively estimates the mixture regression parameters robustly assuming the error distribution is t-distribution
#' @param x explanatory variables matrix with rows for each observation
#' @param y response variable vector
#' @param bet initial values for beta
#' @param sig initial value for sigma; if the component variance is different then use a vector, otherwise a scaler would be fine.
#' @param pr initial value of proportion for each component, vector
#' @param m number of component, default is 2
#' @param maxv maximum degree of freedom, default is 15
#' @param acc stopping criteria, default is 10^(-5)
#'
#' @return
#'theta: estimated parameters beta(intercept, slopes), sigma, proportion for each component
#'likelihood: likelihood of the estimated parameter
#'dif:  difference of likelihood of last iteration and second to last iteration
#'run: number of iteration to converge
#'vdegree: estimated degree of the t distribution
#'degreerange: all degree of freedom use to estimate the parameters
#'vlikelihoodseq: all likelihood associated with degrees of freedom
#'
mixlintv<-function(x,y,bet,sig,pr,m=2,maxv=15,acc=10^(-5))
{ est=mixlinone(x,y,bet,sig,pr,m);
fv=0; a=dim(x);p=a[2]+1; lh=rep(0,maxv);
lh[1]=est$likelihood;
for(v in 1:maxv)
{ temp=mixlintonev(x,y,bet,sig,pr,v,m,acc);lh[v+1]=temp$likelihood;
if(lh[v+1]>est$likelihood & min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;fv=v;}
}
est=list(theta=est$theta,likelihood =est$likelihood,dif=est$dif,run=est$run);
est$vdegree=fv;est$degreerange=c(1,maxv);est$vlikelihoodseq=lh;
est}

#
#'Robust Mixture Regression Based on T-distribution
#'
#'mixlint  adaptively estimates the mixture regression parameters robustly assuming the error distribution is t-distribution
#' @param x explanatory variables matrix with rows for each observation
#' @param y response variable vector
#' @param m number of component, default is 2
#' @param maxv maximum degree of freedom default is 30
#' @param numini number of initial value use to estimate the parameter, default is 20
#' @param acc stopping criteria, default is 10^(-5)
#'
#' @return
#'theta: estimated parameters beta(intercept, slopes), sigma, proportion for each component
#'likelihood: likelihood of the estimated parameter
#'dif:  difference of likelihood of last iteration and second to last iteration
#'run: number of iteration to converge
#'vdegree: estimated degree of the t distribution
#'degreerange: all degree of freedom use to estimate the parameters
#'vlikelihoodseq: all likelihood associated with degrees of freedom
#'lh: likelihood list for all initial values
#' @importFrom robustbase covMcd
#' @export
#'
#' @examples
#' library(fpc)
#' data('tonedata')
#' y=tonedata$tuned
#' x=tonedata$stretchratio
#' k=160;x[151:k]=0;y[151:k]=5
#'\donttest{est_t=mixlint(x,y,2,numini=20,acc=0.1);}
#'
mixlint<-function(x,y,m=2,maxv=30,numini=20,acc=10^(-5))
{ n=length(y); x=matrix(x,nrow=n)
a=dim(x);p=a[2]+1; n1=2*p;
bet= matrix(rep(0,m*p),nrow=m);sig=0;
for(j in seq(m))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=sig/n1/m; lh=rep(0,numini);
est=mixlintv(x,y,bet,sig,pr,m,maxv,acc);lh[1]=est$likelihood;
for(i in seq(numini-1))
{sig=0;
for(j in seq(m))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=sig/n1/m;
temp=mixlintv(x,y,bet,sig,pr,m,maxv,acc);lh[i+1]=temp$likelihood;
if(lh[i+1]>est$likelihood& min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;}}
est$degreerange=c(1,maxv);
if(maxv<30 | acc>10^(-5)){bet=est$theta[,1:p];sig=est$theta[,p+1];pr=est$theta[,p+2];
est=mixlintv(x,y,bet,sig[1],pr,m,30,10^(-5));}
est$inilikelihoodseq=lh;
est}

mixlintw<-function(x,y,m=2,maxv=10,numini=20,acc=10^(-5))
{ n=length(y); x=matrix(x,nrow=n)
w=covMcd(x);w=w$mcd.wt;
x=x[w==1,];y=y[w==1];
est=mixlint(x,y,m,maxv,numini,acc);
est
}

#%%%%%%%%%%%%%%%%%%%%
#mixregt with unknown degree of freedoms
#%%%%%%%%%%%%%%%%%%%%%%%%55
#%%start from one initial value
mixlintonev1<-function(x,y,bet,sig,pr,m=2,v=3,maxv=30,acc=10^(-5)){
  run=0; n=length(y); if(length(v)==1){v=rep(v,m);}
  X=cbind(rep(1,n),x); lh=-10^10; a=dim(x);p=a[2]+1;
  r=matrix(rep(0,m*n),nrow=n);pk=r;u=r;logu=r;
  if(length(sig)>1)  #the component variance are different
  {
    for(j in seq(m))
      r[,j]=(y-X%*%bet[j,])/sig[j];
    #E-steps
    repeat
    { # prest=c(sig,bet,pr);
      run=run+1;prelh=lh;
      for(j in seq(m))
      {
        pk[,j]=pr[j]*pmax(10^(-300),dent(r[,j]*sig[j],0,sig[j],v[j]));
        u[,j]=(v[j]+1)/(v[j]+r[,j]^2)
        #logu[,j]=log(u[,j])+(digamma((v[j]+1)/2)-log((v[j]+1)/2))
      }
      lh= sum(log(apply(pk,1,sum)));pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
      dif=lh-prelh; if(dif<acc|run>500){break}
      #M-step
      np=apply(pk,2,sum);pr=np/n;
      for(j in seq(m))
      {   w=diag(pk[,j]*u[,j]);
      bet[j,]= ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y;
      sig[j]=sqrt(sum((y-X%*%bet[j,])^2*pk[,j]*u[,j])/sum(pk[,j]));
      r[,j]= (y-X%*%bet[j,])/sig[j];  q=rep(0,maxv);
      for(k in seq(maxv))
      {q[k]=-log(gamma(0.5*k))+0.5*k*log(0.5*k)+0.5*k*(sum((log(u[,j])-u[,j])*pk[,j])/np[j]+digamma(0.5*(v[j]+1))-log(0.5*(v[j]+1)));
      }
      ordq=order(-q);v[j]=ordq[1];
      }
    }
  }
  else{
    for(j in seq(m))
      r[,j]=( y-X%*%bet[j,])/sig;
    #E-steps
    repeat
    {  run=run+1; prelh=lh;
    for(j in seq(m))
    {
      pk[,j]=pr[j]*dent(y, X%*%bet[j,],sig,v[j]);
      u[,j]=(v[j]+1)/(v[j]+r[,j]^2)
    }
    lh= sum(log(apply(pk,1,sum)));
    pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
    dif=lh-prelh;
    if(dif<acc|run>500){break}
    #M-step
    np=apply(pk,2,sum);pr=np/n;
    for(j in seq(m))
    {  w=diag(pk[,j]*u[,j]);
    bet[j,]= ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y;
    r[,j]=y-X%*%bet[j,]; q=rep(0,maxv);
    for(k in seq(maxv))
    {q[k]=-log(gamma(0.5*k))+0.5*k*log(0.5*k)+0.5*k*(sum((log(u[,j])-u[,j])*pk[,j])/np[j]+digamma(0.5*(v[j]+1))-log(0.5*(v[j]+1)));
    }
    ordq=order(-q);v[j]=ordq[1];
    }
    sig=sqrt(sum(pk*r^2*u)/sum(pk)) ;r=r/sig;
    }
    sig=sig*rep(1,m);
  }
  theta= matrix(c(bet,sig,pr),nrow=m);
  est=list(theta=theta,vdegree=v, likelihood =lh,degreerange=c(1,maxv),run=run,dif=dif);est
}

##Start from mulitple degrees
mixlintone<-function(x,y,bet,sig,pr,m=2,maxv=30,acc=10^(-5),numdeg=100){
  lh=rep(0,numdeg); u=runif(m,0,1);u=pmax(1,round(maxv*u));
  est=mixlintonev1(x,y,bet,sig,pr,m,u,maxv,acc);
  lh[1]=est$likelihood; a=dim(x);p=a[2]+1;
  for(v in seq(numdeg)){
    u=runif(m,0,1);u=pmax(1,round(maxv*u));
    #u=pmax(1,-sort(round(-10*u)));
    temp=mixlintonev1(x,y,bet,sig,pr,m,u,maxv,acc);
    lh[v+1]=temp$likelihood;
    if(lh[v+1]>est$likelihood& min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;}
  }
  if(max(est$vdegree)==maxv&maxv<30){
    est$vdegree[est$vdegree==maxv]=30
    temp=mixlintonev(x,y,bet,sig,pr,est$vdegree,m,acc);
    lh[v+2]=temp$likelihood;
    if(lh[v+2]>est$likelihood& min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;}
  }
  est$vlikelihoodseq=lh;
  est
}

###Start from multiple initial values
mixlintueqv<-function(x,y,m=2,maxv=30,numini=20,acc=10^(-5))
{ n=length(y); x=matrix(x,nrow=n)
a=dim(x);p=a[2]+1; n1=2*p;
bet= matrix(rep(0,m*p),nrow=m);sig=0;
for(j in seq(m))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=sig/n1/m; lh=rep(0,numini);
est=mixlintone(x,y,bet,sig,pr,m,maxv,acc);lh[1]=est$likelihood;
for(i in seq(numini-1))
{sig=0;
for(j in seq(m))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind]-X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=sig/n1/m;
temp=mixlintone(x,y,bet,sig,pr,m,maxv,acc);lh[i+1]=temp$likelihood;
if(lh[i+1]>est$likelihood& min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;}}
if(maxv<30 | acc>10^(-5)){bet=est$theta[,1:p];sig=est$theta[,p+1];pr=est$theta[,p+2];
est=mixlintone(x,y,bet,sig[1],pr,m,30,10^(-5));}
est$inilikelihoodseq=lh;
est}

mixlintwueqv<-function(x,y,m=2,maxv=30,numini=20,acc=10^(-5))
{ n=length(y); x=matrix(x,nrow=n)
w=covMcd(x);w=w$mcd.wt;
x=x[w==1,];y=y[w==1];
est=mixlintueqv(x,y,m,maxv,numini,acc);
est}

# mixlint  adaptively estimates the mixture regression parameters robustly assuming the error #distribution is t-distribution
mixlintv1<-function(x,y,m=2,v,numini=20,acc=10^(-5))
{ n=length(y); x=matrix(x,nrow=n)
a=dim(x);p=a[2]+1; n1=2*p;
bet= matrix(rep(0,m*p),nrow=m);sig=0;
for(j in seq(m))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=sig/n1/m; lh=rep(0,numini);
est=mixlintonev(x,y,bet,sig,pr,v,m,acc);lh[1]=est$likelihood;
for(i in seq(numini-1))
{sig=0;
for(j in seq(m))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=sig/n1/m;
temp=mixlintonev(x,y,bet,sig,pr,v,m,acc);lh[i+1]=temp$likelihood;
if(lh[i+1]>est$likelihood& min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;}}
if(acc>10^(-5)){bet=est$theta[,1:p];sig=est$theta[,p+1];pr=est$theta[,p+2];
est=mixlintonev(x,y,bet,sig[1],pr,v,m,10^(-5));}
est$inilikelihoodseq=lh;
est}

mixlintueqv1<-function(x,y,m=2,v,numini=20,acc=10^(-5))
{ n=length(y); x=matrix(x,nrow=n)
a=dim(x);p=a[2]+1; n1=2*p;
bet= matrix(rep(0,m*p),nrow=m);sig=0;
for(j in seq(m))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=sig/n1/m; lh=rep(0,numini);
est=mixlintonev1(x,y,bet,sig,pr,m,v,acc);lh[1]=est$likelihood;
for(i in seq(numini-1))
{sig=0;
for(j in seq(m))
{ind=sample(1:n,n1); X=cbind(rep(1,n1),x[ind,]);
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind]-X%*%bet[j,])^2);}
pr=rep(1/m,m);sig=sig/n1/m;
temp=mixlintonev1(x,y,bet,sig,pr,m,v,acc);lh[i+1]=temp$likelihood;
if(lh[i+1]>est$likelihood& min(temp$theta[,p+2])>min(0.05,est$theta[,p+2])){est=temp;}}
if(acc>10^(-5)){bet=est$theta[,1:p];sig=est$theta[,p+1];pr=est$theta[,p+2];
est=mixlintonev1(x,y,bet,sig[1],pr,m,v,10^(-5));}
est$inilikelihoodseq=lh;
est}



### Function permutations
permutations <-function (n, r, v = 1:n){
  if (r == 1)
    matrix(v, n, 1)
  else if (n == 1)
    matrix(v, 1, r)
  else {
    X <- NULL
    for (i in 1:n) X<-rbind(X,cbind(v[i],permutations(n-1, r-1, v[-i])))
    X}}




#######################################
hard_penalty<-function(gamma,lambda)
{
  ((lambda^2)/2)*(gamma!=0)
}

hard_ipod<-function(theta, lambda)
{
  theta*(abs(theta)>lambda)
}

hard_plh<-function(y,X,pi,beta,sigma,gamma,lambda,m=2)
{lh=0
for (r in seq(m))
{
  lh=lh+pi[r,]*pmax(10^(-300),dnorm(y-gamma*sigma[r,]-X%*%beta[,r],mean=0,sd=sigma[r,]))
}
##the objective function =loglikehood-hard penalty
sum(log(lh))-sum(hard_penalty(gamma=gamma,lambda=lambda))
}
#########################################





#########################################
soft_penalty<-function(gamma,lambda)
{
  abs(gamma)*lambda
}

soft_ipod<-function(theta, lambda)
{
  a = abs(theta)-lambda
  ifelse(a<0,0,sign(theta)*a)
}


soft_plh<-function(y,X,pi,beta,sigma,gamma,lambda,m=2)
{lh=0
for (r in seq(m))
{
  lh=lh+pi[r,]*pmax(10^(-300),dnorm(y-gamma*sigma[r,]-X%*%beta[,r],mean=0,sd=sigma[r,]))
}
##the objective function =loglikehood-hard penalty
sum(log(lh))-sum(soft_penalty(gamma=gamma,lambda=lambda))
}
###################################################


##objective function for sigma estimation
sigmaf<-function(sig,y,X,p_ij,beta,gamma,m)
{func=0
for (r in seq(m))
{
  func=func+0.5*sum(p_ij[,r]*log(sig^2))+(t(p_ij[,r])%*%(y-gamma*sig-X%*%beta[,r])^2)/(2*sig^2)
}
func
}



