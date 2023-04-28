
#' Fully iterative backfitting estimator (FIB) of single index mixture models
#'
#' @param x independent variable; n by p matrix, each row is an observation
#' @param y dependent variable vector
#' @param h bandwidth
#' @param est_b alpha transpose in the paper, similar to beta in regression model
#' @param ini initial values: list(est_p=proportion of each component,est_mu=mean;est_var= var of each component)
#' @param grid grid point for linear interpolating.
#' @param maxiter maximum iteration the algorithm run, default is 100.
#'
#'
#'@importFrom MASS ginv
#'
#' @return  est_p:estimated proportion of each component \cr
#' est_mu:estimated mean of each component\cr
#' est_var:estimated variance of each component \cr
#' est_b:estimated beta of the regressions\cr
#' iter: number of iteration before converge.
#' @export
#'
#' @examples
#'  xx=NBA[,c(1,2,4)];
#'  yy=NBA[,3];
#'  x=xx/t(matrix(rep(sqrt(diag(var(xx))),length(yy)),nrow=3));x=as.matrix(x);
#'  y=yy/sd(yy);y=as.vector(y)
#'  ini_bs=sinvreg(x,y)
#'  ini_b=ini_bs$direction[,1];
#'  h=0.3442;
#'  est_b=ini_b
#'  est=sim(x[1:50,],y[1:50],h,ini_b, maxiter = 1);
sim<-function(x,y,h=NULL,est_b=NULL,ini=NULL,grid=NULL,maxiter = 100){
  n=length(y);
  x=as.matrix(x)

  if(is.null(h)){
    h=cv_bw(x,y)$h_opt;
  }
  if(is.null(est_b)){
    ini_b=sinvreg(x,y)
    est_b=ini_b$direction[,1];
  }
  est_b=as.matrix(est_b);
  z=x%*%est_b;

  if(is.null(ini)){
    ini=mixlinreg(z,y,z);
  }
  est_p=ini$est_p*matrix(rep(1,n*length(ini$est_p)),nrow=n);
  est_var=ini$est_var*matrix(rep(1,n*length(ini$est_var)),nrow=n);
  est_mu=ini$est_mu;

  iter=1;
  fmin=999999;
  diff=1;

  while(iter<maxiter&&diff>10^(-3)){
    #step 1
    z=x%*%est_b;
    est1=mixnonpar(z,y,z,est_p,est_mu,est_var,h);
    est_p=est1$est_p;
    est_var=est1$est_var;
    est_mu=est1$est_mu;

    obj2=function(b){
      u=x%*%b;
      est_p_u=matrix(rep(0,2*length(u)),ncol=2);
      est_mu_u=matrix(rep(0,2*length(u)),ncol=2);
      est_var_u=matrix(rep(0,2*length(u)),ncol=2);
      for(i in 1:2){
        est_p_u[,i]=approx(z,est_p[,i],u)$y;
        est_mu_u[,i]=approx(z,est_mu[,i],u)$y;
        est_var_u[,i]=approx(z,est_var[,i],u)$y;
      }
      d=matrix(rep(0,2*length(u)),ncol=2);
      for(i in 1:2){
        d[,i]=dnorm(y,est_mu_u[,i],est_var_u[,i]^0.5);
      }
      f=-sum(log(d%*%t(est_p_u)));
      f
    }

    #step 2
    #fminsearch=nlm(obj2,est_b);
    fminsearch=optim(est_b,obj2);
    est_b=fminsearch$par;
    fval=fminsearch$value;

    diff=abs(fval-fmin);
    fmin=fval;
    iter=iter+1;
  }

  est_b=est_b/max(svd(est_b)$d);
  z=x%*%est_b;
  est2=mixnonpar(z,y,z,est_p,est_mu,est_var,h);
  est_p=est2$est_p;
  est_var=est2$est_var;
  est_mu=est2$est_mu;

  res_p=numeric();res_mu=numeric();res_var=numeric();
  if(is.null(grid)){
    res_p=est_p;
    res_mu=est_mu;
    res_var=est_var;
  }else{
    for(i in 1:2){
      res1=approx(z,est_p[,i],grid)$y;
      res2=approx(z,est_mu[,i],grid)$y;
      res3=approx(z,est_var[,i],grid)$y;
      res_p=cbind(res_p,res1);
      res_mu=cbind(res_mu,res2);
      res_var=cbind(res_var,res3);
    }
  }

  out=list(est_p=res_p,est_mu=res_mu,est_var=res_var,est_b=est_b,iter=iter)
  out
}




#out2=simonestep(x,y,h);
#' One step estimation of iterative backfitting estimator (FIB) of single index mixture models
#'
#' @param x independent variable; n by p matrix, each row is an observation
#' @param y dependent variable vector
#' @param h bandwidth
#' @param est_b alpha transpose in the paper, similar to beta in regression model
#' @param ini initial values: list(est_p=proportion of each component,est_mu=mean;est_var= var of each component)
#' @param grid grid point for linear interpolating.
#' @return  est_p:estimated proportion of each component \cr
#' est_mu:estimated mean of each component \cr
#' est_var:estimated variance of each component \cr
#' est_b:estimated beta of the regressions
#'
#' @importFrom MASS ginv
#'
#' @export
#'
#' @examples
#'  xx=NBA[,c(1,2,4)];
#'  yy=NBA[,3];
#'  x=xx/t(matrix(rep(sqrt(diag(var(xx))),length(yy)),nrow=3));x=as.matrix(x);
#'  y=yy/sd(yy);y=as.vector(y)
#'  ini_bs=sinvreg(x,y)
#'  ini_b=ini_bs$direction[,1];
#'  h=0.3442;
#'  est_b=ini_b
#'  #used a smaller sample for a quicker demonstration of the function
#'  set.seed(123)
#'  est_onestep=simonestep(x[1:50,],y[1:50],h,ini_b);
simonestep<-function(x,y,h,est_b=NULL,ini=NULL,grid=NULL){
  n=length(y);x=as.matrix(x);
  if(is.null(est_b)){
    ini_b=sinvreg(x,y)
    est_b=ini_b$direction[,1];
  }
  est_b=as.matrix(est_b);
  z=x%*%est_b;
  if(is.null(ini)){
    ini=mixlinreg(z,y,z);
  }
  est_p=ini$est_p*matrix(rep(1,n*length(ini$est_p)),nrow=n);
  est_var=ini$est_var*matrix(rep(1,n*length(ini$est_var)),nrow=n);
  est_mu=ini$est_mu;

  est=mixnonpar(z,y,z,est_p,est_mu,est_var,h);
  est_p=est$est_p;
  est_var=est$est_var;
  est_mu=est$est_mu;

  res_p=numeric();res_mu=numeric();res_var=numeric();
  if(is.null(grid)){
    res_p=est_p;
    res_mu=est_mu;
    res_var=est_var;
  }else{
    for(i in 1:2){
      res1=approx(z,est_p[,i],grid)$y;
      res2=approx(z,est_mu[,i],grid)$y;
      res3=approx(z,est_var[,i],grid)$y;
      res_p=cbind(res_p,res1);
      res_mu=cbind(res_mu,res2);
      res_var=cbind(res_var,res3);
    }
  }

  out=list(est_p=res_p,est_mu=res_mu,est_var=res_var,est_b=est_b)
  out
}




#dimension reduction based on sliced inverse regression
#' Dimension Reduction Based on Sliced Inverse Regression
#'
#'Use in the examples for sim, simonestep functions
#'
#' @param x independent variable matrix
#' @param y response variable vector
#' @param numslice number of slice, default is 10
#'
#' @return
#' direction:direction vector
#' eigenvalue: eigenvalues for reduced x
#' matx: reduced x
#' numrep:number of observations in each slice
#' @export
#'
#' @examples
#' #see example section of sim function.
sinvreg<-function(x,y,numslice=NULL){
  n=length(y);
  if(is.null(numslice)){
    numslice=10;
  }

  ind=sort(y,index.return = T)$ix;
  mx=apply(x,2,mean);
  numrep=round(n/numslice+10^-6);
  matx=0;

  xbar=numeric()
  for(i in 1:(numslice-1)){
    a=((i-1)*numrep+1);
    b=i*numrep;
    xbar1=apply(x[ind[a:b],],2,mean);
    xbar=rbind(xbar,xbar1);
    matx=matx+numrep*t(matrix(xbar[i,]-mx,nrow=1))%*%(xbar[i,]-mx);
  }

  xbar=rbind(xbar,apply(x[ind[((numslice-1)*numrep+1):n],],2,mean));
  matx=matx+(n-(numslice-1)*numrep)*t(matrix(xbar[numslice,]-mx,nrow=1))%*%(xbar[numslice,]-mx);
  matx=MASS::ginv(var(x)*(n-1))%*%matx;

  ss=eigen(matx);
  sv=ss$vec;
  sd=ss$val;
  temp=abs(sd);
  b=sort(-temp,index.return = T)$ix;
  sd=temp[b];
  sv=sv[,b];

  out=list(direction=sv,eigenvalue=sd,matx=matx,numrep=numrep)
  out
}




mixlinreg<-function(x,y,u){
  n=length(x);

  #v1=4-sin(2*pi*u);v2=1.5+cos(3*pi*u);
  U=cbind(rep(1,length(u)),u);
  #beta1=ginv(t(U)%*%U)%*%t(U)%*%v1;beta2=ginv(t(U)%*%U)%*%t(U)%*%v2;
  #est_beta=cbind(beta1,beta2);est_var=c(0.1,0.1);est_p=c(0.5,0.5);
  X=cbind(rep(1,n),x);

 # est=emest_linear(X,y,2,est_beta,est_var,est_p);
 # est_beta=est$est_beta;est_var=est$est_var;est_p=est$est_p;
  est=regmixEM(y,x)
  est_beta=est$beta;est_var=(est$sigma)^2;est_p=est$lambda;
  est_mu=X%*%est_beta;
  est_mu_x0=U%*%est_beta;

  out=list(est_p=est_p,est_mu=est_mu,est_var=est_var,est_mu_x0=est_mu_x0)
  out
}



mixnonpar<-function(x,y,u,est_p,est_mu,est_var,h){
  n=dim(est_mu)[1];
  numc=dim(est_mu)[2];
  m=length(u);
  r=matrix(numeric(n*numc),nrow=n);
  f=r;
  est_mu_u=matrix(numeric(m*numc),nrow=m);est_p_u=est_mu_u;est_var_u=est_mu_u;
  rsd=1;likelihood=1;iter=0;acc=10^(-3);
  repeat
  {
    iter=iter+1;
    for(i in 1:numc){
      f[,i]=dnorm(y,est_mu[,i],est_var[,i]^0.5);
    }
    for(i in 1:numc){
      r[,i]=est_p[,i]*f[,i]/apply(f*est_p,1,sum);
      W=exp(-(((t(matrix(rep(t(x),n),nrow=n))-matrix(rep(x,n),ncol=n))^2)/(2*h^2)))/(h*sqrt(2*pi));
      R=t(matrix(rep(r[,i],n),ncol=n));
      RY=t(matrix(rep(r[,i]*y,n),ncol=n));
      est_p[,i]=apply(W*R,1,sum)/apply(W,1,sum);
      est_mu[,i]=apply(W*RY,1,sum)/apply(W*R,1,sum);
      RE=(t(matrix(rep(y,n),ncol=n))-matrix(rep(est_mu[,i],n),ncol=n))^2;
      est_var[,i]=apply(W*RE*R,1,sum)/apply(W*R,1,sum);
    }
    rsd=likelihood-sum(log(apply(f*est_p,1,sum)));
    likelihood=sum(log(apply(f*est_p,1,sum)));
    if(abs(rsd)<acc|iter>200){break}
  }

  for(i in 1:2){
    est_p_u[,i]=approx(x,est_p[,i],u)$y;
    est_mu_u[,i]=approx(x,est_mu[,i],u)$y;
    est_var_u[,i]=approx(x,est_var[,i],u)$y;
  }

  out=list(est_p=est_p_u,est_mu=est_mu_u,est_var=est_var_u)
  out
}



#cv_se
#cv_one
cv_bw<-function(x,y){
  x=as.matrix(x);y=as.matrix(y)
  n=length(y);
  numc=2;
  h0=seq(from=0.06,to=0.16,length.out=10);
  h_l=length(h0);
  J=10;
  cv=matrix(rep(0,h_l*J),nrow=h_l);
  size=floor(n/J);
  f=matrix(numeric(size*numc),ncol=numc);r=f;
  true_b=as.matrix(c(1,1,1)/sqrt(3));

  for(hh in 1:h_l){

    for(j in 1:J){

      if(j==1){
        test=matrix(c(y[((j-1)*size+1):(j*size)],x[((j-1)*size+1):(j*size),]),nrow=size);#nrow=size
        train=matrix(c(y[(j*size+1):n],x[(j*size+1):n,]),nrow=n-size);
      }
      if(j!=1&&j!=10){
        test=matrix(c(y[((j-1)*size+1):(j*size)],x[((j-1)*size+1):(j*size),]),nrow=size);#nrow=size
        train=matrix(c(y[1:((j-1)*size)],y[(j*size+1):n],x[1:((j-1)*size),],x[(j*size+1):n,]),nrow=n-size);
      }
      if(j==10){
        test=matrix(c(y[((j-1)*size+1):(j*size)],x[((j-1)*size+1):(j*size),]),nrow=size);#nrow=size
        train=matrix(c(y[1:((j-1)*size)],x[1:((j-1)*size),]),nrow=n-size);
      }

      ini=sinvreg(train[,-1],train[,1]);
      if(sum(true_b*ini$direction[,1])<0){
        ini_b=-ini$direction[,1];
      }else{
        ini_b=ini$direction[,1];
      }
      ini_b=as.matrix(ini_b);

      est=sim(train[,-1],train[,1],h0[hh],ini_b);
      #est=simonestep(train[,-1],train[,1],h0[hh],ini_b);
      est_p=est$est_p;
      est_mu=est$est_mu;
      est_var=est$est_var;
      est_b=est$est_b;

      z0=test[,-1]%*%est_b;
      z=train[,-1]%*%est_b;

      res_p=numeric();res_mu=numeric();res_var=numeric();
      for(i in 1:numc){
        res1=approx(z,est_p[,i],z0)$y;
        res2=approx(z,est_mu[,i],z0)$y;
        res3=approx(z,est_var[,i],z0)$y;
        res_p=cbind(res_p,res1);
        res_mu=cbind(res_mu,res2);
        res_var=cbind(res_var,res3);
      }

      for(c in 1:numc){
        f[,c]=dnorm(test[,1],res_mu[,c],res_var[,c]^0.5);
        f[,c]=Re(f[,c]*is.finite(f[,c]));
      }
      for(c in 1:numc){
        r[,c]=res_p[,c]*f[,c]/apply(f*res_p,1,sum);
      }

      cal=(apply(r*res_mu,1,sum)-test[,1])^2;
      cv[hh,j]=mean(cal[!is.na(cal)]);
    }

  }
  cv_opt=apply(cv,1,sum);
  m=max(cv_opt);
  I=which.min(cv_opt);
  h_opt=h0[I];
  CV=list(h_opt=h_opt,cv_opt=cv_opt)
}






#emest_linear<-function(X,y,numc,est_beta,est_var,est_p){
#  n=length(y);
#  f=matrix(numeric(n*numc),nrow=n);r=f;
#  rsd=1;likelihood=1;iter=0;acc=10^(-3);
#  repeat
#  {
#    iter=iter+1;
#    for(i in 1:numc){
#      f[,i]=dnorm(y,X%*%est_beta[,i],est_var[i]^0.5*rep(1,n));
#    }
#
#    for(i in 1:numc){
#      est_p=as.numeric(est_p)
#      r[,i]=est_p[i]*f[,i]/apply(f*est_p,1,sum);#(f%*%est_p);apply(f*est_p,1,sum);
#      #s=r[,i];yy=y*s;xx=X*matrix(rep(s,2),nrow=n);b=coefficients(lm(yy~xx));est_beta[,i]=c(b[2],b[3]);
#      S=diag(r[,i]);est_beta[,i]=ginv(t(X)%*%S%*%X)%*%t(X)%*%S%*%y;
#      est_var[i]=sum(r[,i]*(y-X%*%est_beta[,i])^2)/sum(r[,i]);
#    }
#
#    est_p=cbind(mean(r[,1]),mean(r[,2]));
#
#    rsd=likelihood-sum(log(f%*%t(est_p)));
#    likelihood=sum(log(f%*%t(est_p)));
#    if(abs(rsd)<acc|iter>20){break}
#  }

#  out=list(est_p=est_p,est_beta=est_beta,est_var=est_var,likelihood=likelihood)
#  out
#}







