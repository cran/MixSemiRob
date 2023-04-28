





#' Maximum Smoothed Likelihood Estimation for a Class of Semiparametric Pareto Mixture Densities
#'
#' @param y observations as vector.
#' @param n number of observations length(y)
#' @param numc number of components
#' @param est_alpha initial value for alpha, the scale parameter of pareto distribution
#' @param est_beta initial value for beta, the shape parameter of pareto distribution
#' @param est_p initial values for proportions, in vector
#' @param est_h density for nonparametric component.
#' @param h bandwidth h for kernel estimation
#' @param maxiter maximum iteration the algorithm run, default is 100
#'
#' @return
#' est_beta: estimated beta
#' est_p: estimated proportion
#' est_h: estimated density for nonparametric component
#' est_mu: estimated component means
#' likelihood: likelihood of observations
#' @export
#'
#' @examples
#'
#'data(ROE)
#'y=as.matrix(ROE[1:1000,]) #used a smaller sample for a quicker demonstration of the function
#'n=length(y);
#'h= 1.06*(n*(1-0.102))^(-1/5)* min(sd(y), quantile(y,0.25)/1.34);
#'est_alpha=0.065;
#'est_beta=2.911;
#'est_p=c(0.102,0.898);
#'est_mu=0.1172;
#'est_h = 3.33* dt(3.33*(y-est_mu),2);#est_h=dnorm(y-est_mu,0)
#'numc=2;
#'out =paretomix1(y,n,numc,est_alpha,est_beta,est_p,est_h,h,maxiter=1);
paretomix1<-function(y,n,numc,est_alpha,est_beta,est_p,est_h,h,maxiter = 100){

  r=matrix(rep(0,numc*n),nrow=n);
  f=r;

  rsd = 1;
  likelihood = 1;
  iter = 1;

  while(abs(rsd)>10^-5 && iter<maxiter){

    f[,1]=est_beta*est_alpha^est_beta*y^(-est_beta-1)*(y>est_alpha);
    f[,2]=est_h;
    f[is.nan(f)]=0;

    r[,1]=est_p[1]*f[,1]/(est_p[1]*f[,1]+est_p[2]*f[,2]);
    r[,2]=est_p[2]*f[,2]/(est_p[1]*f[,1]+est_p[2]*f[,2]);

    index1=which(est_h==max(est_h));
    est_mu=y[index1];
    new_id=y>est_alpha;

    y2 = y[new_id==1];
    r2 = r[new_id==1,1];
    est_beta = sum(r2)/sum((log(y2)-log(est_alpha))*r2);

    for(j in 1:n){
      t1=(y[j]-(y-est_mu))/h;
      t2=(-y[j]-(y-est_mu))/h;
      temp0=sum((2*n*h)^-1*(2*pi)^-0.5*exp(-0.5*t1^2)*r[,2])+sum((2*n*h)^-1*(2*pi)^-0.5*exp(-0.5*t2^2)*r[,2]);
      est_h[j]=(est_p[2])^-1*temp0;
    }

    est_p=apply(r,2,mean);

    rsd = likelihood-sum(log(apply(f*est_p,1,sum)));
    likelihood = sum(log(apply(f*est_p,1,sum)));
    iter = iter+1;
  }

  out=list(est_beta=est_beta,est_p=est_p,est_h=est_h,est_mu=est_mu,likelihood=likelihood)
  out
}





#Currently not working
# paretomix<-function(y,n,numc,est_alpha,est_beta,est_p,est_h,est_mu,h){
#
#   est1=paretomix1(y,n,numc,est_alpha,est_beta,est_p,est_h,h);
#   est_beta=est1$est_beta;
#   est_p=est1$est_p;
#   est_h=est1$est_h;
#
#   get_likelihood=function(x){
#     #f=paretomix2(y,n,numc,x,est_beta,est_p,est_h,est_mu,h)$likelihood;
#     #neg_likelihood =-f;
#     #neg_likelihood #yan 2022 nov 6
#     -paretomix2(y,n,numc,x,est_beta,est_p,est_h,est_mu,h)$likelihood;
#
#   }
#   est_alpha=optim(est_alpha,get_likelihood)$par
#   #est_alpha=nlm(get_likelihood,c(est_alpha,est_alpha))$estimate;
#
#   est=paretomix1(y,n,numc,est_alpha,est_beta,est_p,est_h,h);
#   est
# }
#
# paretomix2<-function(y,n,numc,est_alpha,est_beta,est_p,est_h,est_mu,h){
#
#   r=matrix(rep(0,numc*n),nrow=n);
#   f=r;
#
#   rsd = 1;
#   likelihood = 1;
#   iter = 1;
#
#   while(abs(rsd)>10^-5 && iter<100){
#
#     f[,1]=est_beta*est_alpha^est_beta*y^(-est_beta-1)*(y>est_alpha);
#     f[,2]=est_h;
#     f[is.nan(f)]=0;
#
#     r[,1]=est_p[1]*f[,1]/(est_p[1]*f[,1]+est_p[2]*f[,2]);
#     r[,2]=est_p[2]*f[,2]/(est_p[1]*f[,1]+est_p[2]*f[,2]);
#
#     #estimate beta
#     new_id=y>est_alpha;
#     y2 = y[new_id==1];
#     r2 = r[new_id==1,1];
#     est_beta = sum(r2)/sum((log(y2)-log(est_alpha))*r2);
#
#     #estmate h()
#     for(j in 1:n){
#       t1=(y[j]-(y-est_mu))/h;
#       t2=(-y[j]-(y-est_mu))/h;
#       temp0=sum((2*n*h)^-1*(2*pi)^-0.5*exp(-0.5*t1^2)*r[,2])+sum((2*n*h)^-1*(2*pi)^-0.5*exp(-0.5*t2^2)*r[,2]);
#       est_h[j]=(est_p[2])^-1*temp0;
#     }
#
#     #estimate mu
#     spl=function(x){
#       f=spline(y,-est_h,xout=x)$y;
#       f
#     }
#     est_mu=nlm(spl,est_mu)$estimate;
#
#     #estimate p
#     est_p=apply(r,2,mean);
#
#     rsd = likelihood-sum(log(apply(f*est_p,1,sum)));
#     likelihood = sum(log(apply(f*est_p,1,sum)));
#     iter = iter+1;
#   }
#
#   out=list(est_beta=est_beta,est_p=est_p,est_h=est_h,est_mu=est_mu,likelihood=likelihood)
#   out
# }
#
#
#
#
