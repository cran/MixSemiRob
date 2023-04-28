#' Root Selection Method for Univariate Finite Normal Mixture Models
#'
#'
#'
#' @param y observations vector
#' @param n number of observations: length(y)
#' @param k number of components
#' @param int number of initial values for EM algorithm.
#' @param ini initial values of parameter list(p=p,m=m,s=s); Default is estimated from y.
#' @param true true parameters. true = c(pr\[1\],mu\[1\],sd\[1\],pr\[2\],mu\[2\],sd\[2\],...)
#' if not specified, then will estimated from initical value using normal mixture model.
#' @return
#' estks: estimated KS statistic
#' estcvm: Cramér–von Mises statistic
#' estkui: Kuiper's statistic
#' estwatson: Watson statistic
#' estad: Anderson Darling statistic
#' cpks: 1:result is significant based on KS statistic; 0: o.w.
#' cpcvm:1:result is significant based on Cramér–von Mises statistic; 0: o.w.
#' cpkui:1:result is significant based on Kuiper's statistic; 0: o.w.
#' cpwatson:1:result is significant based on Watson statistic; 0: o.w.
#' cpad:1:result is significant based on  Anderson Darling statist; 0: o.w.
#'
#' @export
#'
#'
#'@examples
#' int <- 7; n <- 100; k=3
#' pr<- c(0.3,0.3,0.4);mu <- c(0,2,5);sd <- c(1,0.5,1)
#' y1 <- rnorm(pr[1]*n,mu[1],sd[1]); y2 <- rnorm(pr[2]*n,mu[2],sd[2]);
#' y3 <- rnorm(pr[3]*n,mu[3],sd[3])
#' y <- c(y1,y2,y3); rm(y1,y2,y3)
#' out=root_selection(y,n,k,int)
#'
#'@importFrom stats  approx coefficients cov dnorm fft integrate kmeans lm median nlm optimize pnorm qchisq quantile rbinom rnorm runif sd spline var
#'@importFrom mixtools normalmixEM
root_selection<-function(y,n,k,int,ini = NULL,true = NULL){

  e=seq(1:n);

  ### initial value
  if(is.null(ini)){
    p <- matrix(rep(1/k),int,k)
    m <- matrix(rnorm(k*int,mean(y),sd(y)),int,k)
    s <- matrix(rep(sd(y)),int,k)
  }else{
    p=ini$p;m=ini$m; s=ini$s;
  }


  #### true value
  if(is.null(true)){
    thetaah <- matrix(0,int,(3*k+1))
    cd <- c()
    for(j in 1:int)
    {
      res <- normalmixEM (y, lambda = p[j,], mu = m[j,], sigma = s[j,],epsilon = 1e-06, maxit = 5000, maxrestarts = 1000)
      rl<-res$lambda; rm<-res$mu; rs<-res$sigma
      est<-c(rl[order(rm)],rm[order(rm)],rs[order(rm)])
      estt = t(matrix(est, nrow = k))
      for(kk in 1:(3*k)){
        thetaah[j,kk]=estt[kk]
      }
      thetaah[j,3*k+1] <- res$loglik

      la1 <-0
      for(kk in 1:k){
        lakk=thetaah[j,3*kk-2]*pnorm(sort(y),thetaah[j,3*kk-1],thetaah[j,3*kk])
        la1<-la1+lakk; la1
      }
      cd[j] <- sum((la1-(2*e-1)/(2*n))^2)+(1/(12*n))
    }
    lab <- which.min(cd)
    thetaa <- c(thetaah[lab,1:(3*k)])

    true = thetaa
  }


  #theta
  ### Label Swiching
  nor <- function(x) {sqrt(sum(x^2))}

  theta <- matrix(0,int,3*k+1)
  for(i in 1:int)
  {
    res <- normalmixEM (y, lambda = p[i,], mu = m[i,], sigma = s[i,],epsilon = 1e-06, maxit = 5000)

    if(k==2){
      est <- matrix(0,2,6)
      est[1,] <- c(res$lambda[1],res$mu[1],res$sigma[1],res$lambda[2],res$mu[2],res$sigma[2])
      est[2,] <- c(res$lambda[2],res$mu[2],res$sigma[2],res$lambda[1],res$mu[1],res$sigma[1])
      distance <- c(nor(est[1,]-true),nor(est[2,]-true))
      label <- which.min(distance)
      for(kk in 1:6){
        theta[i,kk] <- est[label,kk]
      }
    }else if(k==3){
      est <- matrix(0,6,9)
      est[1,] <- c(res$lambda[1],res$mu[1],res$sigma[1],res$lambda[2],res$mu[2],res$sigma[2],res$lambda[3],res$mu[3],res$sigma[3])
      est[2,] <- c(res$lambda[1],res$mu[1],res$sigma[1],res$lambda[3],res$mu[3],res$sigma[3],res$lambda[2],res$mu[2],res$sigma[2])
      est[3,] <- c(res$lambda[2],res$mu[2],res$sigma[2],res$lambda[1],res$mu[1],res$sigma[1],res$lambda[3],res$mu[3],res$sigma[3])
      est[4,] <- c(res$lambda[2],res$mu[2],res$sigma[2],res$lambda[3],res$mu[3],res$sigma[3],res$lambda[1],res$mu[1],res$sigma[1])
      est[5,] <- c(res$lambda[3],res$mu[3],res$sigma[3],res$lambda[1],res$mu[1],res$sigma[1],res$lambda[2],res$mu[2],res$sigma[2])
      est[6,] <- c(res$lambda[3],res$mu[3],res$sigma[3],res$lambda[2],res$mu[2],res$sigma[2],res$lambda[1],res$mu[1],res$sigma[1])

      distance <- c(nor(est[1,]-true),nor(est[2,]-true),nor(est[3,]-true),nor(est[4,]-true),nor(est[5,]-true),nor(est[6,]-true))
      label <- which.min(distance)
      for(kk in 1:9){
        theta[i,kk]=est[label,kk]
      }
    }else if(k>=4){
      rl<-res$lambda; rm<-res$mu; rs<-res$sigma
      est<-c(rl[order(rm)],rm[order(rm)],rs[order(rm)])
      estt = t(matrix(est, nrow = k))
      for(kk in 1:(3*k)){
        theta[i,kk]=estt[kk]
      }
    }

    theta[i,3*k+1] <- res$loglik
  }


  #######################
  ##  Compare methods  ##
  #######################
  cp1 <- 0;cp2 <- 0;cp3 <- 0;cp4 <- 0;cp5 <- 0;
  pr=c(); mu=c(); sd=c()
  for(kk in 1:k){
    pr[kk]=true[3*kk-2]; mu[kk]=true[3*kk-1]; sd[kk]=true[3*kk]
  }

  x <- sort(y)
  f1=0
  for(kk in 1:k){
    f2=pr[kk]*pnorm(x,mu[kk],sd[kk]);
    f1=f1+f2; f1
  }
  logl <- sum(log(f1))

  f_emx <- t(matrix(e/n,n,int))
  f_emx1 <- t(matrix((e-1)/n,n,int))
  f_emz <- t(matrix((2*e-1)/(2*n),n,int))

  fhatx <- matrix(0,int,n)
  for (j in 1:int)
  {
    fhatx[j,] <- theta[j,1]*pnorm(x,theta[j,2],theta[j,3])+theta[j,4]*pnorm(x,theta[j,5],theta[j,6])+theta[j,7]*pnorm(x,theta[j,8],theta[j,9])
  }

  fhatw <- matrix(0,int,n)
  w <- sort(y,decreasing = TRUE)
  for (j in 1:int)
  {
    fhatw[j,] <- theta[j,1]*pnorm(w,theta[j,2],theta[j,3])+theta[j,4]*pnorm(w,theta[j,5],theta[j,6])+theta[j,7]*pnorm(w,theta[j,8],theta[j,9])
  }

  dplus <- apply(abs(f_emx-fhatx),1,max)
  dminus <- apply(abs(fhatx-f_emx1),1,max)


  ### method KS ###
  kol <- apply(cbind(dplus,dminus),1,max)
  met1 <- theta[which.min(kol),3*k+1]
  set1 <- theta[which.min(kol),1:(3*k)]
  if(2*(met1-logl) <= qchisq(0.95,3*k-1)){cp1 <- cp1+1}


  ### method CVM ###
  cvm <- rowSums((fhatx-f_emz)^2)+(1/(12*n))
  met2 <- theta[which.min(cvm),3*k+1]
  set2 <- theta[which.min(cvm),1:(3*k)]
  if(2*(met2-logl) <= qchisq(0.95,3*k-1)){cp2 <- cp2+1}


  ### method Kuiper ###
  kui <- dplus+dminus
  met3 <- theta[which.min(kui),3*k+1]
  set3 <- theta[which.min(kui),1:(3*k)]
  if(2*(met3-logl) <= qchisq(0.95,3*k-1)){cp3 <- cp3+1}


  ### method Watson ###
  watson <- cvm - n*((apply(fhatx,1,mean)-0.5)^2)
  met4 <- theta[which.min(watson),3*k+1]
  set4 <- theta[which.min(watson),1:(3*k)]
  if(2*(met4-logl) <= qchisq(0.95,3*k-1)){cp4 <- cp4+1}


  ### method AD ###
  ad <- -(rowSums (sweep((log(fhatx)+log(1-fhatw)),MARGIN=2,(2*e-1),'*')))/n-n
  met5 <- theta[which.min(ad),3*k+1]
  set5 <- theta[which.min(ad),1:(3*k)]
  if(2*(met5-logl) <= qchisq(0.95,3*k-1)){cp5 <- cp5+1}

  out=list(estks=set1,estcvm=set2,estkui=set3,estwatson=set4,estad=set5,cpks=cp1,cpcvm=cp2,cpkui=cp3,cpwatson=cp4,cpad=cp5)
}


