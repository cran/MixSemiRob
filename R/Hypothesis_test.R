



#' Hypothesis testing for finite mixture models
#'
#'Test how many component in the data based on KS statistic.
#'
#' @param x observations
#' @param n length of x
#' @param alp alpha level of the test
#' @param k test 2 to k components
#' @param ite time of bootstrap
#' @param init initial value of number of components
#'
#' @return
#' ks:KS statistic,
#' cvm: Cramér–von Mises statistic,
#' kui:Kuiper's statistic,
#' wat: Watson statistic,
#' ad:Anderson Darling statistic
#' result: test result based on KS statistic
#' 1: k component is significant according to KS statistic; 0:o.w.
#'
#' @export
#'
#'@importFrom mixtools rnormmix
#' @examples
#' n<- 100
#' mu <- c(-2.5,0); sd <- c(0.8,0.6); w <- c(0.3,0.7)
#' n1=rbinom(n,1,0.3)
#'  x=c(rnorm(sum(n1),mu[1],sd[1]),rnorm(n-sum(n1),mu[2],sd[2]))
#' \donttest{out=Hypothesis_test(x,n,alp = 0.10,k = 4,ite = 2,init = 5)}
Hypothesis_test<-function(x,n,alp = 0.10,k = 4,ite = 5,init = 5){

  # ite -- bootstrap;  # init -- initial value

  e = seq(1:n);
  e_cdf <- e/n;
  e_cdf1 <-(e-1)/n;
  e_cd <- (2*e-1)/(2*n);

  pro=rep(0,k);

  ### test 1 component ###

  cdf   <- pnorm(sort(x),mean(x),sd(x))
  DP <- max(abs(e_cdf-cdf))
  DM <- max(abs(cdf-e_cdf1))
  ks <- max(DP,DM)
  cvm <- sum((cdf-e_cd)^2)+(1/(12*n))
  kui <- DP+DM
  wat <- cvm - n*((mean(cdf)-0.5)^2)
  ad <- -n-((sum((2*(1:n)-1)*(log(cdf)+log(1-sort(cdf,decreasing = TRUE))))/n))
  real <- c(ks,cvm,kui,wat,ad)

  ks<- c();cvm<- c();
  kui<-c();wat<-c();
  ad<-c();

  ### Bootstrap ###
  for (i in 1:ite)
  {
    xboot <- rnorm(n, mean(x),sd(x))
    cdf <- pnorm(sort(xboot),mean(xboot),sd(xboot))

    DP <- max(abs(e_cdf-cdf))
    DM <- max(abs(cdf-e_cdf1))
    ks[i] <- max(DP,DM)
    cvm[i] <- sum((cdf-e_cd)^2)+(1/(12*n))
    kui[i] <- DP+DM
    wat[i] <- cvm[i] - n*((mean(cdf)-0.5)^2)
    ad[i] <- -n-((sum((2*(1:n)-1)*(log(cdf)+
                                     log(1-sort(cdf,decreasing = TRUE))))/n))
  }
  if(real[1]<sort(ks)[(1-alp)*ite])
  {
    pro[1] <- pro[1]+1
  }
  if(real[1]>=sort(ks)[(1-alp)*ite])
  {
    ### test k component ###
    for(kk in 2:k){
      out=test_k_comp(x,kk,e,ite,init);

      if(out$real[1]<sort(out$ks)[(1-alp)*ite])
      {
        pro[kk] <- pro[kk]+1;
        ks=out$ks;cvm=out$cvm;kui=out$kui;wat=out$wat;ad=out$ad;
      }
      if(out$real[1]<sort(out$ks)[(1-alp)*ite]) {break}
    }
  }


  res=list(ks=ks,cvm=cvm,kui=kui,wat=wat,ad=ad,result=pro)

}




test_k_comp <- function(x,k,e,ite,init){

  ##########################
  #### test k component ####
  ##########################
  n=length(x) #added sep18 2022 Xin Shen to fix the error message of 'no visible binding for global variable ‘n’'
  nor <- function(x) {sqrt(sum(x^2))}

  e_cdf<- e/n;
  e_cdf1<-(e-1)/n;
  e_cd<- (2*e-1)/(2*n);

  m <- matrix(rnorm(k*init,mean(x),sd(x)),init,k)
  s <- matrix(rep(sd(x)),init,k)
  p <- matrix(rep(1/k),init,k)

  #### EM ######
  thetaah <- matrix(0,init,(3*k+1))
  cd <- c()
  for(j in 1:init)
  {
    res <- normalmixEM (x, lambda = p[j,], mu = m[j,],
                        sigma = s[j,],epsilon = 1e-06,
                        maxit = 5000, maxrestarts = 1000)

    ### Label Swiching ###
    rl<-res$lambda; rm<-res$mu; rs<-res$sigma
    est<-c(rl[order(rm)],rm[order(rm)],rs[order(rm)])
    estt = t(matrix(est, nrow = k))

    for(kk in 1:(3*k)){
      thetaah[j,kk]=estt[kk]
    }
    thetaah[j,3*k+1] <- res$loglik

    la1 <-0
    for(kk in 1:k){
      lakk=thetaah[j,3*kk-2]*pnorm(sort(x),thetaah[j,3*kk-1],thetaah[j,3*kk])
      la1<-la1+lakk
      la1
    }

    cd[j] <- sum((la1-e_cd)^2)+(1/(12*n))
  }

  lab <- which.min(cd)

  theta <- c(thetaah[lab,1:(3*k)])

  cdf <-0
  for(kk in 1:k){
    cdfkk=theta[3*kk-2]*pnorm(sort(x),theta[3*kk-1],theta[3*kk])
    cdf<-cdf+cdfkk
    cdf
  }


  DP <- max(abs(e_cdf-cdf))
  DM <- max(abs(cdf-e_cdf1))
  ks <- max(DP,DM)
  cvm <- sum((cdf-e_cd)^2)+(1/(12*n))
  kui <- DP+DM
  wat <- cvm - n*((mean(cdf)-0.5)^2)
  ad <- -n-((sum((2*(1:n)-1)*(log(cdf)+log(1-sort(cdf,decreasing = TRUE))))/n))
  real <- c(ks,cvm,kui,wat,ad)

  ks<- c(); cvm<- c(); kui<-c(); wat<-c(); ad<-c();
  true <- theta

  for (i in 1:ite)
  {
    lambda=numeric(); mu=numeric();sigma=numeric();
    for(kk in 1:k){
      lambda[kk]=true[3*kk-2]
      mu[kk]=true[3*kk-1]
      sigma[kk]=true[3*kk]
    }

    xb<- rnormmix(n, lambda, mu, sigma)

    m <- matrix(rnorm(k*init,mean(xb),sd(xb)),init,k)
    s <- matrix(rep(sd(xb)),init,k)
    p <- matrix(rep(1/k),init,k)

    #### EM ######
    thetaa <- matrix(0,init,3*k+1)
    cd <- c()
    for(j in 1:init)
    {
      res <- normalmixEM (xb, lambda = p[j,], mu = m[j,],
                          sigma = s[j,],epsilon = 1e-06,
                          maxit = 5000, maxrestarts = 1000)

      ### Label Swiching ###
      if(k==2){
        est <- matrix(0,2,6)
        est[1,] <- c(res$lambda[1],res$mu[1],res$sigma[1],
                     res$lambda[2],res$mu[2],res$sigma[2])
        est[2,] <- c(res$lambda[2],res$mu[2],res$sigma[2],
                     res$lambda[1],res$mu[1],res$sigma[1])

        distance <- c(nor(est[1,]-true),nor(est[2,]-true))
        label <- which.min(distance)

        for(kk in 1:(3*k)){
          thetaa[j,kk] <- est[label,kk]
        }
      }

      if(k==3){
        est <- matrix(0,6,9)
        est[1,] <- c(res$lambda[1],res$mu[1],res$sigma[1],
                     res$lambda[2],res$mu[2],res$sigma[2],
                     res$lambda[3],res$mu[3],res$sigma[3])
        est[2,] <- c(res$lambda[1],res$mu[1],res$sigma[1],
                     res$lambda[3],res$mu[3],res$sigma[3],
                     res$lambda[2],res$mu[2],res$sigma[2])
        est[3,] <- c(res$lambda[2],res$mu[2],res$sigma[2],
                     res$lambda[1],res$mu[1],res$sigma[1],
                     res$lambda[3],res$mu[3],res$sigma[3])
        est[4,] <- c(res$lambda[2],res$mu[2],res$sigma[2],
                     res$lambda[3],res$mu[3],res$sigma[3],
                     res$lambda[1],res$mu[1],res$sigma[1])
        est[5,] <- c(res$lambda[3],res$mu[3],res$sigma[3],
                     res$lambda[1],res$mu[1],res$sigma[1],
                     res$lambda[2],res$mu[2],res$sigma[2])
        est[6,] <- c(res$lambda[3],res$mu[3],res$sigma[3],
                     res$lambda[2],res$mu[2],res$sigma[2],
                     res$lambda[1],res$mu[1],res$sigma[1])


        distance <- c(nor(est[1,]-true),nor(est[2,]-true),
                      nor(est[3,]-true),nor(est[4,]-true),
                      nor(est[5,]-true),nor(est[6,]-true))
        label <- which.min(distance)

        for(kk in 1:(3*k)){
          thetaa[j,kk] <- est[label,kk]
        }
      }

      if(k>=4){
        rl<-res$lambda; rm<-res$mu; rs<-res$sigma
        est<-c(rl[order(rm)],rm[order(rm)],rs[order(rm)])
        estt = t(matrix(est, nrow = k))

        for(kk in 1:(3*k)){
          thetaa[j,kk]=estt[kk]
        }
      }

      thetaa[j,3*k+1] <- res$loglik

      la2 <-0
      for(kk in 1:k){
        lakk=thetaa[j,3*kk-2]*pnorm(sort(xb),thetaa[j,3*kk-1],thetaa[j,3*kk])
        la2<-la2+lakk
        la2
      }

      cd[j] <- sum((la2-e_cd)^2)+(1/(12*n))

    }

    lab <- which.min(cd)

    theta <- c(thetaa[lab,1:(3*k)])

    cdf <-0
    for(kk in 1:k){
      cdfkk=theta[3*kk-2]*pnorm(sort(xb),theta[3*kk-1],theta[3*kk])
      cdf<-cdf+cdfkk
      cdf
    }


    DP <- max(abs(e_cdf-cdf))
    DM <- max(abs(cdf-e_cdf1))
    ks[i] <- max(DP,DM)
    cvm[i] <- sum((cdf-e_cd)^2)+(1/(12*n))
    kui[i] <- DP+DM
    wat[i] <- cvm[i] - n*((mean(cdf)-0.5)^2)
    ad[i] <- -n-((sum((2*(1:n)-1)*(log(cdf)+log(1-sort(cdf,decreasing = TRUE))))/n))
  }

  out=list(real=real,ks=ks,cvm=cvm,kui=kui,wat=wat,ad=ad)
  out
}

