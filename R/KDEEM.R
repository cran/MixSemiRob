############################################################################
# Yao Algorithm 2 kdeem
############################################################################
#' Kernel Density-based EM type algorithm with Least Square Estimation for Semiparametric Mixture Regression with Unspecified Homogenous Error Distributions
#'
#' For kernel density-based semiparametric mixture regression for homogeneous error distributions
#' Using weighted least square estimation for beta in M step, where weight is the posterior probability of observation comes from each component.
#' And Gaussian kernel for g(.) component estimation
#'
#' @param x explanatory variables in matrix from with row as observations
#' @param y response variable
#' @param M number of component
#' @param ini initial values for parameters. ini = list(beta = regression coefficients,
#' posterior = matrix: probability of each observation is from component m)
#' the initial value can be obtained by applying traditional MLE with EM algorithm and normal assumption.
#'
#'
#' @return
#' P: posterior probability of each observation is from component m
#' Beta: estimated regression coefficients
#' Tao: Estimated precision parameter, the inverse of standard deviation.
#' PI: estimated component proportion
#' H: bandwidth of the kernel estimation
#' @importFrom ucminf ucminf
#' @importFrom stats bw.SJ
#' @importFrom Rlab rbern
#' @export
#' @examples
#' # See example in KDEEM
KDEEM.LSE<-function(x,y,M,ini){
  n=length(y);  X=cbind(1,x); Dimen=dim(X)[2]
  beta.matrix.new=ini$beta
  beta.matrix.old=matrix(0,Dimen,M)
  P=ini$posterior;
  PI=colMeans(P)

  Eps.matrix=matrix(0,n,M)
  for(j in 1:M){Eps.matrix[,j]=y-X%*%beta.matrix.new[,j]}
  ############################################################################
  # use EM algorithm to calculate the mle
  numiter=0;
  GR.matrix=matrix(0,n,M);

  while ( norm(beta.matrix.new-beta.matrix.old,'2')>1e-2)
  {
    beta.matrix.old=beta.matrix.new;

    Tao=rep(0,M)
    R.matrix=matrix(0,n,M)
    for(k in 1:M){
      Tao[k]=( sum(P[,k]*(Eps.matrix[,k]^2))/sum(P[,k])   )^(-0.5) ;
      R.matrix[,k]=Eps.matrix[,k]*Tao[k] ;
    }
    ###############################################
    Residual.1=R.matrix[,1][which(P[,1]>=P[,2])]
    Residual.2=R.matrix[,2][which(P[,1]<=P[,2])]
    H=bw.SJ(c(Residual.1,Residual.2));
    ###############################################

    for(j in 1:M){
      for(i in 1:n){
        ###############################
        GR.matrix[i,j]=( t(P[,1])%*%dnorm((R.matrix[,1]-R.matrix[i,j])/H)/H+t(P[,2])%*%dnorm((R.matrix[,2]-R.matrix[i,j])/H)/H )/n
        ###############################
      }
    }


    for(j in 1:M){
      PI[j]=mean(P[,j]);
      P[,j]=(PI[j]*GR.matrix[,j]*Tao[j])/(PI[1]*GR.matrix[,1]*Tao[1]+PI[2]*GR.matrix[,2]*Tao[2]);
      Wj=diag(P[,j]);
      beta.matrix.new[,j]=solve(t(X)%*%Wj%*%X)%*%t(X)%*%Wj%*%y;
      Eps.matrix[,j]=y-X%*%beta.matrix.new[,j];
    }
    numiter=numiter+1;
    if(numiter>200){break}
  } # End of the loop while.
  ############################################################################
  return(list(P=P,Beta=beta.matrix.new,Tao=Tao,PI=PI,H=H))
}
#KDEEM.LSE





############################################################################
# Yao Algorithm 3
############################################################################
#' Kernel Density-based EM type algorithm for Semiparametric Mixture Regression with Unspecified Homogenous Error Distributions
#'
#' For kernel density-based semiparametric mixture regression for homogeneous error distributions
#' @param x explanatory variables in matrix from with row as observations
#' @param y response variable
#' @param M number of component
#' @param ini initial values for parameters. ini = list(Beta = regression coefficients,
#' P = matrix: probability of each observation is from component m based on initial beta and H,
#' H = bandwidth for kernel estimation)
#' @param maxiter maximum iteration for the algorithm
#'
#' @return
#' P: posterior probability of each observation is from component m
#' Beta: estimated regression coefficients
#' PI: estimated component proportion
#' H: bandwidth of the kernel estimation
#' @export
#'
#' @examples
#' # See example in KDEEM
KDEEM.H<-function(x,y,M,ini,maxiter)
{
  n=length(y);  X=cbind(1,x); Dimen=dim(X)[2]
  beta.matrix.new=ini$Beta
  beta.matrix.old=matrix(0,Dimen,M)
  P=ini$P;H=ini$H;

  Eps.matrix=matrix(0,n,M) #residual matrix
  PI=c()
  for(k in 1:M){
    PI[k]=mean(P[,k]);
    Eps.matrix[,k]=y-X%*%beta.matrix.new[,k];
    }

  numiter=0
  while ( norm(beta.matrix.new-beta.matrix.old,'f')>1e-2 & numiter<maxiter )
    {
    numiter=numiter+1;
    beta.matrix.old=beta.matrix.new;

    grid=list();
    for(j in 1:M){grid[[j]]=seq(min(Eps.matrix[,j])-3*H,max(Eps.matrix[,j])+3*H,by=H/10);}
    ###############################OK
    gfungrid=list()
    for(j in 1:M){
      g.vector=vector()
      for(i in 1:length(grid[[j]])){
           g.vector[i]=( t(P[,j])%*%dnorm((Eps.matrix[,j]-grid[[j]][i])/H)/H )/sum(P[,j]);
      }
       gfungrid[[j]]=g.vector;
    }
    ###############################

    #M step
    for(j in 1:M){
      estg<-function(x){ out=approx(grid[[j]],gfungrid[[j]],x)$y
                         out[x>max(grid[[j]])]=0
                         out[x<min(grid[[j]])]=0
                         return(out)
                         }
      fn<-function(bet){out=-P[,j]%*%log(estg((y-X%*%bet)))}

      beta.matrix.new[,j]=ucminf(beta.matrix.new[,j],fn)$par
      Eps.matrix[,j]=y-X%*%beta.matrix.new[,j]
      P[,j]=PI[j]*estg(Eps.matrix[,j]);
    }
    P=P/matrix(rep(rowSums(P),M),ncol=M)
    PI=colMeans(P);
  }
  # print(beta.matrix.new)
  # print(numiter)
  ############################################################################
  return(list(P=P,Beta=beta.matrix.new,PI=PI,H=H))
}
#KDEEM.H


















############################################################################
# Yao Algorithm 4
############################################################################
#' Kernel Density-based EM type algorithm for Semiparametric Mixture Regression with Unspecified Error Distributions
#'
#' Unlike KDEEM.LSE and KDEEM.H, KDEEM is used for unspecified error distributions, can be applied for unspecified
#' error distribution (both homogeneous and heterogenous)
#' @param x explanatory variables in matrix from with row as observations
#' @param y response variable
#' @param M number of component
#' @param ini initial values for parameters. ini = list(Beta = regression coefficients,
#' P = matrix: probability of each observation is from component m based on initial beta and H,
#' Tao =  precision parameter, the inverse of standard deviation
#' PI = component proportion
#' H = bandwidth for kernel estimation)
#' the initial value can be obtained by applying traditional MLE with EM algorithm and normal assumption.
#' @param maxiter maximum iteration for the algorithm
#'
#' @return
#' P: posterior probability of each observation is from component m
#' Beta: estimated regression coefficients
#' Tao: Estimated precision parameter, the inverse of standard deviation.
#' PI: estimated component proportion
#' H: bandwidth of the kernel estimation
#' @export
#'
#' @examples
#' Samplesize=300;
#'M=2 # M is the number of groups.
#'Dimen=2
#'Beta.true.matrix<-matrix(c(-3,3,3,-3),Dimen,M);
#'PI.true=c(0.5,0.5);
#'#simulate data
#'x=runif(Samplesize)
#'X=cbind(1,x);
#'n=Samplesize;
#'Group.ID=Rlab::rbern(n, prob=0.5);
#'Error=rnorm(n,0,1);
#'n1=sum(Group.ID);
#'n2=n-n1;
#'
#'y=rep(0,Samplesize);
#'err=rep(0,n);
#'
#'for(i in 1:n){
#'  if(Group.ID[i]==1)
#'  {err[i]=Error[i];
#'  y[i]=X[i,]%*%Beta.true.matrix[,1]+err[i];}
#'  else
#'  {err[i]=0.5*Error[i];
#'  y[i]=X[i,]%*%Beta.true.matrix[,2]+err[i];}
#'}
#'
#'# Get initial value from MLE
#'Result.MLE= mixtools::regmixEM(y, x, arbmean = TRUE, arbvar = TRUE, epsilon = 1e-4);
#'
#'Result.KDEEM.LSE = KDEEM.LSE(x,y,M,Result.MLE)
#'
#'Result.KDEEM.H = KDEEM.H(x,y,M,Result.KDEEM.LSE,maxiter=3)
#'
#'Result.KDEEM = KDEEM(x,y,M,Result.KDEEM.LSE,maxiter=3)
KDEEM<-function(x,y,M,ini,maxiter)
{
  n=length(y);  X=cbind(1,x); Dimen=dim(X)[2]
  beta.matrix.new=ini$Beta
  beta.matrix.old=matrix(0,Dimen,M)
  P=ini$P;H=ini$H;Tao=ini$Tao

  Eps.matrix=matrix(0,n,M) #residual matrix
  R.matrix=matrix(0,n,M)
  PI=c()
  for(k in 1:M){
    PI[k]=mean(P[,k]);
    Eps.matrix[,k]=y-X%*%beta.matrix.new[,k]
    R.matrix[,k]=Eps.matrix[,k]*Tao[k]; #standarized residual
  }
  numiter=0
  while ( norm(beta.matrix.new-beta.matrix.old,'f')>1e-2 & numiter<maxiter )
    { numiter=numiter+1;
    beta.matrix.old=beta.matrix.new;
    grid=seq(min(R.matrix)-3*H,max(R.matrix)+3*H,by=H/10)
    ###############################
    gfungrid=c();
    for(i in 1:length(grid))
      gfungrid[i]=( t(P[,1])%*%dnorm((R.matrix[,1]-grid[i])/H)/H+t(P[,2])%*%dnorm((R.matrix[,2]-grid[i])/H)/H )/n
    ###############################
    estg<-function(x){   out=approx(grid,gfungrid,x)$y
                         out[x>max(grid)]=0
                         out[x<min(grid)]=0
                         return(out)
                         }
    #M step
    for(j in 1:M){
      fn<-function(bet)
      {out=-P[,j]%*%log(estg((y-X%*%bet)*Tao[j]))}
      beta.matrix.new[,j]=ucminf(beta.matrix.new[,j],fn)$par
      Eps.matrix[,j]=y-X%*%beta.matrix.new[,j]
      Tao[j]=( sum(P[,j]*(Eps.matrix[,j]^2))/sum(P[,j]))^(-0.5) ;
      R.matrix[,j]=Eps.matrix[,j]*Tao[j];
      P[,j]=PI[j]*Tao[j]*estg(R.matrix[,j]);
    }
    P=P/matrix(rep(rowSums(P),M),ncol=M)
    PI=colMeans(P);
  }
  # print(beta.matrix.new)
  # print(numiter)
  ############################################################################
  return(list(P=P,Beta=beta.matrix.new,Tao=Tao,PI=PI,H=H))
}
#KDEEM



