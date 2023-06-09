% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/freq.R
\name{mixnorm}
\alias{mixnorm}
\title{Normal Likelihood Based Method for Label Switching}
\usage{
mixnorm(x, k = NULL, ini = NULL, acc = NULL)
}
\arguments{
\item{x}{The observations. See details below.}

\item{k}{k components (2 by default); k can also be a vector. See details below.}

\item{ini}{List of initial value of mu, pi, sigma (when k>2). \cr Example: ini = list(mu=c(0,1),pi=c(0.5,0.5), sigma=c(1,1)),
If left as NULL, initial values will be calculated from random sample from observations x by default.}

\item{acc}{Stopping criteria. When x is one dimensional, the default value of acc is
0.001; When x is multidimensional the stopping criteria is 10^-5.}
}
\value{
List of estimations:
mu: estimated means of each group;
pi: estimated proportion of each group;
sigma: estimated sd of each group;
p: the classification probability matrix;
lh: the likelihood value of the estimations;
id: which method we use for initial value;
inilh: likelihood of five random starting point. The program choose the point with largest
likelihood to continue run the estimation process until converge;
iniid: Which one of the initial point has the largest likelihood;
}
\description{
Solving the label switching problem by using Em algorithm to estimate normal
mixture model with equal variance. The function support one dimensional and high
dimensional observations.
}
\details{
If the observation is one dimensional

When k (k=2 by default) is a scalar:
Input x: the data vector, k: k components (2 by default);
Output: 3 by k matrix, each column corresponding to one component and 1st line:mu,2nd:pi,3rd:sigma.

When k is a vector:
Input k:the standard deviation vector of k=length(k) components.
Output: 2 by k matrix, each column corresponding to one component and 1st line:mu, 2nd:pi.

If the observation is high dimensional

When k (k=2 by default) is a scalar:
Input x: the data matrix with each row corresponding to a observation, the dimension of data s is column dimension of x
Output: (k+s+k) by s matrix, first k rows are component means, next s rows are sigma matrix, the last k rows are pi estimation

When k is a vector:
Input k:the standard deviation matrix of all the components with number of components k=k/s
Output: (k+k) by s matrix,first k rows are component means, next k rows are pi estimation
}
\examples{
#See example in complhfrequency
}
