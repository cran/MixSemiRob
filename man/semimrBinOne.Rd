% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semiregBinom.R
\name{semimrBinOne}
\alias{semimrBinOne}
\title{Semiparametric Mixture of Binomial Regression with a Degenerate Component
with Constant Proportion and Time-Varying Success Probability (One-step Backfitting)}
\usage{
semimrBinOne(t, x, N, tg = NULL, tune = 1, tol = 1e-02)
}
\arguments{
\item{t}{a vector of time variable along which \eqn{p(t)} varies.}

\item{x}{a vector of observed number of successes. The length of \code{t} and \code{x} must be the same.}

\item{N}{a scalar, specifying the number of trials for the Binomial distribution.}

\item{tg}{grid points of time used in the kernel regression for the estimation of \eqn{p(t)}.
Default is NULL, and 100 equally spaced grid points will automatically generated
using the minimum and maximum values of \code{t}.}

\item{tune}{a scalar, specifying the percentage of data included in local estimation.
related to the bandwidth selection and local estimation. Default is 1.}

\item{tol}{stopping criteria for the algorithm.}
}
\value{
A list containing the following elements:
\item{pt}{estimated time-varying success probabilities for the first component.}
\item{w}{estimated constant proportion for the first component.}
\item{h}{bandwidth for the kernel regression. The bandwidth calculation can be found in Section 4 of Cao and Yao (2012).}
}
\description{
`semimrBinOne' implements the one-step backfitting method (Cao and Yao, 2012)
for semiparametric estimation of a mixture of binomial distributions with one degenerate component,
with constant proportion and time-varying success probability.
}
\details{
The semiparametric mixture of binomial regression model is as follows:
\deqn{w \times B(N,p(t))+(1-w)\times B(N,0),}
where \eqn{B(N,p)} is the probability mass function of a binomial distribution
with the number of trials \eqn{N} and the success probability \eqn{p}.
Here, the second component is a degenerate distribution with mass 1 on 0.
The time-varying success probability \eqn{p(t)} for the binomial components
are estimated by the kernel regression using one-step estimation for faster computation with some bandwidth.
}
\examples{
nobs = 50
tobs = seq(from = 0, to = 1, length.out = nobs)
pi1Tru = 0.4
ptTru = 0.3 * (1.5 + cos(2 * pi * tobs))
nfine = nobs
tfine = seq(from = 0, to = 1, length.out = nfine)
b = rbinom(nobs, size = 1, pi1Tru)
yobs = apply(X = matrix(ptTru), 1, rbinom, n = 1, size = 7)
yobs = ifelse(b == 1, 0, yobs)
ftonestep = semimrBinOne(t = tobs, x = yobs, N = 7, tg = tfine)
}
\references{
Cao, J. and Yao, W. (2012). Semiparametric mixture of binomial regression with a degenerate component.
Statistica Sinica, 27-46.
}
\seealso{
\code{\link{semimrBin}}, \code{\link{semimrBinFull}}
}
