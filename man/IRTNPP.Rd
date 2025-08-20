\name{IRTNPP}
\alias{IRTNPP}
\title{MCMC Sampling for IRT Model Ability Parameters using Normalized Power Prior}

\description{
Conduct posterior sampling for IRT model ability parameters with normalized power prior.

For the power parameter \eqn{\delta}, a Metropolis-Hastings algorithm with either 
independence proposal, or a random walk proposal on its logit scale is used.

For the model parameters \eqn{\beta}, a Metropolis-Hastings algorithm with either 
normal proposal, or uniform proposal is used.
}

\usage{
IRTNPP(y, dseq, prior_mu, prior_sd, MCsize, disa, difa1, difa2,
       cut, prior_beta, prior_delta, disb, difb1, difb2,
       prop_delta, rw_delta, rw_n_beta, rw_u_beta, ind_delta,
       prop_beta, n_sample, burnin, thin)
}

\arguments{
\item{y}{a vector that contains historical data and current data, where the first half  
         consists of historical data and the second half consists of current data.}
\item{dseq}{numeric vector or scalar between 0 and 1. The value of \eqn{\delta}.}
\item{prior_mu}{the prior mean of each ability parameter \eqn{\beta}.}
\item{prior_sd}{the prior standard deviation of each ability parameter \eqn{\beta}.}
\item{MCsize}{positive integer. Sample size of importance sampling.}
\item{disa}{a matrix of item discriminability parameters in historical data.}
\item{difa1}{a vector of the first difficulty parameter of items in historical data.}
\item{difa2}{a vector of the second difficulty parameter of items in historical data.}
\item{cut}{critical value between 0 and 1. If \eqn{\delta} is less than or equal to this
           value, select the initial prior of \eqn{\beta} as the proposed distribution
           of importance sampling; Otherwise, select the posterior distribution
           with historical data of \eqn{\beta} as the proposed distribution.}
\item{prior_beta}{list. Parameters of normal prior for \eqn{\beta}.}
\item{prior_delta}{list. Parameters of beta prior for for \eqn{\delta}.}
\item{disb}{a matrix of item discriminability parameters in current data.}
\item{difb1}{a vector of the first difficulty parameter of items in current data.}
\item{difb2}{a vector of the second difficulty parameter of items in current data.}
\item{prop_delta}{character. The class of proposal distribution for \eqn{\delta}.}
\item{rw_delta}{numeric. The stepsize(variance of the normal distribution) for the random walk 
                proposal of logit \eqn{\delta}. Only applicable if prop_delta = 'RW'.}
\item{prop_beta}{character. The class of proposal distribution for \eqn{\beta}.}
\item{rw_n_beta}{numeric vector. Standard deviation of proposed distribution of for \eqn{\beta}.}
\item{rw_u_beta}{numeric. rw_u_beta*2 is the interval length of uniform distribution.}
\item{ind_delta}{numeric vector. Two parameters when the proposed distribution of 
                 \eqn{\delta} is beta distribution.}
\item{n_sample}{positive integer. Specifies the number of posterior samples in the output.}
\item{burnin}{positive integer. The output will only show MCMC samples after bunrin.}
\item{thin}{positive integer. The thinning parameter in MCMC sampling.}
}

\details{
This function needs three additional R packages: KernSmooth, msm, mvtnorm.

This function needs two additional R functions: makePositiveDefinite, Metro_Hastings.

The outputs include the posterior estimates of the ability parameters of the IRT model and 
power parameter, as well as the acceptance rates in sampling \eqn{\delta} and \eqn{\beta}.
}

\value{
A vector consisting of 5 parts: 
the acceptance rate in MCMC sampling for \eqn{\beta} and \eqn{\delta} using Metropolis-Hastings algorithm, 
the posterior mean of \eqn{\beta} and \eqn{\delta}, 
the posterior standard deviation of \eqn{\beta} and \eqn{\delta}, 
the posterior median of \eqn{\beta} and \eqn{\delta}, and 
the posterior mode of power parameter \eqn{\delta}.
}

\references{
Chalmers, R.P. (2012).
mirt: A multidimensional item response theory package for the R environment.
\emph{Journal of Statistical Software 48:1–29.}

Matteucci, M., Veldkamp, B. (2015).
The approach of power priors for ability estimation in IRT models.
\emph{Qual Quant 49:917–926.}

Han, Z., Zhang, Q., Wang, M., Ye, K., Chen, M.H. (2023).
On efficient posterior inference in normalized power prior Bayesian analysis.
\emph{Biometrical Journal  65:2200194.}
}

\author{
Qiang Zhang \email{zqzjf0408@163.com}
}
