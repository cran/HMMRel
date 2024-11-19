\name{fit.hmmR}
\alias{fit.hmmR}
\title{
Non-parametric fitting of a HMM using the Baum-Welch algorithm.
}
\description{
This function adapts the EM algorithm to fit the transition matrix of the hidden Markov chain as well as the emission probability matrix
of a HMM.
}
\usage{
fit.hmmR(Y,P0,M0,alpha0,max.iter=50,epsilon=1e-9,Nx,Ny)
}

\arguments{
  \item{Y}{A sequence of observations consisting of signals of system performance.}
  \item{P0}{A square matrix of dimension \code{Nx} with the transition probabilities of the hidden MC at the first iteration of the algorithm.}
  \item{M0}{A matrix of dimension \code{Nx} \eqn{\times} \code{Ny} with the emission probabilities at the first step of the algorithm.}
  \item{alpha0}{A vector of size \code{Nx} with the initial distribution of the hidden Markov chain at the first iteration of the algorithm.}
  \item{max.iter}{An integer value with the maximum number of iterations in the iterative algorithm. Default value is \code{max.ite}=50.}
  \item{epsilon}{A numeric value with the tolerance in the iterative algorithm. Default value is \code{epsilon}=1e-9.}
  \item{Nx}{An integer value with the maximum number of states of the hidden Markov chain.}
  \item{Ny}{An integer value with the maximum number of signals emitted by the system.}
}
\details{
\itemize{
\item{The argument \code{alpha0} representing the initial distribution of the hidden MC is fixed, and defined by the user. The argument \code{Nx} is the size of the state space of the hidden MC. As default, the set of numbers 1,...,\code{Nx} is the state space of the hidden MC.}
\item{\code{Ny} is the size of the alphabet of signals emitted.  As default, the set of numbers 1,...,\code{Ny} is the signal-alphabet.}
\item{The successive iterations of the algorithm can be traced and information is accessible from the outcome of this function.}
}}

\value{Among other information, this function provides the following values:
\item{P}{A square matrix with \code{Nx} rows with the estimated transition probabilities.}
\item{M}{A matrix of \code{Nx} rows and \code{Ny} columns, with the estimated emission probabilities.}
\item{n.iter}{An integer indicating the number of iterations performed in the algorithm.}
\item{tol}{A numeric value with the achieved tolerance value in the algorithm.}
\item{Fm}{A matrix of dimension \code{n} \eqn{\times} \code{Nx} with the estimated forward probability values at the last iteration of the algorithm, where \code{n} is the size of the observed vector \code{Y}. }
\item{Bm}{A matrix of dimension \code{Nx} \eqn{\times} \code{n} with the estimated forward probability values at the last iteration of the algorithm.}
\item{AIC}{The estimated value of the Akaike statistics. The number of parameters to be estimated is \code{nparam}=\code{Nx}*(\code{Nx}-1)+\code{Nx}*(\code{Ny}-1).}
}

\references{Gamiz, M.L., Limnios, N., and Segovia-Garcia, M.C. (2023). Hidden Markov models in reliability and maintenance. European Journal of Operational Research, 304(3), 1242-1255.
}

\author{M.L. Gamiz, N. Limnios, and M.C. Segovia-Garcia (2024)
}

\seealso{
See \code{def.hmmR} to define an object HMM, and \code{sim.hmmR} to simulate a random path from a given HMM object.}

\examples{
model<-'other'
rate<-NA
p<-NA
P<-matrix(c(0.7,0.3,1,0),2,2,byrow=TRUE)
M<-matrix(c(0.6,0.4,0,0,0,1),2,3,byrow=TRUE)
alpha<-c(1,0)
Nx<-2
Ny<-3
n.up<-1
n.green<-2
hmm0<-def.hmmR(model,rate,p,alpha,P,M,Nx,Ny,n.up,n.green)
set.seed(1969)
datos<-sim.hmmR(hmmR=hmm0,n=10)
estim<-fit.hmmR(Y=datos$Yn,P0=P,M0=M,alpha0=alpha,max.iter=50,epsilon=1e-9,Ny=3,Nx=2)
estim$P;P
estim$M;M
}
