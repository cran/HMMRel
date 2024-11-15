\name{Rcalc.hmmR}
\alias{Rcalc.hmmR}

\title{
Calculate the reliability of a system based on HMM.
}

\description{
For a given time \code{t} this function returns the value of the probability that the system does not fail in the interval \code{(0,t]}.
It gives the probability that the system survives and is still working beyond time \code{t}.}
\usage{
Rcalc.hmmR(hmmR,t)
}

\arguments{
   \item{hmmR}{ A Hidden Markov Model.}
   \item{t}{ A value of time, it must be an integer equal or greater than 0.}
}

\details{
The state space is split into two subsets, i.e. \code{states}=\code{up} \eqn{\cup} \code{down}. The subset \code{up} contains the states of good functioning, while the subset \code{down} contains the failure states.
The signals aphabet is split into two subsets, i.e. \code{signals}= \code{green} \eqn{\cup} \code{red}.
A \code{green}-signal indicates good performance of the system, while a \code{red}-signal alerts of something wrong in the system.
This function returns the probability that the system has not entered the set of \code{down} states or any signal from the \code{red} subset of signals has been emitted at any time before \code{t}.
}
\value{
This function returns the probability that the system is working through a state in the \code{up} subset, and a \code{green} signal is being received.
If \code{t}=0, then the returned value is 1.
}
\references{Gamiz, M. L., Limnios, N., and Segovia-Garcia, M.C. (2023). Hidden Markov models in reliability and maintenance. European Journal of Operational Research, 304(3), 1242-1255.
}
\author{M.L. Gamiz, N. Limnios, and M.C. Segovia-Garcia (2024)
}

\seealso{
See \code{def.hmmR} to define a HMM object.
}

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
hmm0<-def.hmmR(model=model,rate=NA,p=NA,alpha=alpha,P=P,M=M,Nx=Nx,Ny=Ny,n.up=n.up,n.green=n.green)
times<-0:30
Rt<-Rcalc.hmmR(hmmR=hmm0,t=times)
oldpar <- par(mar = c(5, 5, 10, 10))
plot(times,Rt,type='s',ylim=c(0,1),ylab='',xlab='time',main='Reliability based on HMM')
grid()
par(oldpar)
}

