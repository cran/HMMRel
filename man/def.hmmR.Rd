\name{def.hmmR}
\alias{def.hmmR}

\title{
Define a HMM object for Reliability Analysis.
}

\description{
This function creates a list with all the elements that describe a HMM in the context of Reliability and Maintenance.
}

\usage{
def.hmmR(model,rate,p,alpha,P,M,Nx,Ny,n.up,n.green)
}

\arguments{
\item{model}{ A character string to choose which HMM model is considered. Possible values are "KooN", "shock", "other"}
\item{rate}{ A positive real number indicating the failure rate of one unit of the system.}
\item{p}{ A real number in the interval (0,1) indicating the probability that the system receives one shock during two consecutive inspections.}
\item{alpha}{ A vector of size \code{Nx} with the initial law of the hidden Markov chain.}
\item{P}{ A square matrix of dimension \code{Nx} with the transition probabilities between the hidden states.}
\item{M}{ A matrix of dimension \code{Nx} \eqn{\times} \code{Ny} with the emission probabilities.}
\item{Nx}{ An integer indicating the total number of states in the system. By default the states are labelled: 1,...,\code{Nx}.}
\item{Ny}{ An integer indicating the total number of signals received.  By default the signals are labelled: 1,...,\code{Ny}.}
\item{n.up}{ An integer lower than \code{Nx} indicating the number of (hidden) operative states in the system. The first \code{n.up} values in the state set denote operative states.}
\item{n.green}{ An integer lower then \code{Ny} indicating the number of signals of good performance. The first \code{n.green} signals are read as good performance of the system.}
}

\details{
\itemize{
\item When \code{model}=\code{"KooN"} the argument \code{Nx} is the maximum number of units in the system. There must be  \code{K}=\code{n.up} operative units for the system to function.
If \code{K}=1 a parallel system is built. If \code{K}=N a series system is built.
\item When \code{model}=\code{"shock"} the argument \code{Nx} is the maximum number of shocks that the system can accumulate before breakdown.
}}

\value{
A list with the elements of the HMM.
\item{states}{A set of \code{Nx} characters or integers decribing the hidden states of the system.}
\item{signals}{A set of \code{Ny} characters or integers decribing the possible signals observed.}
\item{P}{A square matrix with \code{Nx} rows with the transition probabilities between the hidden states.}
 \item{M}{A matrix of dimension \code{Nx} \eqn{\times} \code{Ny} with the emission probabilities.}
}
\references{Gamiz, M.L., Limnios, N., and Segovia-Garcia, M.C. (2023). Hidden Markov models in reliability and maintenance. European Journal of Operational Research, 304(3), 1242-1255.
}
\author{M.L. Gamiz, N. Limnios, and M.C. Segovia-Garcia (2024)
}

\seealso{
See also \code{sim.hmmR} to simulate data from a given HMM.
}

\examples{
## Define a HMM object describing a repairable system
## The system can be in one of 3 states: 2 up states and 1 down state.
## 3 different signals can be received: 2 good performance signals (green)
##  and 1 signal of failure (red)
P<-matrix(c(8,2,1,0,6,4,6,2,2)/10,3,3,byrow=TRUE)
M<-matrix(c(7,3,0,4,3,3,0,4,6)/10,3,3,byrow=TRUE)
Nx<-3; Ny<-3
n.up<-2; n.green<-2
alpha<-c(1,0,0)
hmm1<-def.hmmR(model='other',rate=NA,p=NA,alpha=alpha,P=P,M=M,Nx=Nx,Ny=Ny,
               n.up=n.up,n.green=n.green)
hmm1
}
