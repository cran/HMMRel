\name{sim.hmmR}
\alias{sim.hmmR}

\title{
Simulate sequence of states and signals of functioning of a system modelled by a HMM.
}

\description{
This function simulates a sample path from a 2-dimensional HMM. It returns the hidden sequence of states and signals. At each time, the hidden state of the system is simulated from the HMM as well as the associated signal that informs on the system performance at that time.
}

\usage{
sim.hmmR(hmmR,n)
}

\arguments{
  \item{hmmR}{ A Hidden Markov Model.}
  \item{n}{ An integer number indicating the length of the sequence of states and signals to be simulated.}
}

\value{
The function \code{sim.hmmR} returns a list with the following information:
\item{Xn }{The sequence of simulated hidden states.}
\item{Yn }{The sequence of observed signals.}
\item{P}{The transition probability matrix of the hidden Markov chain (MC).}
\item{alpha}{The initial distribution of the hidden MC.}
\item{M}{The emission probability matrix.}
\item{states}{The set of hidden states of the system.}
\item{signal}{The alphabet corresponding to the observations.}
}

\references{Gamiz, M.L., Limnios, N., Segovia-Garcia, M.C. (2023). Hidden Markov models in reliability and maintenance. European Journal of Operational Research, 304(3), 1242-1255.
}

\author{M.L. Gamiz, N. Limnios, and M.C. Segovia-Garcia (2024)
}
\seealso{
See \code{def.hmmR} to define a HMM object.
}

\examples{
## Define a HMM object describing a repairable system
P<-matrix(c(8,2,1,0,6,4,6,2,2)/10,3,3,byrow=TRUE)
M<-matrix(c(7,3,0,4,3,3,0,4,6)/10,3,3,byrow=TRUE)
hmm1<-def.hmmR(model='other',rate=NA,p=NA,alpha=c(1,0,0),P=P,M=M,Nx=3,Ny=3,n.up=2,n.green=2)
sim.hmmR(hmmR=hmm1,n=20)
}
