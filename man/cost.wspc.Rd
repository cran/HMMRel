\name{cost.wcspc}
\alias{cost.wspc}
\title{
Maintenance Policy based on Warning Signal Probability Criteria.
}

\description{
Preventive maintenance based on Warning Signal Probability Criteria (CSPC).
}

\usage{
cost.wspc(prob,hmmR,n.up1,n.green1,cost.C,cost.P,t.max)
}
\arguments{
  \item{prob}{ A real number in the interval (0,1).}
  \item{hmmR}{ A hidden Markov Model.}
  \item{n.up1}{ An integer value for indicating the total number of optimal performance states of the hidden MC. }
  \item{n.green1}{ An integer value for indicating the total number of safe signals. A safe signal indicates an optimal system performance.}
  \item{cost.C}{ A positive real number denoting the cost value in monetary units incurred by a corrective maintenance action.}
  \item{cost.P}{ A positive real number denoting the cost value in monetary units incurred by a preventive maintenance action.}
  \item{t.max}{ A time value for the maximum time the system will be in use. After that time the system will not operate anymore.}
}

\details{
Preventive maintenance policies based on Warning Signal Probability criteria (CSPC) is considered.
Roughly speaking, a preventive maintenance action is carried out at a time when the probability that system is not occupying an optimal operational
state or a warning signal is received is above a prespecified value \code{prob}.
The subset of operative states \code{up} is in turn split into two subsets:
optimal states or \code{up1} and operative but critical states or \code{up2}, where \code{up}=\code{up1}\eqn{\cup}\code{up2}.
Similarly, the set of \code{green} signals is split into two subsets: \code{safe} signals and \code{warning} signals.
\code{n.green1} is the size of subset \code{safe}.
For a given probability value (\code{prob}) this function first calculates the optimal inspection time
\deqn{\code{t.insp}=\min \{t>0: \Pr( Y(t) \in \text{warning}, Y(u) \in \text{safe}, \forall u \leq t )\}.}
The system is inspected every \code{t.insp} units.of time. At the time of inspection, any of three situations can be found:
\enumerate{
\item the system is in failure, then the system is returned to operational conditions (\code{up1}), and a cost of \code{cost.C} monetary-units is implied;
\item the system is in a state of \code{up2}, then a preventive maintenance action is carried out, returning the system to a state in \code{up1} and implying a cost of \code{cost.P} monetary-units; and
\item the system is found in a state of \code{up1}, then no maintenance action is carried out and there is no associated cost.
}}

\value{
\item{time.insp }{The time at which preventive maintenance is carried out.}
\item{t.max }{The maximum time that the system is being used.}
\item{n.insp}{The total number of inspections that are carried out during the system lifetime, i.e. in the interval (0, \code{t.max}).}
\item{total.cost}{The total cost incurred by all maintenance actions (corrective and preventive) developed in the system.}
}

\references{Gamiz, M.L., Limnios, N., and Segovia-Garcia, M.C. (2023). Hidden Markov models in reliability and maintenance. European Journal of Operational Research, 304(3), 1242-1255.
}

\author{M.L. Gamiz, N. Limnios, and M.C. Segovia-Garcia (2024)
}

\seealso{
See \code{cost.cspc} for the implementation of the CSPC algorithm for maintenance policy.
}
\examples{
model<-'other'
rate<-p<-NA
P<-matrix(c(8,2,1,0,6,4,6,2,2)/10,3,3,byrow=TRUE)
M<-matrix(c(7,3,0,4,3,3,0,4,6)/10,3,3,byrow=TRUE)
Nx<-3; Ny<-3
n.up<-2; n.green<-2
alpha<-c(1,0,0)
hmm1<-def.hmmR(model=model,rate=NA,p=NA,alpha=alpha,P=P,M=M,Nx=Nx,Ny=Ny,n.up=n.up,n.green=n.green)
prob<-0.8;
n.up1<-n.green1<-1;cost.C<-10;cost.P<-5;t.max<-50
cost1<-cost.wspc(prob=prob,hmmR=hmm1,n.up1=n.up1,n.green1=n.green1,
         cost.C=cost.C,cost.P=cost.P,t.max=t.max)
cost1
#
v.prob<-seq(0.1,0.99,length=100)
v.cost1<-inspection.time<-double(100)
for(i in 1:100)
{cost<-cost.wspc(prob=v.prob[i],hmmR=hmm1,n.up1=n.up1,n.green1=n.green1,
cost.C=cost.C,cost.P=cost.P,t.max=t.max)
v.cost1[i]<-cost$total.cost
#inspection.time[i]<-cost$time.insp
}
oldpar<-par(mar=c(5,5,10,10))
plot(v.prob,v.cost1,type='s',main='WSPC Algorithm for Maintenance Policy',
      xlab='Probability of critical state',
     ylab='Cost of maintenance')
grid()
par(oldpar)
}

