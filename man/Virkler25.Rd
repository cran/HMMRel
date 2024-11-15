\name{Virkler25}
\alias{Virkler25}
\docType{data}

\title{
Fatigue crack growth in materials: Virkler dataset (tests 1 to 25)
}
\description{The data consist of an aluminum alloy specimen that was tested to investigate fatigue crack propagation.
Starting from an initial crack of length 9 mm for a particular item in test, the number of cycles for the size of the crack to reach a predetermined value was recorded successively. That is, it is registered the number of cycles every time an increment of size 0.2 mm in length occurs. The experiment finishes once a critical size of the crack is reached, meaning the failure of the item. The data were first published in Virkler et al. (1979) where
there were 68 specimens tested to grow the initial crack of 9 mm to the final crack of 50 mm. The first 25 tests are included.}
\format{
  A data frame with 26 variables:
  \describe{
    \item{CrackLength}{Length of the crack in the material.}
    \item{CycleCount1}{Cycle count for the first test.}
    \item{CycleCount2}{Cycle count for the second test.}
    \item{CycleCount3}{Cycle count for the third test.}
    \item{CycleCount4}{Cycle count for the fourth test.}
    \item{CycleCount5}{Cycle count for the fifth test.}
    \item{CycleCount6}{Cycle count for the sixth test.}
    \item{CycleCount7}{Cycle count for the seventh test.}
    \item{CycleCount8}{Cycle count for the eighth test.}
    \item{CycleCount9}{Cycle count for the ninth test.}
    \item{CycleCount10}{Cycle count for the tenth test.}
    \item{CycleCount11}{Cycle count for the eleventh test.}
    \item{CycleCount12}{Cycle count for the twelfth test.}
    \item{CycleCount13}{Cycle count for the thirteenth test.}
    \item{CycleCount14}{Cycle count for the fourteenth test.}
    \item{CycleCount15}{Cycle count for the fifteenth test.}
    \item{CycleCount16}{Cycle count for the sixteenth test.}
    \item{CycleCount17}{Cycle count for the seventeenth test.}
    \item{CycleCount18}{Cycle count for the eighteenth test.}
    \item{CycleCount19}{Cycle count for the nineteenth test.}
    \item{CycleCount20}{Cycle count for the twentieth test.}
    \item{CycleCount21}{Cycle count for the twenty-first test.}
    \item{CycleCount22}{Cycle count for the twenty-second test.}
    \item{CycleCount23}{Cycle count for the twenty-third test.}
    \item{CycleCount24}{Cycle count for the twenty-fourth test.}
    \item{CycleCount25}{Cycle count for the twenty-fifth test.}
     }
}

\references{
Virkler, D. A., Hillberry, B. M., and Goel, P. K. (1979). The statistical nature of fatigue
crack propagation. Journal of Engineering Materials and Technology, 101 , 148–153 .
}

\examples{
data(Virkler25)
i<-1 ## choose specimen number 1
ti<-Virkler25[,(i+1)]/10000  ## cycles at which the cracksize increases 0.2 units.
yt<-Virkler25[,1] ##
### the system is observed every 2000 cycles (0.2 unit times)
### observations:
t.obs<-seq(0,max(ti),by=0.2)
yi<-approx(x=ti,y=yt,xout=t.obs,method='linear',rule=2)$y
yi<-diff(yi)
#discretize the observations:
yi<-kmeans(yi,4)$cluster ## consider an alphabet of 4 signals
Nx<-2; # consider 2 hidden states
Ny<-4
alpha0<-c(1,0)
estim<-fit.hmmR(Y=yi,P0=NA,M0=NA,alpha0=alpha0,max.iter=50,epsilon=1e-9,Nx=Nx,Ny=Ny)
estim$P
estim$M
}
\keyword{datasets}
