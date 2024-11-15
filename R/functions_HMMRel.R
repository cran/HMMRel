##################

def.hmmR<-function(model,rate,p,alpha,P,M,Nx,Ny,n.up,n.green)
{

if (model=='KooN'){
   p=NA
   states<-1:(Nx) ## state=number of components working
   if(is.na(alpha[1])){alpha<-c(1,rep(0,Nx-1))}
   P<-matrix(0,Nx,Nx);q<-1-exp(-rate)
   for (i in 0:(Nx-1)) ##
     {for (j in i:Nx)
     {
       P[i,j]<-(factorial(Nx-i)/(factorial(j-i)*factorial(Nx-j)))*q^(j-i)*(1-q)^(Nx-j)}}
   P[Nx,1]<-1
   down<-(n.up+1):Nx
   up<-1:n.up;
   K<-n.up
   KooN<-list(states=states,K=K,alpha=alpha,P=P,rate=rate)
}

if(model=='shock'){
  rate=NA
  states<-1:Nx ## state=wearing level of the system, 1=optimal; Nx=failure
  if(is.na(alpha[1])){alpha<-c(1,rep(0,Nx-1))}
  P<-matrix(0,Nx,Nx)
  diag(P)[1:(Nx-1)]<-1-p;
  for(i in 1:(Nx-1)){P[i,i+1]<-p}
  P[Nx,1]<-1
  up<-1:(Nx-1)
  down<-Nx
}

if(model=='other'){
  Nx<-dim(P)[1]
  ## states must be sorted to have E={U|D}
  states<-1:Nx ## state=number of shocks that have arrived
  up<-1:n.up
  down<-states[-up] ## up=set of sub-indices indicating up-states
  id<-c(up,down)
  P<-P[id,id]
  if(is.na(alpha[1])){alpha<-c(1,rep(0,Nx-1))}else{alpha<-alpha[id]}
  rate=NA;p<-NA}

### Signal space and emission matrix:
signals<-1:Ny
green<-1:n.green
red<-signals[-green]
if(is.na(sum(M))){
  M<-matrix(0,Nx,Ny)
  M[1,1]<-1
  M[Nx,Ny]<-1
  M[2:(Nx-1),1:(Ny)]<-1/(Ny) # assume that the signals have equal probability
}else{
  id<-c(green,red)
  M<-M[,id]}

if(model=='KooN'){M[down,]<-0;M[down,Ny]<-1}


if(model=='KooN')
{
  hmmR<-list(states=states,up=up,down=down,signals=signals,green=green,red=red,alpha=alpha,P=P,M=M,rate=rate)
  message('System K-out-of-N: ', K,'-out-of-',Nx)
  message('Component failure-rate: ',rate)
  return(hmmR)
}


if(model=='shock')
{
  hmmR<-list(states=states,up=up,down=down,signals=signals,green=green,red=red,alpha=alpha,P=P,M=M)
  message('Shock deteriorating system. Shock probability, p= ',p)
  return(hmmR)
}

if(model=='other')
{hmmR<-list(states=states,up=up,down=down,signals=signals,green=green,red=red,alpha=alpha,P=P,M=M)
message('Repairable System')
return(hmmR)
}
}## end of def.hmmR

##################

sim.hmmR<-function(hmmR,n) #sytem K-out-of-ncom
{
  # hmmR is a list with all the elements of the HMM: P, M, alpha, up, green, etc
  # n= sample size: Y_1, Y_2, ... , Y_n

  P<-hmmR$P
  M<-hmmR$M
  alpha<-hmmR$alpha
  Nx<-length(hmmR$states)   ## total number of states
  Ny<-length(hmmR$signals)   ## total number of signals

  states<-hmmR$states;     ## states are sorted: up are first and down last. Same for signals
  signals<-hmmR$signals;


  Xn<-c();Yn<-c()

  Xn[1]<-states[sample(x=1:Nx,size=1,replace=F,prob=alpha)]
  for(k in 2:n)
  { x.k<-sample(x=1:Nx,size=1,replace=F,prob=P[Xn[k-1],])
  Xn<-c(Xn,states[x.k])}
  for(k in 1:n){
    y.k<-sample(x=1:Ny,size=1,replace=F,prob=M[Xn[k],])
    Yn<-c(Yn,signals[y.k])
  }
  return(list(Xn=Xn,Yn=Yn,P=P,M=M,alpha=alpha,states=states,signals=signals))
}## end of sim.hmmR

##################

fit.hmmR<-function(Y,P0,M0,alpha0,max.iter=50,epsilon=1e-9,Nx,Ny)
{
  n<-length(Y)
  states<-1:Nx;
  signals<-1:Ny


  ## INITIAL VALUES:
  m<-1;
  if(is.na(sum(M0))){
    M<-matrix(0,Nx,Ny)
    M[1,1]<-1
    M[Nx,Ny]<-1
    M[2:(Nx-1),1:(Ny)]<-1/(Ny) # assume that the signals have equal probability
  }else{M<-M0}
  M.all<-list(M)
  if(is.na(sum(P0))){
    P<-matrix(rep(1,Nx)/Nx,Nx,Nx,byrow=T)
  }else{P<-P0}
  P.all<-list(P)
  if(is.na(sum(alpha0))){alpha<-rep(0,Nx);alpha[1]<-1}else{alpha<-alpha0}


  ## BAUM-WELCH ALGORITHM:
  ## FORWARD-BACKWARD:
  # 1. Forward probabilities at the m-th iteration:
  ## Organize the forward probabilities in a matrix Fm of dimensions (n+1)x Nx, such that Fm[k,i]=F^(m)_k(i)
  Fm<-matrix(0,n,Nx)

  for(i in 1:Nx){Fm[1,i]<-alpha[i]*M[i,Y[1]]}

  for(k in 2:n)
  {
    Fm[k,]<-Fm[k-1,]%*%P
    Fm[k,]<-Fm[k,]*M[,Y[k]]
  }
  Fm[Fm<0]<-0;


  # 2. Backward probabilities at the m-th iteration:
  ## Organize the backward probabilities in a matrix Bm of dimensions Nx x (n+1), such that Bm[i,k]=B^(m)_k(i)
  Bm<-matrix(0,Nx,n)
  Bm[,n]<-1 ### for convenience!!
  for(k in (n-1):1)
  {Aux<-M[,Y[k+1]]*Bm[,k+1]
  Bm[,k]<-P%*%Aux
  }
  Bm[Bm<0]<-0;

  ## likelihood:
  Py.m<-sum(Bm[,1])   ### equally:  Py.m<-sum(Fm[n+1,])

  ## E-STEP:
  ## Estimate the probabilities of occupation: P(X_k=i|Y) for all k=1,...,n; and i=1,...,Nx
  ##   P^(m)(X_k=i|Y)=(B^(m)_k(i)*F^(m)_k(i))/sum_j(B^(m)_0(j))
  Px.m<-Bm*t(Fm)#/Py.m ### do not need to divide by Py.m because it will cancel, but we need Py.m as the likelihood value
  # matrix[Nx, n]

  ## Estimate the transition probabilities: P(X_{k-1}=i,X_k=j|Y) for all k=1,..., n; i,j=1,..., Nx
  Pxx.m<-array(NA,dim=c(Nx,Nx,n-1))
  for(k in 1:(n-1))
  {Fk<-matrix(Fm[k,],Nx,Nx,byrow=F)
  Bk<-matrix(Bm[,k+1],Nx,Nx,byrow=T)
  Myk<-matrix(M[,Y[k+1]],Nx,Nx,byrow=T)
  Pxx.m[,,k]<-Fk*P*Myk*Bk}

  #########
  ## M-STEP:
  # 1. new.P:
  # new.P[i,j]<-sum(P.xx[i,j,])/sum(P.x[,i])!!!
  num<-apply(Pxx.m,1:2,sum,na.rm=T)
  den<-rowSums(Px.m[,-n],na.rm=T)
  new.P<-num/den; new.P[den==0]<-0
  new.P[new.P<0]<-0
  new.P[new.P>1]<-1
  new.P<-new.P/rowSums(new.P,na.rm=T)
  # 2. new.M:  #non parametric approach to estimate matrix M!!!!
  Isa<-matrix(0,n,Ny)
  for(a in 1:Ny){Isa[,a]<-as.numeric(Y==a)}
  P.aux<-Px.m%*%Isa
  new.M<-P.aux/rowSums(Px.m,na.rm=T)
  new.M[is.na(new.M)]<-0
  new.M[new.M<0]<-0;new.M[new.M>1]<-1
  new.M<-new.M/rowSums(new.M,na.rm=T)

  # Check the convergence:
  tol.P<-sum((new.P-P)^2,na.rm=T)/sum(P^2,na.rm=T);
  tol.M<-sum((new.M-M)^2,na.rm=T)/sum(M^2,na.rm=T);
  new.tol<-tol.P +tol.M
  v.tol<-c(new.tol);

  if(tol.P<epsilon){new.P<-P} ### no cambio los parametros si ya son muy parecidos
  if(tol.M<epsilon){new.M<-M}
  v.tol.P<-c(tol.P);v.tol.M<-c(tol.M)

  # Keep the estimations:
  M.all[[length(M.all)+1]]<-new.M
  P.all[[length(P.all)+1]]<-new.P
  v.logLik<-c(log(Py.m))

  control<-0
  if(sum(Px.m)==0){control<-1}else{
    while((m<max.iter)&(new.tol>2*epsilon)) ### m= the iteration number
    {
      #update the parameters:
      m<-m+1
      P<-new.P
      M<-new.M

      ## FORWARD-BACKWARD:
      # 1. Forward probabilities at the m-th iteration:
      ## Organize the forward probabilities in a matrix Fm of dimensions (n+1)x Nx, such that Fm[k,i]=F^(m)_k(i)
      Fm<-matrix(0,n,Nx)

      for(i in 1:Nx){Fm[1,i]<-alpha[i]*M[i,Y[1]]}

      for(k in 2:n)
      {
        Fm[k,]<-Fm[k-1,]%*%P
        Fm[k,]<-Fm[k,]*M[,Y[k]]
      }
      Fm[Fm<0]<-0;


      # 2. Backward probabilities at the m-th iteration:
      ## Organize the backward probabilities in a matrix Bm of dimensions Nx x (n+1), such that Bm[i,k]=B^(m)_k(i)
      Bm<-matrix(0,Nx,n)
      Bm[,n]<-1 ### for convenience!!
      for(k in (n-1):1)
      {Aux<-M[,Y[k+1]]*Bm[,k+1]
      Bm[,k]<-P%*%Aux
      }
      Bm[Bm<0]<-0;

      ## likelihood:
      Py.m<-sum(Bm[,1])   ### equally:  Py.m<-sum(Fm[n+1,])

      ## E-STEP:
      ## Estimate the probabilities of occupation: P(X_k=i|Y) for all k=1,...,n; and i=1,...,Nx
      ##   P^(m)(X_k=i|Y)=(B^(m)_k(i)*F^(m)_k(i))/sum_j(B^(m)_0(j))
      Px.m<-Bm*t(Fm)#/Py.m ### do not need to divide by Py.m because it will cancel, but we need Py.m as the likelihood value
      # matrix[Nx, n]

      ## Estimate the transition probabilities: P(X_{k-1}=i,X_k=j|Y) for all k=1,..., n; i,j=1,..., Nx
      Pxx.m<-array(NA,dim=c(Nx,Nx,n-1))
      for(k in 1:(n-1))
      {Fk<-matrix(Fm[k,],Nx,Nx,byrow=F)
      Bk<-matrix(Bm[,k+1],Nx,Nx,byrow=T)
      Myk<-matrix(M[,Y[k+1]],Nx,Nx,byrow=T)
      Pxx.m[,,k]<-Fk*P*Myk*Bk}

      #########
      ## M-STEP:
      # 1. new.P:
      # new.P[i,j]<-sum(P.xx[i,j,])/sum(P.x[,i])!!!
      num<-apply(Pxx.m,1:2,sum,na.rm=T)
      den<-rowSums(Px.m[,-n],na.rm=T)
      new.P<-num/den; new.P[den==0]<-0
      new.P[new.P<0]<-0
      new.P[new.P>1]<-1
      new.P<-new.P/rowSums(new.P,na.rm=T)
      # 2. new.M:  #non parametric approach to estimate matrix M!!!!
      Isa<-matrix(0,n,Ny)
      for(a in 1:Ny){Isa[,a]<-as.numeric(Y==a)}
      P.aux<-Px.m%*%Isa
      new.M<-P.aux/rowSums(Px.m,na.rm=T)
      new.M[is.na(new.M)]<-0
      new.M[new.M<0]<-0;new.M[new.M>1]<-1
      new.M<-new.M/rowSums(new.M,na.rm=T)


      # Check the convergence:
      tol.P<-sum((new.P-P)^2,na.rm=T)/sum(P^2,na.rm=T);
      tol.M<-sum((new.M-M)^2,na.rm=T)/sum(M^2,na.rm=T);
      new.tol<-tol.P +tol.M
      v.tol<-c(v.tol,new.tol)
      if(tol.P<epsilon){new.P<-P}
      if(tol.M<epsilon){new.M<-M}
      v.tol.P<-c(v.tol.P,tol.P);v.tol.M<-c(v.tol.M,tol.M)

      # Keep the estimations:
      M.all[[length(M.all)+1]]<-new.M
      P.all[[length(P.all)+1]]<-new.P
      v.logLik<-c(v.logLik,log(Py.m))
    }
    n.par<-Nx*(Nx-1)+Nx*(Ny-1)
    log.lik<--2*log(Py.m)
    aic<-2*n.par+log.lik
    result<-list(P=round(new.P,5),M=round(new.M,5),n.iter=m,
                 v.tol=v.tol,AIC=aic,v.logLik=v.logLik,P.backward=Bm,P.forward=Fm,
                 P.all=P.all,M.all=M.all,v.tol.P=v.tol.P,v.tol.M=v.tol.M)
  }
  if(control==1){message('No solution')}else{return(result)}


}## end of fit.hmmR

###################

Rcalc.hmmR<-function(hmmR,t=0)
{
  # hmmR = list(states,signals,up,alpha=alpha,P=P,M=M,down,green,red,...)
  # alpha=initial law
  # P  = transition matrix; Nx<-nrow(P)
  # up  = number or subset of Up states
  # M  = emission matrix; Ny<-ncol(M)
  # green= number or subset of good (green) signals

  ## The set states is sorted such that Up states are first, and down states are last.
  ## Same with signals:
  nt<-length(t)
  if(t[nt]==0){res<-1}else{
    t<-t[-1];nt<-nt-1
  n.up<-length(hmmR$up)
  n.green<-length(hmmR$green)

  alpha<-hmmR$alpha
  alpha.1<-alpha[1:n.up]

  P<-hmmR$P;
  P.U<-matrix(P[1:n.up,1:n.up],nrow=n.up,ncol=n.up) ### asumme that the states are sorted states = U + D

  M<-hmmR$M;
  M.1<-matrix(M[1:n.up,1:n.green],nrow=n.up,ncol=n.green) ### assume that the signals are sorted signals = green + red

  ### Let us build matrix Psi.U ### transition matrix of the 2-dim MC (X,Y)
  #                             ### restricted to the up states and green signals.
  Psi.1<-c()
  for(s in 1:n.green)
  {Psi.1<-cbind(Psi.1,P.U%*%diag(M.1[,s],nrow=n.up,ncol=n.up))}

  s<-1;Psi.U<-Psi.1
  while(s<n.green){s<-s+1;Psi.U<-rbind(Psi.U,Psi.1)}

  #### The reliability calculation:
  ### The initial law of the 2-dim process:
  ### Alpha[i,y]=P(X0=i;Y0=y)=alpha[i]*M[i,y]
  Alpha.1<-t(as.vector(M.1*alpha.1))

  Reliab.t<-c(1);n<-0;reliab.n<-Alpha.1
  while(n<nt){
    n<-n+1;reliab.n<-reliab.n%*%Psi.U
    Reliab.t<-c(Reliab.t,rowSums(reliab.n))
  }
  Reliab.t[Reliab.t<0]<-0
  Reliab.t[Reliab.t>1]<-1
  res<-c(Reliab.t)
  if(nt==1){return(res[nt+1])}else{return(res)}}

}## end of Rcalc.hmmR

###################

cost.cspc<-function(prob,hmmR,n.up1,cost.C,cost.P,t.max)
{
  P<-hmmR$P
  states<-hmmR$states;Nx<-length(states)
  up<-hmmR$up
  down<-states[-up]
  n.down<-length(down)
  cost.C<-rep(cost.C,n.down)

  n.up<-length(up)
  cost.P<-rep(cost.P,n.up)

  alpha<-hmmR$alpha

# Algorithm CSPC: Critical State Probability Criterion
  P.U1<-as.matrix(P[1:n.up1,1:n.up1]) ### asumme that the states are sorted states = U + D
  P.U2<-as.matrix(P[1:n.up1,(n.up1+1):n.up])
  alpha.11<-as.vector(alpha[1:n.up1])

  #The distribution law of first entry time in U2:
       potU<-c(alpha.11);i<-0;f.U2<-c();F.U2<-c()
       while(i<t.max)
       {i<-i+1
       potU<-potU%*%P.U1
       f.U2<-c(f.U2,sum(potU%*%P.U2,na.rm=T))}
       F.U2<-cumsum(f.U2)
       F.U2<-F.U2/F.U2[t.max]

    #t.U2=min{n:F.U2(t.U2)>prob}
    t.U2<-0;p0<-0
    while(p0<prob){t.U2<-t.U2+1; p0<-F.U2[t.U2]}

  # the system is preventively maintained at t=t.U2,
  # the system is returned to state a state in up, with probability alpha.1
  # the system will then have another PM inspection every t.U2 steps.
  # for time.pm<-c(t.U2,2*t.U2,3*t.U2,...)

    #1. Calculate the total cost incurred during a complete inspection period: (0,tU2]
    n.pm<-t.max%/%t.U2
    r.pm<-t.max%%t.U2
    n<-1;E.pm<-0;E.cm<-0;last.cost<-0;nr<-1
    potP<-diag(1,Nx,Nx)
    while(n<=t.U2){
      ## the occupation probabilities
      potP<-potP%*%P; pn<-alpha%*%potP

      #E[C(n)]= P(Xn= D)*C_{cm} + I_{Nc(q)=n}*P(Xn=U2)*C_{pm}
      e.cm<-pn[down]%*%cost.C ## expected corrective maintenance cost
      is.pm<-sum(n==t.U2)
      e.pm<-is.pm*(pn[up]%*%cost.P) ## expected preventive maintenance cost

      E.cm<-E.cm+e.cm;
      E.pm<-E.pm+e.pm
      if(r.pm>0){if(nr<r.pm){last.cost<-last.cost+e.cm;nr<-nr+1}}
      n<-n+1
    }

    if(r.pm==0){last.cost<-0}
    total.cost<-n.pm*(E.cm+E.pm)+last.cost

res<-list(time.insp=t.U2,n.insp=n.pm,total.cost=total.cost)
return(res)
}## end of cost.cspc


##################

cost.wspc<-function(prob,hmmR,n.up1,n.green1,cost.C,cost.P,t.max)
{
  P<-hmmR$P
  alpha<-hmmR$alpha
  states<-hmmR$states;Nx<-length(states)
  up<-hmmR$up;n.up<-length(up)
  down<-states[-up];n.down<-length(down)

  M<-hmmR$M
  signals<-hmmR$signals;Ny<-length(signals)
  green<-hmmR$green;n.green<-length(green)
  red<-signals[-green]

  ## vector of cost:
  Down<-c()
  for(s in 1:n.green){Down<-c(Down,Nx*(s-1)+down)}
  Down<-c(Down,(Nx*n.green+1):(Nx*Ny))
  n.Down<-length(Down)
  cost.C<-rep(cost.C,n.Down)
  States<-1:(Nx*Ny)
  Up<-States[-Down]
  n.Up<-length(Up)
  cost.P<-rep(cost.P,n.Up)

  ## Algorithm WSPC: Warning Signal Probability Criterion
  ## build transition matrix (Psi) of the 2-dim process with
  ## states: {{1:Nx}x{1}, {1:Nx}x{2},..., {1:Nx}x{Ny}}

  Psi0<-c()
  for(s in 1:Ny){
    aux<-P%*%diag(M[,s])
    Psi0<-cbind(Psi0,aux)
  }
  s<-1;Psi<-Psi0
  while(s<Ny){s<-s+1;Psi<-rbind(Psi,Psi0)}

  #the submatrices:
  Psi.G1<-Psi[1:(Nx*n.green1),1:(Nx*n.green1)]
  Psi.G2<-Psi[1:(Nx*n.green1),(Nx*n.green1+1):(Nx*n.green)]


  ## initial law of the 2-dim process: Alpha
  Alpha<-as.vector(alpha*M)
  Alpha.1<-Alpha[1:(Nx*n.green1)] ## restricted to non-critical signals

  #The distribution law of time to first entrance in G2 (set of critical signals):
  potG<-c(Alpha.1);i<-0;f.G2<-c()
  while(i<t.max)
  { i<-i+1
    potG<-potG%*%Psi.G1
    f.G2<-c(f.G2,sum(potG%*%Psi.G2,na.rm=T))}
  F.G2<-cumsum(f.G2)
  F.G2<- F.G2/F.G2[t.max]

  #tG2=min{n:F.G2(n)>prob}
  t.G2<-0;p0<-0
  while(p0<prob){t.G2<-t.G2+1; p0<-F.G2[t.G2]}

  # the system is preventively maintained at t=t.G2,
  # the system is returned to an up state
  # the system will then have another PM inspetion every t.G2 steps.
  # Inspection times: time.pm<-c(t.G2,2*t.G2,3*t.G2,...)

  # Calculate the total cost incurred during a complete inspection period: (0,t.G2)
  n.pm<-t.max%/%t.G2
  r.pm<-t.max%%t.G2
  n<-nr<-1
  E.pm<-E.cm<-last.cost<-0
  potPsi<-diag(1,Nx*Ny,Nx*Ny)
  while(n<=t.G2){
    ## the occupation probabilities
    potPsi<-potPsi%*%Psi;psi.n<-Alpha%*%potPsi
    e.cm<-psi.n[Down]%*%cost.C ## expected corrective maintenance cost
    is.pm<-sum(n==t.G2)
    e.pm<-is.pm*(psi.n[Up]%*%cost.P) ## expected preventive maintenance cost
    E.cm<-E.cm+e.cm
    E.pm<-E.pm+e.pm
    if(r.pm>0){if(nr<r.pm){last.cost<-last.cost+e.cm;nr<-nr+1}}
    n<-n+1
   }

  if(r.pm==0){last.cost<-0}
  total.cost<-n.pm*(E.cm+E.pm)+last.cost

  res<-list(time.insp=t.G2,n.insp=n.pm,total.cost=total.cost)

  return(res)
}## end of cost.wspc




