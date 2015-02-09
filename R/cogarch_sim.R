cogarch_sim<-function(t=0:10,beta=1,eta=0.05,phi=0.03,Lp="cp",rate=1,distribution="normal",mean=0,var=1,sigma=1,nu=0.5,theta=1,gs=0.01)
{# Simulation of a COGARCH(1,1) process driven by a compound Poisson or a Variance-Gamma Lévy process
  #Input: t - fixed time grid
  #       beta,eta,phi - Input parameters for the cogarch process, have to be positive  
  #       Lp - driving Lévy process: either "cp"=compound Poisson or ""vg"= Variance-Gamma 
  #       rate - intensity of the compound Poisson Levy process
  #       distribution - jump size distribution, e.g "normal"
  #       mean,var - mean and variance of specified jump distribution
  #       sigma, nu, theta - parameters for Variance-Gamma process (for E(L_1)=0 and E(L_1^2)=1 choose theta=0 and sigma=1)
  #       gs - time grid size
  #Output: G - COGARCH(1,1) process
  #        vol - volatility process
  #        Lt - driving Lévy process
  #        delta_Lt - jumps of Lévy process
  #        randomjumptimes - random jumptimes
  
  #use the Euler-Maruyama scheme in order to get a numerical solution for a SDE
  #See the Book: Peter E. Kloeden: "Numerical solution of stochastic differential equations"
  t<-sort(t)
  n_t<-length(t)
  
  #find starting values
  if(Lp=="cp"){
    index<-1
  }else if(Lp=="vg"){
    index<-2
  }else{
    stop("No valid Levy process specified! Select 'cp' or 'vg'.")
  }
  switch(index,
{output_cp<-compoundPoisson(0:10,rate,distribution,mean,var) #calls the function compoundPoisson with output c(randomjumptimes,randomjumpsizes,Lt)
 rt<-output_cp[,1] #randomjumptimes
 #delta_rt<-rt[2:N]-rt[1:(N-1)] #randomjumpintervals, delta_ti
 delta_rt<-output_cp[,4]#randomjumpintervals, delta_ti
 delta_Lt<-output_cp[,2] #randomjumpsizes, delta_Lti
 Lt<-output_cp[,3] #L_t
},
{output_vg<-vargamma(0:10,sigma,nu,theta,gs) #calls the function varGamma with output timesequence and V
 Lt<-output_vg[,2]
 rt<-output_vg[,1]
 nv<-length(Lt)
 nt<-length(rt)
 
 #increments of the process
 delta_Lt<-Lt[2:nv]-Lt[1:(nv-1)]
 delta_rt<-rt[2:nt]-rt[1:(nt-1)]
}
  )
  
  
  n<-length(delta_Lt)
  
  delta_Lt_0<-c(0,delta_Lt[1:(n-1)]) #used for sigma_t-, i.e. the sigma_t without the jump at time t!
  
  #start value volatility, for state process set Y=0 in the beginning and therefore for the observation process: V = beta/eta + phi*Y = beta/eta
  voli<-beta/eta
  
  #state space process
  Yi<-0
  g<-0
  
  for (k in 1:n){
    
    Yi<-Yi-eta*Yi*delta_rt[k]+voli*delta_Lt_0[k]^2  # delta_Lt hat null im ersten Vektoreintrag, daher erhalte sigma_t-
    voli<-beta/eta+phi*Yi 
    g<-g+sqrt(voli)*delta_Lt[k] 
    
  }
  
  #-----------------------------------------------
  #obtain values via Euler approximation
  
  switch(index,
{output_cp_new<-compoundPoisson(t,rate,distribution,mean,var) #calls function compoundPoisson
 
 rt_new<-output_cp_new[,1] #randomjumptimes
 #NN<-length(rt_new)
 #delta_rt_new<-rt_new[2:NN]-rt_new[1:(NN-1)] #randomjumpintervals, delta_ti
 delta_rt_new<-output_cp_new[,4]
 delta_Lt_new<-output_cp_new[,2] #randomjumpsizes, delta_Lti
 Lt_new<-output_cp_new[,3] #L_t
},
{output_vg_new<-vargamma(t,sigma,nu,theta,gs) #calls the function varGamma with output timesequence and V
 
 Lt_new<-output_vg_new[,2]
 rt_new<-output_vg_new[,1]
 nvn<-length(Lt_new)
 ntn<-length(rt_new)
 
 #increments of the process
 delta_Lt_new<-Lt_new[2:nvn]-Lt_new[1:(nvn-1)]
 delta_rt_new<-rt_new[2:ntn]-rt_new[1:(ntn-1)]
 
}
  )
  
  
  nn<-length(delta_Lt_new)
  #add again 0 as first entry of delta_Lt_new in order to get sigma_t- !
  delta_Lt_new_0<-c(0,delta_Lt_new[1:(nn-1)]) #length nn
  
  #volatility
  vol<-vector(length=(nn))
  vol[1]<-voli
  
  #state space
  Y<-vector(length=(nn))
  Y[1]<-Yi
  
  for (k in 1:nn){
    
    Y[k+1]<-Y[k]-eta*Y[k]*delta_rt_new[k]+vol[k]*delta_Lt_new_0[k]^2
    vol[k+1]<-beta/eta+phi*Y[k]
  }
  
  delta_G<-vector(length=length(vol))  #length(vol)=nn+1 mit nn=length(delta_Lt_new)
  delta_G[1]<-g
  s<-sqrt(vol[1:(length(vol)-1)]) #length=nn
  delta_G[2:length(vol)]<-s*delta_Lt_new[1:nn] #length=nn
  G<-cumsum(delta_G) #length=nn+1
  
  switch(index,
{output<-cbind(G,vol,c(0,Lt_new),c(0,delta_Lt_new),c(0,rt_new))
 colnames(output)<-c("G","vol","Lt","delta_Lt","randomjumptimes")
},
{output<-cbind(G,vol,Lt_new,c(0,delta_Lt_new),rt_new) 
 colnames(output)<-c("G","vol","Lt","delta_Lt","randomjumptimes")
}
  )
  return(output)
  
}