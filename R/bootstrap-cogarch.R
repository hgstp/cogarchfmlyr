psi.2 = function(eta,phi) -2*eta+2*phi+3*phi^2

psi.4=function(eta,phi) -4*eta+105*phi^4 + 60*phi^3+18*phi^2+4*phi

-2*eta+(105*phi^3+60*phi^2+15*phi+2)*phi

beta.sim=.1
eta.sim = .1
phi.sim= .05


psi.2(eta.sim,phi.sim)
psi.4(eta.sim,phi.sim)


G.sim.cp=cogarch_sim(t=0:201000,beta=beta.sim,eta=eta.sim, phi=phi.sim)


#G.sim=cogarch_sim(t=0:5000,beta=beta.sim,eta=eta.sim,phi=phi.sim,Lp="vg",sigma=1,nu=1,theta=0)



G.cp = as.matrix(prevTick(1:200000,G.sim.cp[,5],G.sim.cp[,1]))

#time = seq(101,500001, by=100)

#G= as.matrix(G.sim[time,1])

#theta = (betha,eta,phi)

(theta.hat =est_cogarch(G.cp)[1:3])

#(theta.hat =est_cogarch(G)[1:3])

#n= dim(G)[1]

b= 4751

overlap = floor(b/4)

samples = floor((n-b)/overlap)+1

theta.hat.boot = matrix(rep(0,samples*3),nrow=samples,ncol=3)

for( k in 1:samples){
  theta.hat.boot[k,]=est_cogarch( as.matrix(G.cp[((k-1)*overlap+1):(((k-1)*overlap)+b)]))[1:3]
}


#for( k in 1:4000) theta.hat.boot[k,]=est_cogarch( as.matrix(G[k:(k+b),1]))[1:3]

theta.hat

length(which(theta.hat.boot[,1]==0))

theta.hat.boot = theta.hat.boot[-which(theta.hat.boot[,1]==0),]


par(mfrow=c(1,3))
hist(sqrt(b)*(theta.hat.boot[,1]-theta.hat[1]),freq=F, main=beta.sim,xlab="beta")

hist(sqrt(b)*(theta.hat.boot[,2]-theta.hat[2]),freq=F,main=eta.sim,xlab="eta")
hist(sqrt(b)*(theta.hat.boot[,3]-theta.hat[3]),freq=F,main=phi.sim,xlab="phi")


par(mfrow=c(1,3))
hist(sqrt(b)*(theta.hat.boot[,1]-beta.sim),freq=F, main=beta.sim,xlab="beta")

hist(sqrt(b)*(theta.hat.boot[,2]-eta.sim),freq=F,main=eta.sim,xlab="eta")
hist(sqrt(b)*(theta.hat.boot[,3]-phi.sim),freq=F,main=phi.sim,xlab="phi")


# Bias

bias.beta = sqrt(b)/sqrt(200000)*(mean(theta.hat.boot[,1])-theta.hat[1])
bias.eta = sqrt(b)/sqrt(200000)*(mean(theta.hat.boot[,2])-theta.hat[2])
bias.phi = sqrt(b)/sqrt(200000)*(mean(theta.hat.boot[,3])-theta.hat[3])


theta.hat + c(bias.beta,bias.eta,bias.phi)

theta.hat
