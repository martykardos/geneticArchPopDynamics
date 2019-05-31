
library(quantPop)
###################################################################
# Individual-based simulations with variable initial allele
# frequency among quantitative trait loci. See function help files
# for parameter definitions 
###################################################################
reps <- 20      # number of simulation replicates to run

### make some objects to store simulation output
NMat <- NULL  # population size   
phenMat <- NULL # mean phenotype per generation

for(j in 1:reps){
  logQuantIndiv_varP0(N0=500,K=1000,t=80,lambda=1.5,f=4,sexType="herm",nLoci=100,phen_0=100,phen_opt=110,
                      fit_sd=6,h2_0=0.6,Vp_0=10,beta1=0.5,beta2=0.5,freqBounds=c(0.05,0.95),nMajor=1,propMajor=1/100,
                      freqMajor=runif(1,min=0.05,max=0.95),immRate=0)
  NMat <- rbind(NMat,NVec)
  phenMat <- rbind(phenMat,phenVec)
  print(j)
}

# Plot the simulation results
par(mfrow=c(2,1))

# population size
plot(c(0,ncol(NMat)),c(0,1300),xlab="Generation",ylab="Population size",type="n")
for(i in 1:nrow(NMat)){
  lines(1:ncol(NMat),NMat[i,],col="darkblue")
}

# heritability
plot(c(0,ncol(phenMat)),c(98,115),xlab="Generation",ylab="Heritability",type="n")
for(i in 1:nrow(phenMat)){
  lines(1:ncol(phenMat),phenMat[i,],col="darkblue")
}

############################################################################
# individual-based simulations with a stochastic, temporal, linear increase
# in the phenotypic optimum
############################################################################
reps <- 10       # number of simulation replicates
NMat <- NULL    # store population sizes
phenMat <- NULL # store mean phenotype
gensToOpt <- 20
for(j in 1:reps){
  logQuantIndiv_varP0_stocPhenOpt(N0=500,K=1000,t=80,lambda=1.5,f=4,sexType="herm",nLoci=100,phen_0=100,phen_opt=110,timeToOpt=gensToOpt,optSd=2,  
                      fit_sd=6,h2_0=0.6,Vp_0=10,beta1=0.5,beta2=0.5,freqBounds=c(0.05,0.95),nMajor=1,propMajor=1/100,
                      freqMajor=runif(1,min=0.05,max=0.95),immRate=0)
  NMat <- rbind(NMat,NVec)
  phenMat <- rbind(phenMat,phenVec)
  print(j)
}


# Plot the simulation results
par(mfrow=c(2,1))

# population size
plot(c(0,ncol(NMat)),c(0,1300),xlab="Generation",ylab="Population size",type="n")
for(i in 1:nrow(NMat)){
  lines(1:ncol(NMat),NMat[i,],col="darkblue")
}

# heritability
plot(c(0,ncol(phenMat)),c(95,115),xlab="Generation",ylab="Heritability",type="n")
for(i in 1:nrow(phenMat)){
  lines(1:ncol(phenMat),phenMat[i,],col="darkblue")
}


