
library (quantPop)
###############################################################################
# plot major and polygenic architecture effects
# across a range of starting allele frequencies when 
# there are one, two, and 100 loci and starting heritability = 0.6
###############################################################################
startN <- 0.5                    # starting population size (density)
carryCap <- 1                    # carrying capacity
t <- 80                          # number of generations to analyze the population dynamics
intGrowth <- 1.5                 # lambda for a perfectly adapted population at N -> 0.
startPheno <- 100                # beginning phenotype mean
optPheno <- 110                  # optimum phenotype
sdFitFun <- 6                    # standard deviation of the Gaussian fitness function
h2Start <- 0.6                   # starting narrow sense heritability
phenVar <- 10                    # starting phenotype variance
startFreq <- c(0.1,0.3,0.5,0.7,0.9)     # starting allele frequencies of interest


######### 1 major locus
major1NMat <- NULL
major1PhenMat <- NULL
major1H2Mat <- NULL
major1PMat <- NULL
for(i in 1:length(startFreq)){
  logQuant(N0=startN,K=carryCap,t=t,lambda=intGrowth,nLoci=1,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,Vp_0=phenVar,p0=startFreq[i])
  major1NMat <- rbind(major1NMat,NVec)
  major1PhenMat <- rbind(major1PhenMat,phenVec)
  major1H2Mat <- rbind(major1H2Mat,h2Vec)
  major1PMat <- rbind(major1PMat,freqVec)
}

######### 2 major loci
major2NMat <- NULL
major2PhenMat <- NULL
major2H2Mat <- NULL
major2PMat <- NULL
for(i in 1:length(startFreq)){
  logQuant(N0=startN,K=carryCap,t=t,lambda=intGrowth,nLoci=2,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,Vp_0=phenVar,p0=startFreq[i])
  major2NMat <- rbind(major2NMat,NVec)
  major2PhenMat <- rbind(major2PhenMat,phenVec)
  major2H2Mat <- rbind(major2H2Mat,h2Vec)
  major2PMat <- rbind(major2PMat,freqVec)
}

###### POLYGENIC
polyNMat <- NULL
polyPhenMat <- NULL
polyH2Mat <- NULL
polyPMat <- NULL
for(i in 1:length(startFreq)){
  logQuant(N0=startN,K=carryCap,t=t,lambda=intGrowth,nLoci=100,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,Vp_0=phenVar,p0=startFreq[i])
  polyNMat <- rbind(polyNMat,NVec)
  polyPhenMat <- rbind(polyPhenMat,phenVec)
  polyH2Mat <- rbind(polyH2Mat,h2Vec)
  polyPMat <- rbind(polyPMat,freqVec)
}


########################################
#################### Make the figure
########################################

# population size

# one major locus
par(mfrow=c(3,3),xpd=TRUE,mar=c(4,5,2,2))
plot(1:t,1:t,ylim=c(0,carryCap),type="n",cex.axis=0.8, main ="1 major locus",cex.main=1.4,
     cex.lab=1.5,ylab=expression(paste(italic(""*N*"")," / ",italic(""*K*""))),xlab="")
text(x=-25,y=1.1,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:nrow(major1NMat)){
  lines(1:t,major1NMat[i,1:t],lty=lineVec[i])
}

# two major loci
plot(1:t,1:t,ylim=c(0,carryCap),type="n",cex.axis=0.8,main ="2 major loci",cex.main=1.4,
     cex.lab=1.5,ylab="",xlab="")
text(x=-25,y=1.1,labels="B",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:nrow(major2NMat)){
  lines(1:t,major2NMat[i,1:t],lty=lineVec[i])
}

# polygenic
plot(1:t,1:t,ylim=c(0,carryCap),type="n",cex.axis=0.8,main ="100 small-effect loci",cex.main=1.4,
     cex.lab=1.5,ylab="",xlab="")
text(x=-25,y=1.1,labels="C",cex=2)
for(i in 1:nrow(polyNMat)){
  lines(1:t,polyNMat[i,1:t],lty=lineVec[i])
}
legend(x=20,y=0,lty=lineVec,legend=c(expression(paste(italic(""*p*"")[0]," = 0.1",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.3",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.5",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.7",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.9",sep=""))),
       bty="n",xjust=FALSE,yjust=FALSE,cex=0.9)

############# heritability

# one major locus
plot(1:t,1:t,ylim=c(0,0.8),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab=expression(italic(""*h*"")^2),xlab="")
for(i in 1:nrow(major1H2Mat)){
  lines(1:t,major1H2Mat[i,1:t],lty=lineVec[i])
}

# two major loci
plot(1:t,1:t,ylim=c(0,0.8),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="")
for(i in 1:nrow(major1H2Mat)){
  lines(1:t,major2H2Mat[i,1:t],lty=lineVec[i])
}

# 100 loci
plot(1:t,1:t,ylim=c(0,0.8),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="")
for(i in 1:nrow(polyH2Mat)){
  lines(1:t,polyH2Mat[i,1:t],lty=lineVec[i])
}

############ phenotype
# one major locus
plot(1:t,1:t,ylim=c(100,112),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="Mean phenotype",xlab="Generation")
for(i in 1:nrow(major1PhenMat)){
  lines(1:t,major1PhenMat[i,1:t],lty=lineVec[i])
}

# two major loci
plot(1:t,1:t,ylim=c(100,112),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="Generation")
for(i in 1:nrow(major2PhenMat)){
  lines(1:t,major2PhenMat[i,1:t],lty=lineVec[i])
}

# 100 loci
plot(1:t,1:t,ylim=c(100,112),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="Generation")
for(i in 1:nrow(polyPhenMat)){
  lines(1:t,polyPhenMat[i,1:t],lty=lineVec[i])
}



###############################################################################
# plot major and polygenic architecture effects
# across a range of starting frequencies  when there are one, two, and 100 loci
# WITH QUASI GENOTYPES TO ACCOUNT FOR LD AND NON-NORMAL PHENOTYPE DISTRIBUTION
###############################################################################

startN <- 0.5
carryCap <- 1
runTime <- 80
intGrowth <- 1.5                 
startPheno <- 100
optPheno <- 110
sdFitFun <- 6
h2Start <- 0.6
phenVar <- 10
startFreq <- c(0.1,0.25,0.5,0.75,0.9)


######### 1 major locus
major1NMat <- NULL
major1PhenMat <- NULL
major1H2Mat <- NULL
major1PMat <- NULL
for(i in 1:length(startFreq)){
  logQuantPseudoGenos(N0=startN,K=carryCap,t=runTime,lambda=intGrowth,nLoci=1,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,
           Vp_0=phenVar,p0=startFreq[i],pseudoN=50000)
  major1NMat <- rbind(major1NMat,NVec)
  major1PhenMat <- rbind(major1PhenMat,phenVec)
  major1H2Mat <- rbind(major1H2Mat,h2Vec)
  major1PMat <- rbind(major1PMat,freqVec)
}



######### 2 major loci
major2NMat <- NULL
major2PhenMat <- NULL
major2H2Mat <- NULL
major2PMat <- NULL
for(i in 1:length(startFreq)){
  logQuantPseudoGenos(N0=startN,K=carryCap,t=runTime,lambda=intGrowth,nLoci=2,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,
           Vp_0=phenVar,p0=startFreq[i],pseudoN=50000)
  major2NMat <- rbind(major2NMat,NVec)
  major2PhenMat <- rbind(major2PhenMat,phenVec)
  major2H2Mat <- rbind(major2H2Mat,h2Vec)
  major2PMat <- rbind(major2PMat,freqVec)
}

###### POLYGENIC   (THIS WILL TAKE A WHILE BECAUSE OF THE LARGE NUMBER OF LOCI AND PSEUDO-INDIVIDUALS)
polyNMat <- NULL
polyPhenMat <- NULL
polyH2Mat <- NULL
polyPMat <- NULL
for(i in 1:length(startFreq)){
  logQuantPseudoGenos(N0=startN,K=carryCap,t=runTime,lambda=intGrowth,nLoci=100,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,
           Vp_0=phenVar,p0=startFreq[i],pseudoN = 50000)
  polyNMat <- rbind(polyNMat,NVec)
  polyPhenMat <- rbind(polyPhenMat,phenVec)
  polyH2Mat <- rbind(polyH2Mat,h2Vec)
  polyPMat <- rbind(polyPMat,freqVec)
}


########################################
#################### Make the figure
########################################

# population size

# one major locus
par(mfrow=c(3,3),xpd=TRUE,mar=c(4,5,2,2))
plot(1:runTime,1:runTime,ylim=c(0,carryCap),type="n",cex.axis=0.8, main ="1 major locus",cex.main=1.4,
     cex.lab=1.5,ylab=expression(paste(italic(""*N*"")," / ",italic(""*K*""))),xlab="")
text(x=-35,y=1.1,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:nrow(major1NMat)){
  lines(1:runTime,major1NMat[i,1:runTime],lty=lineVec[i])
}

# two major loci
plot(1:t,1:runTime,ylim=c(0,carryCap),type="n",cex.axis=0.8,main ="2 major loci",cex.main=1.4,
     cex.lab=1.5,ylab="",xlab="")
text(x=-35,y=1.1,labels="B",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:nrow(major2NMat)){
  lines(1:t,major2NMat[i,1:runTime],lty=lineVec[i])
}

plot(1:runTime,1:runTime,ylim=c(0,carryCap),type="n",cex.axis=0.8,main ="100 small-effect loci",cex.main=1.4,
     cex.lab=1.5,ylab="",xlab="")
text(x=-35,y=1.1,labels="C",cex=2)
for(i in 1:nrow(polyNMat)){
  lines(1:t,polyNMat[i,1:runTime],lty=lineVec[i])
}
legend(x=20,y=0,lty=lineVec,legend=c(expression(paste(italic(""*p*"")[0]," = 0.1",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.25",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.5",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.75",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.9",sep=""))),
       bty="n",xjust=FALSE,yjust=FALSE,cex=0.9)

############# heritability

# one major locus
plot(1:runTime,1:runTime,ylim=c(0,1),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab=expression(italic(""*h*"")^2),xlab="")
for(i in 1:nrow(major1H2Mat)){
  lines(1:runTime,major1H2Mat[i,1:runTime],lty=lineVec[i])
}

# two major loci
plot(1:t,1:t,ylim=c(0,1),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="")
for(i in 1:nrow(polyH2Mat)){
  lines(1:runTime,major2H2Mat[i,1:runTime],lty=lineVec[i])
}

# 100 loci
plot(1:runTime,1:runTime,ylim=c(0,1),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="")
for(i in 1:nrow(polyH2Mat)){
  lines(1:runTime,polyH2Mat[i,1:runTime],lty=lineVec[i])
}


############ phenotype

# one major locus
plot(1:runTime,1:runTime,ylim=c(100,112),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="Mean phenotype",xlab="Generation")
for(i in 1:nrow(major1PhenMat)){
  lines(1:runTime,major1PhenMat[i,1:runTime],lty=lineVec[i])
}

# two major loci
plot(1:runTime,1:runTime,ylim=c(100,112),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="Generation")
for(i in 1:nrow(majorPhenMat)){
  lines(1:runTime,major2PhenMat[i,1:runTime],lty=lineVec[i])
}

# 100 loci
plot(1:runTime,1:runTime,ylim=c(100,112),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="Generation")
for(i in 1:nrow(polyPhenMat)){
  lines(1:runTime,polyPhenMat[i,1:runTime],lty=lineVec[i])
}


