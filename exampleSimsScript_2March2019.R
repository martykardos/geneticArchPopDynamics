####### R script to simulate populations with approximate coral and large mammal life histories,
####### with a polygenic selected trait, and also with the selected trait having a large-effect locus.
####### The code below will generate a figure consistent with Figure 1 in our paper.

setwd("~/Desktop/work/rapidAdaptation/MS/sampleScripts")    # directory to save the output
reps <- 20               # number of repetitions to simulate

#######**********************
####### LARGE MAMMAL
#######**********************
####### simulate selection on a polygenic trait and population dynamics in closed populations with large mammal life history
####### for 100 generations
library(pedR)
for (i in 1:reps){
  pedSimHardSelecFast(chrNum=10,chrLengs=rep(50,10),physLengs=rep(100,10),map=NULL,simLociNum=200,beta1=0.5,
                      beta2=0.5,burnin=1,gens=100,lastGenSelec = 100,genoImport=NULL,QTLFreqs=c(0.01,0.99),h2=0.6,numLgQTL=1,
                      numSmQTL=99,propMajor = 1/100,VP=10,phenStart = 100,phenFit = 110,survProbs = 0.8,survProbSlope=0.4/10,carryCap = 300,
                      popSize = 150,famSize = 4)
  write.table(qtlInfo, file=paste("qtlInfo_noMajorGene_100Gens_mammal",i,sep=""),quote=FALSE,row.names=FALSE)
  write.table(quantTraitDat, file=paste("quantDat_noMajorGene_100Gens_mammal",i,sep=""),quote=FALSE,row.names=FALSE)
  write.table(pedObject, file=paste("pedObject_noMajorGene_100Gens_mammal",i,sep=""),quote=FALSE,row.names=FALSE)
}

####### simulate selection on a trait with a large-effect locus and population dynamics in closed populations with large mammal life history
####### for 100 generations
library(pedR)
for (i in 1:reps){
  pedSimHardSelecFast(chrNum=10,chrLengs=rep(50,10),physLengs=rep(100,10),map=NULL,simLociNum=200,beta1=0.5,
                      beta2=0.5,burnin=1,gens=100,lastGenSelec = 100,genoImport=NULL,QTLFreqs=c(0.01,0.99),h2=0.6,numLgQTL=1,
                      numSmQTL=99,propMajor = 0.9,VP=10,phenStart = 100,phenFit = 110,survProbs = 0.8,survProbSlope=0.4/10,carryCap = 300,
                      popSize = 150,famSize = 4)
  write.table(qtlInfo, file=paste("qtlInfo_majorGene_100Gens_mammal",i,sep=""),quote=FALSE,row.names=FALSE)
  write.table(quantTraitDat, file=paste("quantDat_majorGene_100Gens_mammal",i,sep=""),quote=FALSE,row.names=FALSE)
  write.table(pedObject, file=paste("pedObject_majorGene_100Gens_mammal",i,sep=""),quote=FALSE,row.names=FALSE)
}

#######**********************
####### CORAL
#######**********************
####### simulate selection on a polygenic trait and population dynamics in closed populations with a coral life history
####### for 100 generations

for (i in 1:reps){
  pedSimHardSelecFast(chrNum=10,chrLengs=rep(50,10),physLengs=rep(100,10),map=NULL,simLociNum=200,beta1=0.5,
                      beta2=0.5,burnin=1,gens=100,lastGenSelec = 100,genoImport=NULL,QTLFreqs=c(0.01,0.99),h2=0.6,numLgQTL=1,
                      numSmQTL=99,propMajor = 1/100,VP=10,phenStart = 100,phenFit = 110,survProbs = 0.1,survProbSlope=0.04/10,carryCap = 6000,
                      popSize = 3000,famSize = 26.2)
  write.table(qtlInfo, file=paste("qtlInfo_noMajorGene_100Gens_coral",i,sep=""),quote=FALSE,row.names=FALSE)
  write.table(quantTraitDat, file=paste("quantDat_noMajorGene_100Gens_coral",i,sep=""),quote=FALSE,row.names=FALSE)
  write.table(pedObject, file=paste("pedObject_noMajorGene_100Gens_coral",i,sep=""),quote=FALSE,row.names=FALSE)
}

####### simulate selection on a trait with a large-effect locus and population dynamics in closed populations with a coral life history
####### for 100 generations

for (i in 1:reps){
  pedSimHardSelecFast(chrNum=10,chrLengs=rep(50,10),physLengs=rep(100,10),map=NULL,simLociNum=200,beta1=0.5,
                      beta2=0.5,burnin=1,gens=100,lastGenSelec = 100,genoImport=NULL,QTLFreqs=c(0.01,0.99),h2=0.6,numLgQTL=1,
                      numSmQTL=99,propMajor = 0.9,VP=10,phenStart = 100,phenFit = 110,survProbs = 0.1,survProbSlope=0.04/10,carryCap = 6000,
                      popSize = 3000,famSize = 26.2)
  write.table(qtlInfo, file=paste("qtlInfo_majorGene_100Gens_coral",i,sep=""),quote=FALSE,row.names=FALSE)
  write.table(quantTraitDat, file=paste("quantDat_majorGene_100Gens_coral",i,sep=""),quote=FALSE,row.names=FALSE)
  write.table(pedObject, file=paste("pedObject_majorGene_100Gens_coral",i,sep=""),quote=FALSE,row.names=FALSE)
}


####################################################################################
# collect population and phenotype sizes through time for each simulated population
####################################################################################

##### first corals
gens <- 1:100
polyObjects <- paste("pedObject_noMajorGene_100Gens_coral",1:reps,sep="")
oligObjects <- paste("pedObject_majorGene_100Gens_coral",1:reps,sep="")
poly1Phen_coral <- NULL
poly1Size_coral <- NULL
olig1Phen_coral <- NULL
olig1Size_coral <- NULL

for(i in 1:reps){
    polyDat <- read.table(polyObjects[i],header=TRUE)
    poly1SizeVec <- rep(NA,length(gens))
    poly1PhenVec <- rep(NA,length(gens))
    for(j in 1:length(gens)){
      poly1SizeVec[j] <- sum(polyDat[,4] == j)
      poly1PhenVec[j] <- mean(polyDat[which(polyDat[,4] == j),6])
    }
    poly1Phen_coral <- rbind(poly1Phen_coral,poly1PhenVec)
    poly1Size_coral <- rbind(poly1Size_coral,poly1SizeVec)
    oligDat <- read.table(oligObjects[i],header=TRUE)
    olig1SizeVec <- rep(NA,length(gens))
    olig1PhenVec <- rep(NA,length(gens))
    for(j in 1:length(gens)){
      olig1SizeVec[j] <- sum(oligDat[,4] == j)
      olig1PhenVec[j] <- mean(oligDat[which(oligDat[,4] == j),6])
    }
    olig1Phen_coral <- rbind(olig1Phen_coral,olig1PhenVec)
    olig1Size_coral <- rbind(olig1Size_coral,olig1SizeVec)
  print(i)
}

polyExt_coral <- rep(NA,length(gens))
oligExt_coral <- rep(NA,length(gens))
for(i in 1:length(gens)){
  polyExt_coral[i] <- sum(poly1Size_coral[,i] == 0)
  oligExt_coral[i] <- sum(olig1Size_coral[,i] == 0)
}

####### Now mammals
polyObjects <- paste("pedObject_noMajorGene_100Gens_mammal",1:reps,sep="")
oligObjects <- paste("pedObject_majorGene_100Gens_mammal",1:reps,sep="")

poly1Phen_mammal <- NULL
poly1Size_mammal <- NULL
olig1Phen_mammal <- NULL
olig1Size_mammal <- NULL
for(i in 1:reps){
    polyDat <- read.table(polyObjects[i],header=TRUE)
    poly1SizeVec <- rep(NA,length(gens))
    poly1PhenVec <- rep(NA,length(gens))
    for(j in 1:length(gens)){
      poly1SizeVec[j] <- sum(polyDat[,4] == j)
      poly1PhenVec[j] <- mean(polyDat[which(polyDat[,4] == j),6])
    }
    poly1Phen_mammal <- rbind(poly1Phen_mammal,poly1PhenVec)
    poly1Size_mammal <- rbind(poly1Size_mammal,poly1SizeVec)
    oligDat <- read.table(oligObjects[i],header=TRUE)
    olig1SizeVec <- rep(NA,length(gens))
    olig1PhenVec <- rep(NA,length(gens))
    for(j in 1:length(gens)){
      olig1SizeVec[j] <- sum(oligDat[,4] == j)
      olig1PhenVec[j] <- mean(oligDat[which(oligDat[,4] == j),6])
    }
    olig1Phen_mammal <- rbind(olig1Phen_mammal,olig1PhenVec)
    olig1Size_mammal <- rbind(olig1Size_mammal,olig1SizeVec)
  print(i)
}

polyExt_mammal <- rep(NA,length(gens))
oligExt_mammal <- rep(NA,length(gens))
for(i in 1:length(gens)){
  polyExt_mammal[i] <- sum(poly1Size_mammal[,i] == 0)
  oligExt_mammal[i] <- sum(olig1Size_mammal[,i] == 0)
}

#############################################
# get bootstrap CIs for proportion extinct
#############################################
bootReps <- 1000    # number of bootstrap samples
coralPolyBoots <- matrix(NA,nrow=bootReps,ncol=100)
coralOligBoots <- matrix(NA,nrow=bootReps,ncol=100)
mammalPolyBoots <- matrix(NA,nrow=bootReps,ncol=100)
mammalOligBoots <- matrix(NA,nrow=bootReps,ncol=100)

for(i in 1:bootReps){
  mammalPolyDat <- poly1Size_mammal[sample(1:nrow(poly1Size_mammal),nrow(poly1Size_mammal),replace=TRUE),]
  coralPolyDat <- poly1Size_coral[sample(1:nrow(poly1Size_coral),nrow(poly1Size_coral),replace=TRUE),]
  mammalOligDat <- olig1Size_mammal[sample(1:nrow(olig1Size_mammal),nrow(olig1Size_mammal),replace=TRUE),]
  coralOligDat <- olig1Size_coral[sample(1:nrow(olig1Size_coral),nrow(olig1Size_coral),replace=TRUE),]
  coralPolyBoots [i,] <- colSums(coralPolyDat == 0)/nrow(coralPolyDat )
  mammalPolyBoots [i,] <- colSums(mammalPolyDat == 0)/nrow(mammalPolyDat)
  coralOligBoots [i,] <- colSums(coralOligDat == 0)/nrow(coralOligDat)
  mammalOligBoots [i,] <- colSums(mammalOligDat == 0)/nrow(mammalOligDat)
  print(i)
}

mammalPolyCI <- NULL
mammalOligCI <- NULL
coralPolyCI <- NULL
coralOligCI <- NULL

for(i in 1:100){
  mammalPolyCI <- cbind(mammalPolyCI,quantile(mammalPolyBoots[,i],c(0.025,0.975)))
  coralPolyCI <- cbind(coralPolyCI,quantile(coralPolyBoots[,i],c(0.025,0.975)))
  mammalOligCI <- cbind(mammalOligCI,quantile(mammalOligBoots[,i],c(0.025,0.975)))
  coralOligCI <- cbind(coralOligCI,quantile(coralOligBoots[,i],c(0.025,0.975)))
}


############################################
# plot the results
############################################
# coral on the left side
par(mfrow=c(3,2),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)

# coral population sizes
plot(c(0,max(gens)),c(0,6000 + 1500),type="n",ylab="Population size",xlab="",cex.lab=1.5,main = "Coral")
for(i in 1:nrow(poly1Size_coral)){
  lines(1:max(gens),poly1Size_coral[i,],col=alpha("orange",alpha=0.2))
  lines(1:max(gens),olig1Size_coral[i,],col=alpha("darkblue",alpha=0.2))
}
# get rolling mean population sizes
polySizeMean <- rep(NA,max(gens))
oligSizeMean <- rep(NA,max(gens))
for(i in 1:max(gens)){
  polySizeMean[i] <- mean(poly1Size_coral[,i])
  oligSizeMean[i] <- mean(olig1Size_coral[,i])
}
lines(gens,polySizeMean,lwd=3,col="orange")
lines(gens,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(6000,6000),lty="dashed")

par(xpd=TRUE)
text(x=-30,y=8750,labels="A",cex=2)
par(xpd=FALSE)

# mammal population sizes
plot(c(0,max(gens)),c(0,300 + 150),type="n",ylab="",xlab="",cex.lab=1.5,main ="Mammal")
for(i in 1:nrow(poly1Size_mammal)){
  lines(1:max(gens),poly1Size_mammal[i,],col=alpha("orange",alpha=0.2))
  lines(1:max(gens),olig1Size_mammal[i,],col=alpha("darkblue",alpha=0.2))
}
# get rolling mean
polySizeMean <- rep(NA,max(gens))
oligSizeMean <- rep(NA,max(gens))
for(i in 1:max(gens)){
  polySizeMean[i] <- mean(poly1Size_mammal[,i])
  oligSizeMean[i] <- mean(olig1Size_mammal[,i])
}
lines(gens,polySizeMean,lwd=3,col="orange")
lines(gens,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(300,300),lty="dashed")

par(xpd=TRUE)
text(x=-30,y=525,labels="B",cex=2)
par(xpd=FALSE)
#---------------------------------
##### plot the phenotypes
#---------------------------------

par(mar=c(3,5,2,0.5))
# phenotypes for coral
idealPhen <- 110

plot(c(0,max(gens)),c(100,115),type="n",ylab="Mean phenotype",xlab="",cex.lab=1.5)
for(i in 1:length(polyObjects)){
  lines(1:max(gens),poly1Phen_coral[i,],col=alpha("orange",alpha=0.2))
  lines(1:max(gens),olig1Phen_coral[i,],col=alpha("darkblue",alpha=0.2))
}

# get rolling mean
polyPhenMean <- rep(NA,max(gens))
oligPhenMean <- rep(NA,max(gens))
for(i in 1:max(gens)){
  polyPhenMean[i] <- mean(poly1Phen_coral[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(olig1Phen_coral[,i],na.rm=TRUE)
}
lines(gens,polyPhenMean,lwd=3,col="orange")
lines(gens,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,100),c(110,110),lty="dashed")

#### phenotypes for mammal
plot(c(0,max(gens)),c(100,115),type="n",ylab="",xlab="",cex.lab=1.5)
for(i in 1:length(polyObjects)){
  lines(1:max(gens),poly1Phen_mammal[i,],col=alpha("orange",alpha=0.2))
  lines(1:max(gens),olig1Phen_mammal[i,],col=alpha("darkblue",alpha=0.2))
}

# get rolling mean
polyPhenMean <- rep(NA,max(gens))
oligPhenMean <- rep(NA,max(gens))
for(i in 1:max(gens)){
  polyPhenMean[i] <- mean(poly1Phen_mammal[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(olig1Phen_mammal[,i],na.rm=TRUE)
}
lines(gens,polyPhenMean,lwd=3,col="orange")
lines(gens,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,100),c(110,110),lty="dashed")

#####################################################################
# plot the fraction of extinct population versus time
#####################################################################

par(mar=c(4,5,1,0.5))

# extinct proportion for coral
plot(c(1,max(gens)),c(0,0.8),type="n",xlab="Generation",ylab="Proportion extinct",cex.lab=1.5)
lines(gens,polyExt_coral/nrow(poly1Size_coral),col="orange",lwd=2)
polygon(x=c(1:100,rev(1:100)),y=c(coralPolyCI[1,],rev(coralPolyCI[2,])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
lines(gens,oligExt_coral/nrow(olig1Size_coral),col="darkblue",lwd=2)
polygon(x=c(1:100,rev(1:100)),y=c(coralOligCI[1,],rev(coralOligCI[2,])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)
legend(x=40,y=0.2,xjust=FALSE,yjust=FALSE,legend=c("Polygenic","Large-effect locus"),
       lwd=3,col=c("orange","darkblue"),bty="n")

# extinct proportion for mammal
plot(c(1,max(gens)),c(0,0.8),type="n",xlab="Generation",ylab="",cex.lab=1.5)
lines(gens,polyExt_mammal/nrow(poly1Size_mammal),col="orange",lwd=2)
polygon(x=c(1:100,rev(1:100)),y=c(mammalPolyCI[1,],rev(mammalPolyCI[2,])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
lines(gens,oligExt_mammal/nrow(olig1Size_mammal),col="darkblue",lwd=2)
polygon(x=c(1:100,rev(1:100)),y=c(mammalOligCI[1,],rev(mammalOligCI[2,])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)





