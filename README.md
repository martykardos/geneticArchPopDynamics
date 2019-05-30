# geneticArchPopDynamics
pedR_0.1.5.tar.gz is an R package for forward -time individual-based simulations of selection on a quantitative trait with linkage; this was used for the forward-time simulations that included linkage in our paper. 

quantPop_0.1.5.tar.gz is an R package implementation of the deterministic eco-evolutionary models, and the individual-based simulations of selection on a quantitative trait without linkage used in our paper.

rCode_deterministicPopDynamicsExamples_27May2019.R is an R script with example scenarios of the deterministic models. 

rCode_stochasticPopDynamicsExamples_27May2019.R is an R script with example scenarios of the stochastic, individual-based simulations we used in our paper. 

The function help files define the input parameters. 


Main function outputs:

pedSimHardSelecFast 
NVec: population sizes for each simulated generation
pedObject: data frame with 6 columns including id (identification), mom (mother ID),
dad (father ID), ten (generation), immVec (0 = non-immigrant; 1 = immigrant), indTraitVal (phenotype value)
quantTraitDat: data frame with six columns including generation (col. 1), VaVec (additive genetic variance), VpVec (phenotype variance), h2Vec (heritability, narrow sense), traitMean (mean trait value), sVec (selection differential)
qtlInfo: data frame with one row for each quantitative trait locus. columns are qtly (locus IDs), qtlExpHet (heterozygosity), qtlEffs (parameter a [defined in the paper]),qtlChrs (chromosome where the qtl is located), qtlPhysPos (location of the QTL on the chromosome in Mb), the remaining columns are the allele frequencies each generation simulated. 

logQuantIndiv_varP0_stocPhenOpt
NVec: vector of population sizes for each generation
freqMat: data frame with one row for each locus. The allele frequencies for each generation are in columns.
meanFit: mean fitness in each generation
vpVec: vector of total phenotype variance for each generation
h2Vec: vector of heritability for each generation
phenVec: mean phenotype each generation
qtlInfo: data frame with one row for each locus. columns are allFreq (starting allele frequency), aVec (parameter a for each locus), and vVec (starting additive genetic variance for each locus).

logQuant outputs:
NVec: vector of population sizes each generation
freqVec: allele frequency for the QTLs each generation
meanFit: mean fitness in the population each generation
vpVec: phenotype variance each generation
h2Vec: heritability (narrow sense) each generation
phenVec: mean phenotype each generation
vQTL, Ve, and a are the locus-specific additive genetic variance, environmental variance, and the parameter a (half the phenotype difference between the two homozygous genotypes).

outputs for the logQuantPseudoGenos function are the same as for logQuant.
