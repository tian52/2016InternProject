# graphics.off()
# rm(list=ls(all=TRUE))

require(rjags) 
                       
set.seed(7)

#-- THE DATA.
#-------------------------------------------------------
# let's generate iTRAQ-like dataset
nProt = 1000
nPerExp = 4
itraq1 = t(replicate( nProt, rnorm( nPerExp ) ))
itraq2 = t(replicate( nProt, rnorm( nPerExp, 1, 1 ) ))
itraq<-cbind(itraq1,itraq2)
itraq = itraq - rowMeans(itraq)
p.vals = apply( itraq, 1, function(x) {t.test(x[1:nPerExp], x[-(1:nPerExp)])$p.value})
groupDefs = c( rep(-1,nPerExp), rep(+1, nPerExp) )
#groupDefsm = t(replicate(nProt, groupDefs))

#-- minimalistic t-test model---------------------------

modelstring <- "
model {
# Model
for (i in 1:Ntotal) {
for (j in 1:Mtotal) {
mu[i, j] <- b[i]*x[j];
y[i,j] ~ dnorm(mu[i,j], tau[i])
}
}
for(i in 1:Ntotal){
sigma[i] ~ dunif(L[i], R[i]);
tau[i] <- pow( sigma[i] , -2 );
b[i] ~ dnorm( 0.0 , T[i] );
}
# Priors
for(i in 1:Ntotal){
T[i] ~ dgamma(0.001, 0.001); # mode=0.1,sd=10.0
L[i] ~ dunif(0,10); # left
width[i] ~ dunif(0,10); # width of uniform dist
R[i] <- L[i] + width[i]
}
}
"
writeLines(modelstring,con="model.txt")
#------------------------------------------------------

y= itraq
x= groupDefs

# how do I feed this into JAGS?

#blist <- list()
#ttestProb <- vector(length = nProt)
#mcmcProb <- vector(length = nProt)
#bm <- c(1:10000)

# Consider each row independently:

#  y = itraq
#  x = groupDefsm
# Specify the data in a form that is compatible with BRugs model, as a list:
dataList <- list(
  y= itraq,
  x= groupDefs,
  Ntotal = nProt,
  Mtotal = nPerExp*2
)

#------------------------------------------------------------------------------
#--  RUN THE CHAINS (with burn-in)
jagsModel = jags.model( "model.txt" , data=dataList, 
                        inits=list(b=rep(0,nProt)) , n.chains=2 , 
                        n.adapt=20000, quiet=T )
update( jagsModel , n.iter=20000, n.burnin=10000, progress.bar="none" )
codaSamples = coda.samples( jagsModel , variable.names=c("b") , 
                            n.iter=10000 , thin=1, progress.bar="none" )
#------------------------------------------------------------------------------

head(as.matrix(codaSamples))

# HDI

source("HDIofMCMC.R")

credMass <- 0.95
HDI <- data.frame(-1, 1)
samplev <- data.frame(c(rep(0,20000)), 1:nProt)
for (i in 1:nProt) {
  samplev[,i] <- unlist(codaSamples[,i])
  HDI[i,] = HDIofMCMC( samplev[,i] , credMass )
}

# Report the cases that the 95% HDI does not contain 0
count <-0 
prUnequal <- {}
for (i in 1:nProt) {
  if (HDI[i,1]>0||HDI[i,2]<0) {
    count = count+1;
    prUnequal <- {c(prUnequal,i)};
  }
}
unequalPerc <- count/nProt*100  

cat("There are ",unequalPerc," % proteins detected to have unequal expression levels in two groups.")
cat("The proteins with unequal expression levels in two groups are: ", prUnequal)


#frequentists approach(t-test)
library(genefilter)
itraqm<-as.matrix(itraq)
groupDefsv <- as.factor(groupDefs)
tTest<-rowttests(itraqm, groupDefsv, tstatOnly = FALSE)
tTestunequal <-tTest[tTest[,3]<0.05,]
TprUnequal <- as.numeric(row.names(tTestunequal))
commenPr <- intersect(TprUnequal,prUnequal)
tTestunequal

# Plot results
#plot(codaSamples)

#  Autocorrelation
#autocorr(codaSamples)
#autocorr.diag(codaSamples)
#autocorr.plot(codaSamples)

#  mixin. Gelman-Rubin criterion
#  measures whether there is a significant difference between the variance within several chains 
#  and the variance between several chains by the potential scale reduction factors.
#gelman.diag(codaSamples)
#gelman.plot(codaSamples)

# #-- DERIVE THE P- AND P-LIKE VALUES
# b = as.numeric( codaSamples[[1]][,"b"] )
# ttestProb = t.test( y[x==1], y[x==-1], var.equal=T )$p.value
# mcmcProb = min( c(2*sum(b>0)/length(b), 2*sum(b<0)/length(b))) 
# #------------------------------------------------------------------------------














