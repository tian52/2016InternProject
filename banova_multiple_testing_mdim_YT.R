# graphics.off()
# rm(list=ls(all=TRUE))

require(rjags) 
                       
set.seed(7)



#-- minimalistic t-test model---------------------------
modelstring = "
model {

   for ( i in 1:Ntotal ) {
      y[i] ~ dnorm( mu[i] , tau )
      mu[i] <- b*x[i]
  }

   tau <- pow( sigma , -2 )
   sigma ~ dunif(0,10) # the standard deviation

   b ~ dnorm( 0.0 , 1e-12 )  # flat prior on b
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------


#-- THE DATA.
#-------------------------------------------------------
# let's generate iTRAQ-like dataset
nProt = 10
nPerExp = 4
itraq = t(replicate( nProt, rnorm( nPerExp*2 ) ))
itraq = itraq - rowMeans(itraq)
p.vals = apply( itraq, 1, function(x) {t.test(x[1:nPerExp], x[-(1:nPerExp)])$p.value})
groupDefs = c( rep(0,nPerExp), rep(+1, nPerExp) )
# how do I feed this into JAGS?

blist <- list()
ttestProb <- vector(length = nProt)
mcmcProb <- vector(length = nProt)
bm <- c(1:10000)

# Consider each row independently:

for (j in 1:nProt) {
  y = itraq[j,]
  x = groupDefs
# Specify the data in a form that is compatible with BRugs model, as a list:
dataList <- list(
  y= itraq[j,],
  x= groupDefs,
  Ntotal = nPerExp*2
)

#------------------------------------------------------------------------------
#-- RUN THE CHAINS
jagsModel = jags.model( "model.txt" , data=dataList, inits=list(b=0) , 
                        n.chains=1 , n.adapt=1000, quiet=T )
update( jagsModel , n.iter=1000, progress.bar="none" )
codaSamples = coda.samples( jagsModel , variable.names=c("b") , 
                            n.iter=10000 , thin=1, progress.bar="none" )
#------------------------------------------------------------------------------
#-- DERIVE THE P- AND P-LIKE VALUES
b = as.numeric( codaSamples[[1]][,"b"] )
plot(b,main=j)
bm = rbind(bm,b)
ttestProb[j] = t.test( y[x==1], y[x==-1], var.equal=T )$p.value
mcmcProb[j] = min( c(2*sum(b>0)/length(b), 2*sum(b<0)/length(b))) 
rm(dataList)
#------------------------------------------------------------------------------
}












