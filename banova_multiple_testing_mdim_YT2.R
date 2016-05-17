# graphics.off()
# rm(list=ls(all=TRUE))

require(rjags) 
                       
set.seed(7)

#-- THE DATA.
#-------------------------------------------------------
# let's generate iTRAQ-like dataset
nProt = 10
nPerExp = 4
itraq = t(replicate( nProt, rnorm( nPerExp*2 ) ))
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
            mu[i, j] <- b*x[j];
            y[i,j] ~ dnorm(mu[i,j], tau)
        }
    }

    # Priors
    tau <- pow( sigma , -2 );
    sigma ~ dunif(0,10);
    b ~ dnorm( 0.0 , 1e-12 )
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
#-- RUN THE CHAINS
jagsModel = jags.model( "model.txt" , data=dataList, inits=list(b=0) , n.chains=1 , n.adapt=1000, quiet=T )
update( jagsModel , n.iter=1000, progress.bar="none" )
codaSamples = coda.samples( jagsModel , variable.names=c("b") , 
                            n.iter=10000 , thin=1, progress.bar="none" )
#------------------------------------------------------------------------------
#-- DERIVE THE P- AND P-LIKE VALUES
b = as.numeric( codaSamples[[1]][,"b"] )
ttestProb = t.test( y[x==1], y[x==-1], var.equal=T )$p.value
mcmcProb = min( c(2*sum(b>0)/length(b), 2*sum(b<0)/length(b))) 
#------------------------------------------------------------------------------












