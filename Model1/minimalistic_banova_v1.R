



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
npercell = 20 # 4 # 8 # 100
x = c( rep(1, npercell), rep(-1, npercell))
y = c(rnorm( npercell, 0, 1), rnorm( npercell, 1, 1))
y = y - mean(y) # zero center to be sure that b0 == 0

# Specify the data in a form that is compatible with BRugs model, as a list:
dataList <- list(
 y = y , #/ sd(y) ,  
 x = x ,
 Ntotal = 2*npercell # ,
#      tau = 1  
)
#--

#-- RUN THE CHAINS
jagsModel <- jags.model( "model.txt" , data=dataList , inits=list(b=0) , 
                       n.chains=2 , n.adapt=2000, quiet=T )
update(jagsModel , n.iter=2000, n.burnin=1000, progress.bar="none")
codaSamples <- coda.samples( jagsModel , variable.names=c("b") , 
                           n.iter=10000 , thin=1, progress.bar="none" )

# posterior b distribution
# Note, b is symmetric deviation from 0.
#b <- as.numeric( codaSamples[[1]][,"b"] )

# Plot results
plot(codaSamples)

#  Autocorrelation
autocorr(codaSamples)
autocorr.diag(codaSamples)
autocorr.plot(codaSamples)

#  mixin. Gelman-Rubin criterion
#  measures whether there is a significant difference between the variance within several chains 
#  and the variance between several chains by the potential scale reduction factors.
gelman.diag(codaSamples)
gelman.plot(codaSamples)



