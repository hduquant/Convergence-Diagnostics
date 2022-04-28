library("mvtnorm")
library("MASS")
library("MCMCpack")
library("matrixcalc")
library(doParallel)
require(rjags)

id <- Sys.getenv(x = "SGE_TASK_ID")
ncores <- Sys.getenv(x = "NSLOTS")

N.l<-c(50,100,200,500,1000)
rep<-c(1:1)

grid <- expand.grid(rep,N.l)
grid<-as.data.frame(grid)
colnames(grid)<-c("rep","N")

# id specifies the a specific simulation condition
args <- grid[id, ]

do_one <- function(rep,N){

# k is from 1 to 1000 (1000 replications)
  estimate<-function(k){
#p=5, one factor model
    p=5
    m=1
    ms=m*(m-1)/2
    q=2*p+ms
    ps=p*(p+1)/2
    p2=p*2
    df=ps-q
    p_m=p/m
set.seed(1)
lamb<-sort(sample(seq(70,95,5),p,replace=TRUE))/100

set.seed(k)
lamb_t <- t(lamb)
phi <- 0.5*diag(m)+matrix(0.5,m,m)
lamblamb <- lamb%*%phi%*%lamb_t

psi <- diag(p)-diag(diag(lamblamb))

sig <- lamblamb+psi    #BIG SIGMA
eigen_sig <- eigen(sig)
sig_vec <- eigen_sig$vectors
sig_val <- eigen_sig$values
sig12 <-  sig_vec%*%diag(sqrt(sig_val))%*%t(sig_vec) # sqrt of BIG SIGMA

z <- rmvnorm(N, mean=rep(0,p))
x <- z%*%sig12

### Use jags to draw posterior samples
model= "
model {

for (i in 1:N){
     f[i]~dnorm(0, pre.phi)
for (j in 1:p){
mux[i,j]<-l[j]*f[i]
x[i,j]~dnorm(mux[i,j],pre.p[j])
}
   }
  
  for (i in 1:p){
  pre.p[i]~dgamma(.001,.001)
  pvar[i]<-1/pre.p[i]
  }
pre.phi~dgamma(.001,.001)
  phi<-1/pre.phi



for (i in 1:p){
 l[i]~dunif(0, 1.0E6)}

  }
  "
  # Save model
  writeLines(model, con="model.txt" )
  #------------------------------------------------------------------------------
  # Load data
  dataList = list(
    x=x, # Read in the data of first group
    N = N,
    p=p
  )
  #------------------------------------------------------------------------------
  # Specifying starting values in two independent chains
 load <- rep(1,p)
 pre.phi <- c(1)
 pre.p  <- rep(1,p)
  initsList1 = list(  l = load , pre.phi = pre.phi,  pre.p= pre.p,
                     .RNG.name="base::Wichmann-Hill",.RNG.seed=2018)
  load <- rep(0.05,p)
  pre.phi <- c(50)
  pre.p  <- rep(0.5,p)
  initsList2 = list( l= load , pre.phi = pre.phi,  pre.p= pre.p,
                     .RNG.name="base::Wichmann-Hill",.RNG.seed=2016)
  #------------------------------------------------------------------------------
  parameters = c( "l ","phi", "pvar")     # Specify the estimated parameters 
  adaptSteps =0              # Adaptive period
  burnInSteps = 1000         # Burn-in period
  nChains = 1
  thinSteps=1                # Thinning period
  numSavedSteps=10^5      # The number of kept iterations
  nIter = ceiling( numSavedSteps * thinSteps )
  jagsModel1 = jags.model( "model.txt" , data=dataList , inits=initsList1, 
                           n.chains=nChains , n.adapt=adaptSteps )
  jagsModel2 = jags.model( "model.txt" , data=dataList , inits=initsList2, 
                           n.chains=nChains , n.adapt=adaptSteps )
  update( jagsModel1 , n.iter=burnInSteps)
  update( jagsModel2 , n.iter=burnInSteps)
  codaSamples1 = coda.samples( jagsModel1 , variable.names=parameters, 
                               n.iter=nIter , thin=thinSteps)
  codaSamples2 = coda.samples( jagsModel2 , variable.names=parameters, 
                               n.iter=nIter , thin=thinSteps)
  mcmcChain1 = as.matrix( codaSamples1 )
  mcmcChain2 = as.matrix( codaSamples2 )
 # head(mcmcChain1)
 # head(mcmcChain2)
 
  resultsm1<-mcmc(  mcmcChain1)
  resultsm2<-mcmc(  mcmcChain2)
  
############# summarize indices
  sum11.psrf=sum10.psrf=sum.u.psrf=max.u.psrf=max.psrf=m.psrf=if.m.psrf=c()
  sum.geweke=c()
iter<-c(100,500, 1000,3000,5000,10^4,5*10^4,10^5)
k=1
for (j in iter){
  
  listresult<-mcmc.list(mcmc(resultsm1[(j/2):j,]),mcmc(resultsm2[(j/2):j,]))
  
  gelman<-try(gelman.diag(listresult))
  if (class(gelman)!="try-error" ){
    sum.u.psrf[k]<-sum(  gelman$psrf[,2]>=1.1)
    max.u.psrf[k]<-max(  gelman$psrf[,2])
    max.psrf[k]<-max(  gelman$psrf[,1])
    sum11.psrf[k]<-sum(  gelman$psrf[,1]>=1.1)
    sum10.psrf[k]<-sum(  gelman$psrf[,1]>=1)
    m.psrf[k]<-  gelman$ mpsrf
    if.m.psrf[k]<- m.psrf[k]>1.1} else {
      sum.u.psrf[k]<-NA
      max.u.psrf[k]<-NA
      max.psrf[k]<-NA
      sum11.psrf[k]<-NA
      sum10.psrf[k]<-NA
      m.psrf[k]<-NA
      if.m.psrf[k]<-NA
    }
  
  geweke<-try(geweke.diag(mcmc(resultsm1[(j/2):j,]), frac1=0.1, frac2=0.5))
  if (class(geweke)!="try-error" ){
    geweke<-   geweke$z
    sum.geweke[k]<- sum(abs(geweke)>1.96)}else{
      sum.geweke[k]<-NA
    }
  
 k=k+1

}
result<-c(N,p,m, q, df,sum.u.psrf1,max.u.psrf1,max.psrf1,m.psrf1, if.m.psrf1, sum11.psrf1, sum10.psrf1,
          sum.u.psrf2,max.u.psrf2,max.psrf2,m.psrf2, if.m.psrf2, sum11.psrf2, sum10.psrf2,
          sum.geweke,sum.geweke2)

#############

fn <- paste0(N,"N","result.csv")

# save simulation results
write.table(t(result), file = fn,  sep = ",", append = TRUE,col.names = FALSE)

  }

## 1000 replications
  results<-mclapply(1:1000,estimate, mc.cores=ncores)
 
}

do.call(do_one, args)
