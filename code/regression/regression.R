library("mvtnorm")
library("MASS")
library("MCMCpack")
library("matrixcalc")
library(doParallel)

id <- Sys.getenv(x = "SGE_TASK_ID")
ncores <- Sys.getenv(x = "NSLOTS")

p.l<-c(5,10,50,80,90,100) 
N.l<-c(100,200,500,1000) 
rep<-c(1:1)

grid <- expand.grid(rep,N.l,p.l)
grid<-as.data.frame(grid)
colnames(grid)<-c("rep","N","p")

# id specifies the a specific simulation condition
args <- grid[id, ]

do_one <- function(rep,N,p){

# k is from 1 to 1000 (1000 replications)
  estimate<-function(k){

x<-rmvnorm(N, mean = rep(0,p),sigma = diag(p))
beta<-rep(1,p+1)
x2<-cbind(1,x)
y<-x2%*%beta+rnorm(N, 0,sd=0.5)

betahat<-solve(t(x2)%*%x2)%*%t(x2)%*%y
sigmahat<-1/(N-p-1)*t(y-x2%*%beta)%*%(y-x2%*%beta)

t=10^5
sigmae1<-c()
betae1<-matrix(NA,nrow=t,ncol=p+1)

sigmae2<-c()
betae2<-matrix(NA,nrow=t,ncol=p+1)

####### draw posterior samples from posterior marginal distributions; two chains
for (i in 1:t){
  pre.y1<-rgamma(1,shape=(N-p-1)/2,rate=(N-p-1)*sigmahat/2)
  sigmae1[i]<- 1/pre.y1
  betae1[i,]<-  rmvt(1, delta=betahat,sigma =c(sigmahat)*solve(t(x2)%*%x2), df = N-p-1)
  
  pre.y2<-rgamma(1,shape=(N-p-1)/2,rate=(N-p-1)*sigmahat/2)
  sigmae2[i]<- 1/pre.y2
  betae2[i,]<-  rmvt(1, delta=betahat,sigma =c(sigmahat)*solve(t(x2)%*%x2), df = N-p-1)
}
resultsm1<-mcmc(cbind(betae1,sigmae1))
resultsm2<-mcmc(cbind(betae2,sigmae2))

colnames(resultsm1)[p+2]=colnames(resultsm2)[p+2]="sigmae"

############# summarize indices
sum11.psrf=sum10.psrf=sum.u.psrf=max.u.psrf=max.psrf=m.psrf=if.m.psrf=sum.geweke=sum.geweke2=c()
iter<-c(100,500, 1000,3000,5000,10^4,5*10^4,10^5)
k=1
for (j in iter){
  listresult<-mcmc.list(mcmc(resultsm1[1:j,]),mcmc(resultsm2[1:j,]))
  
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
  
  geweke<-try(geweke.diag(mcmc(resultsm1[1:j,]), frac1=0.1, frac2=0.5))
  if (class(geweke)!="try-error" ){
  geweke<-   geweke$z
  sum.geweke[k]<- sum(abs(geweke)>1.96)}else{
    sum.geweke[k]<-NA
  }
  
 k=k+1

}
result<-c(N,p,p+2, sum.u.psrf,max.u.psrf,max.psrf,sum11.psrf,sum10.psrf,
          m.psrf, if.m.psrf,sum.geweke)

#############

fn <- paste0(p,"p",N,"N","result.csv")

# save simulation results
write.table(t(result), file = fn,  sep = ",", append = TRUE,col.names = FALSE)

  }
  
## 1000 replications
  results<-mclapply(1:1000,estimate, mc.cores=ncores)

}

do.call(do_one, args)
