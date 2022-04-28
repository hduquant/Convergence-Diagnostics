library("mvtnorm")
library("MASS")
library("MCMCpack")
library("matrixcalc")
library(doParallel)

id <- Sys.getenv(x = "SGE_TASK_ID")
ncores <- Sys.getenv(x = "NSLOTS")

p.l<-c(5,10,12,20,50,100)
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
    
    sigma0<-matrix(0.3,nrow=p,ncol=p)
    diag(sigma0)<-1
    mean0<-rep(0,p)
    y<- rmvnorm(N, mean = mean0,sigma = sigma0)
    
    xbar<-colMeans(y)
    
    fun1<-function(i){
      (y[i,]-xbar)%*%t(y[i,]-xbar)
    }
    ee<-mclapply(1:nrow(y),fun1, mc.cores=4)
    S=Reduce('+', ee)
    
    t=10^5
    sigma1<-list()
    mu1<-matrix(NA,nrow=t,ncol=p)
    sigma2<-list()
    mu2<-matrix(NA,nrow=t,ncol=p)
    
####### draw posterior samples from posterior marginal distributions; two chains
    for (i in 1:t){
      sigma1[[i]]<-riwish(N-1,S)
      mu1[i,]<-  rmvt(1, delta=xbar,sigma =S/((N-p)*N), df = N-p)
      
      sigma2[[i]]<-riwish(N-1,S)
      mu2[i,]<-  rmvt(1, delta=xbar,sigma =S/((N-p)*N), df = N-p)
    }
    sigma1.2<-lapply(sigma1,vech)
    sigma1.3<-matrix(unlist(sigma1.2),nrow=t,byrow = T)
    resultsm1<-mcmc(cbind(mu1,sigma1.3))
    
    sigma2.2<-lapply(sigma2,vech)
    sigma2.3<-matrix(unlist(sigma2.2),nrow=t,byrow = T)
    resultsm2<-mcmc(cbind(mu2,sigma2.3))
    
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
    result<-c(N,p,(p*(p+1)/2+p), sum.u.psrf,max.u.psrf,max.psrf,sum11.psrf,sum10.psrf,
              m.psrf, if.m.psrf,sum.geweke)
    
    fn <- paste0(p,"p",N,"N","result.csv")

# save simulation results    
    write.table(t(result), file = fn,  sep = ",", append = TRUE,col.names = FALSE)
    
  }
  
## 1000 replications
  results<-mclapply(1:1000,estimate, mc.cores=ncores)
  
}

do.call(do_one, args)



#############



