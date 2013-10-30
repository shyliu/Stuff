Data=read.table("breast_cancer.txt",header=TRUE)
Y=matrix(0,nrow=nrow(Data))
Y[which(Data$diagnosis=="M"),1]=1
Data$diagnosis=Y[,1]

Data$intercept=rep(1,nrow(Data))
y=matrix(Data$diagnosis)
m=matrix(rep(1,nrow(Data)))
X=as.matrix(Data[,-(ncol(Data)-1)])
for (i in 1:(ncol(X)-1)){
   X[,i]=(X[,i]-mean(X[,i]))/sd(X[,i])
}
beta.0=matrix(0,nrow=ncol(X))
Sigma.0.inv=1000*diag(ncol(X))

library(MASS)
library(coda)

######fit glm and use the covariace matrix in glm as 
######the covariance matrix in proposal density
GLM=glm(diagnosis~area+compactness+concavepts+concavity+
        fracdim+perimeter+radius+smoothness+symmetry+texture,
        data=Data,family=binomial(link = "logit"))
GLM_cov=summary(GLM)$cov.unscaled


#######bayes.logregMG() is the function to generate beta from posterior
bayes.logregMG <- function (m,y,X,beta.0,Sigma.0.inv,
                   niter,burnin,
                   print.every,retune,
                   verbose=TRUE){


######set up initial value for estimated beta and covariance matrix sigma
      beta.cur=beta.0
      sigma.cur=GLM_cov

######creat a function that generate a set of data from multivariate normal distribution
######with designed covariance matrix sigma
      
      jumping.multi.sigma=function(Beta.cur,Sigma.cur){
              mvrnorm(n=1,mu=Beta.cur,Sigma=Sigma.cur)}

######post() is posterior likelihood function 

      post=function(beta.input){
      part1=(-m)*log(1+exp(X%*%beta.input))+(X%*%beta.input)*y
      part2=t(beta.input-beta.0)%*%solve(Sigma.0.inv)%*%(beta.input-beta.0)
      return(sum(part1)-0.5*part2)}

######beta.update() is a function that will update the estimated beta
######one by one(conditional)
######input: current beta with current covariance matrix which 
######is used for proposal density
######output: updated estimate

      beta.update=function(beta.input,sigma.input){
      beta.c=beta.input
      yes=0
      for (b in 1:ncol(X)){
         beta.get=matrix(jumping.multi.sigma(beta.c,sigma.input),ncol=1)
         beta.can=beta.c
         beta.can[b,1]=rnorm(n=1,mean=beta.get[b,1],sd=sqrt(sigma.cur[b,b]))
         compare=post(beta.can)-post(beta.input)
         logU=log(runif(1))
         if (logU<=compare){
         yes=yes+1
         beta.c=beta.can}
      }
      if (yes >= 7){dec=1}
      else {dec=0}
      return(c(beta.c,dec))
      }

######start the draws from zero
######and use beta.update() to update the estimates
######check acceptance rate and tune the covariance matrix for proposal density
######decrease the covariance matrix if acceptance is greater than 0.3
######increase the covariance matrix if acceptance is greater than 0.5

      draws=matrix(c(rep(0,ncol(X)),0),nrow=1)
      sigma.cur=diag(ncol(X))
      for (x in 2:burnin){
      beta.cur=matrix(draws[x-1,-(ncol(X)+1)],ncol=1)
      if (x %% retune == 0){
       accept.rate=mean(draws[(x-retune+1):(x-1),(ncol(X)+1)])
       if (accept.rate>=0.5){
       where=which(draws[(x-retune+1):(x-1),(ncol(X)+1)]==1)
       Cov=cov(draws[x-retune+where,1:ncol(X)])
       sigma.cur=accept.rate*Cov+(1-accept.rate)*sigma.cur}
       if (accept.rate<0.3){sigma.cur=(accept.rate+0.01)*sigma.cur}
       print(list(paste("accept-ratio=",round(accept.rate,3)),
             Sigma=sigma.cur))
       }
      renew=matrix(beta.update(beta.cur,sigma.cur),nrow=1)
      draws=rbind(draws,renew)
      }
######We can abandon estimates in draws()

######After burnin process in draws() we fix covariance matrix for proposal density
######and we continue to sample from posterior distribution
######Moredraws are the estimate we use.

      Moredraws=matrix(c(draws[burnin,]),nrow=1)
      for (t in 2:niter){
       Beta.ite=matrix(Moredraws[t-1,-(ncol(X)+1)],ncol=1)
       newd=matrix(beta.update(Beta.ite,sigma.cur),nrow=1)
       Moredraws=rbind(Moredraws,newd) 

#######print the acceptance rate at print.every
       if(t  %% print.every == 0){
       accept.rate=mean(Moredraws[(t-print.every+1):(t-1),(ncol(X)+1)])
       print(paste("update-accept-ratio=",round(accept.rate,3)))}
      }

######After we have estimates, calculate quantiles and confidence interval 
######for the estimates

      Beta.quantile=c()
      for(j in 1:ncol(X)){
      quant=matrix(quantile(Moredraws[,j],seq(0.01,0.99,0.01)),ncol=1)
      Beta.quantile=cbind(Beta.quantile,quant)
      }

      CI.beta=c()
      for (i in 1:ncol(X)){
      CI=matrix(quantile(Moredraws[,i],c(0.025,0.975)),nrow=1)
      CI.beta=rbind(CI.beta,CI)
      }
      return(list(Beta.quantile,CI.beta,Moredraws))
}   #end of bayes.logregMG()

#############################################################################
######sample posterior distribution from MCMC algorithm
CheckG=bayes.logregMG(m,y,X,beta.0,Sigma.0.inv,
                   niter=200000,burnin=20000,
                   print.every=1000,retune=100,
                   verbose=TRUE)
save(CheckG,file="CheckG.RData")

##############################################################################
C=CheckG[[3]]
New.par=C[C[,12]==1,]
NewData=Data

######generate 200 new dataset by using previous estimate
NewDataList=lapply(1:200,function(B){
  BETA=matrix(New.par[B,],ncol=1) 
  for (yi in 1:nrow(NewData)){
    logit=function(u){exp(u)/(1+exp(u))}
    p=c(logit(X%*%BETA)[yi,1])
    NewData[yi,11]=sample(c(0,1),size=1,prob=c(1-p,p))
    }
    return(list(NewData,BETA))
})

######calculate credible interval for 200 dataset
identity1=rep(0,200)
identity2=rep(0,200)

for(G in 1:200){
  Data.new=NewDataList[[G]][[1]]
  y.new=matrix(Data.new[,11])
  X.new=as.matrix(Data.new[,-11])
  m.new=Data.new[,12]
  
  for (i in 1:(ncol(X.new)-1)){
   X.new[,i]=(X.new[,i]-mean(X.new[,i]))/sd(X.new[,i])
  }
  Cred.inv=bayes.logregMG2(m.new,y.new,X.new,beta.0,Sigma.0.inv,
                   niter=10000,burnin=1000,
                   print.every=1000,retune=100,
                   verbose=TRUE) 
  Find1=findInterval(NewDataList[[G]][[2]][1],Cred.inv[[2]][1,])
  Find2=findInterval(NewDataList[[G]][[2]][2],Cred.inv[[2]][2,])
  if (Find1==1){identity1[G]=G}
  if (Find2==1){identity2[G]=G}
}

identityG1=identity1
identityG2=identity2
NewDataListG=NewDataList

save(identityG1,file="identityG1.RData")
save(identityG2,file="identityG2.RData")
save(NewDataListG,file="NewDataListG.RData")

