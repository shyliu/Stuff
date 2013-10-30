Data=read.csv("blr_data_1001.csv",header=TRUE)
Par=read.csv("blr_pars_1001.csv",header=TRUE)


####################write MCMC function
####################mvrnorm() is in package(MASS)
library(MASS)

#######bayes.logreg() is the function to generate beta from posterior
bayes.logreg <- function (m,y,X,beta.0,Sigma.0.inv,
                   niter=10000,burnin=1000,
                   print.every=10000,retune=100,
                   verbose=TRUE){

######set up initial value for estimated beta and covariance matrix sigma
######draws() will be the time sequence of the iterated estimates

      beta.cur=beta.0
      sigma.cur=Sigma.0.inv
      draws=c()

######creat a function that generate a set of data from multivariate normal distribution
######with designed covariance matrix sigma
      
      jumping.multi.sigma=function(Beta.input,Sigma.input){
              mvrnorm(n=1,mu=Beta.input,Sigma=Sigma.input)}

######post() is posterior likelihood function 

      post=function(beta.input){
      part1=(-m)*log(1+exp(X%*%beta.input))+(X%*%beta.input)*y
      part2=t(beta.input-beta.0)%*%solve(Sigma.0.inv)%*%(beta.input-beta.0)
      return(sum(part1)-0.5*part2)}

######beta.update() is a function that will update the estimated beta
######input: current beta with current covariance matrix which 
######is used for proposal density
######output: updated estimate

      beta.update=function(beta.input,sigma.input){
      beta.can=matrix(jumping.multi.sigma(beta.input,sigma.input),nrow=2)
      compare=post(beta.can)-post(beta.input)
      logU=log(runif(1))
      if (logU<=compare){
       result=c(beta.can,1)
      }
      else
       result=c(beta.input,0)
      return(result)    
      }

######start the draws from mean of covariates 
######and use beta.update() to update the estimates
######check acceptance rate and tune the covariance matrix for proposal density
######decrease the covariance matrix if acceptance is greater than 0.3
######increase the covariance matrix if acceptance is greater than 0.5

      draws=matrix(c(colMeans(X),1),nrow=1)
      for (x in 2:burnin){
      beta.cur=matrix(draws[x-1,1:2],nrow=2)
      if (x %% retune == 0){
       accept.rate=mean(draws[(x-retune+1):(x-1),3])
       if (accept.rate>=0.5){
       where=which(draws[(x-99):(x-1),3]==1)
       Cov=cov(draws[x-retune+where,1:2])
       sigma.cur=accept.rate*Cov+(1-accept.rate)*sigma.cur}
       if (accept.rate<0.3){sigma.cur=accept.rate*sigma.cur}
       print(list(paste("accept-ratio=",round(accept.rate,3)),Sigma=sigma.cur))
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
       Beta.ite=matrix(Moredraws[t-1,1:2],nrow=2)
       newd=matrix(beta.update(Beta.ite,sigma.cur),nrow=1)
       Moredraws=rbind(Moredraws,newd) 
       if(t  %% print.every == 0){
       accept.rate=mean(Moredraws[(t-print.every+1):(t-1),3])
       print(paste("update-accept-ratio=",round(accept.rate,3)))}
      }

######After we have estimates, calculate quantiles and confidence interval 
######for the estimates
      Beta.1.quant=matrix(quantile(Moredraws[,1],seq(0.01,0.99,0.01)),ncol=1)
      Beta.2.quant=matrix(quantile(Moredraws[,2],seq(0.01,0.99,0.01)),ncol=1)
      Beta.quantile=cbind(Beta.1.quant,Beta.2.quant)

      CI.beta1=matrix(quantile(Moredraws[,1],c(0.025,0.975)),nrow=1)
      CI.beta2=matrix(quantile(Moredraws[,2],c(0.025,0.975)),nrow=1)
      CI.beta=rbind(CI.beta1,CI.beta2)
      return(list(Beta.quantile,CI.beta))

}     #end of bayes.logreg()

 
#####################################################################
y=matrix(Data[,1])
m=matrix(Data[,2])
X=as.matrix(Data[,3:4])
p=2
Sigma.0.inv <- diag(rep(1.0,p))
beta.0=matrix(0,nrow=p)

bayes.logreg <- function (m,y,X,beta.0,Sigma.0.inv,
                   niter=10000,burnin=10000,
                   print.every=10000,retune=100,
                   verbose=TRUE)


######check MCMC algorithm###################
NewData=Data
identity1=rep(0,200)
identity2=rep(0,200)

######generate 200 new data set by the estimate generated from MCMC  
   NewDataList=lapply(1:200,function(rep){
    beta.gen=mvrnorm(n=1,mu=matrix(0,nrow=2),Sigma=diag(2))
    for (yi in 1:nrow(NewData)){
    logit=function(u){exp(u)/(1+exp(u))}
    NewData[yi,1]=sum(rbinom(NewData[yi,2],size=1,
                  prob=c(logit(X%*%beta.gen)[yi,1])))
    }
    return(list(NewData,beta.gen))
   })
   
######calculate credible interval for 200 dataset
for(G in 1:200){
  Data.new=NewDataList[[G]][[1]]
  y.new=matrix(Data.new[,1])
  m.new=matrix(Data.new[,2])
  X.new=as.matrix(Data.new[,3:4])

  Cred.inv=bayes.logreg(m.new,y.new,X.new,beta.0,Sigma.0.inv,
                   niter=10000,burnin=1000,
                   print.every=1000,retune=100,
                   verbose=TRUE) 
  Find1=findInterval(NewDataList[[G]][[2]][1],Cred.inv[[2]][1,])
  Find2=findInterval(NewDataList[[G]][[2]][2],Cred.inv[[2]][2,])
  if (Find1==1){identity1[G]=G}
  if (Find2==1){identity2[G]=G}
}


