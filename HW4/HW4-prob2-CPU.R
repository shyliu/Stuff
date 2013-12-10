######mvrnorm() is in package(MASS)
library(MASS)

######start probit_mcmc_cpu()

probit_mcmc_cpu <- function (y,X,beta.0,Sigma.0.inv,
                   niter,burnin,
                   print.every,retune,
                   verbose=TRUE){

######set up initial value for estimated beta and covariance matrix sigma
######draws() will be the time sequence of the iterated estimates

      beta.cur=beta.0
      sigma.cur=Sigma.0.inv
      draws=list()


######creat a function that generate a "beta" 
######from multivariate normal distribution
######with designed covariance matrix sigma
      
      beta.generate=function(Z.input,Sigma.input){
           B.beta.gen=solve(solve(Sigma.input)+t(X)%*%X)
           mu.beta.gen=B.beta.gen%*%(solve(Sigma.input)%*%beta.0+t(X)%*%Z.input)
           Beta=matrix(mvrnorm(n=1,mu=mu.beta.gen,Sigma=B.beta.gen))
           return(Beta)} #end of beta.generate


######creat a function that generate a "Z" 
######from truncated multivariate normal distribution

      z.generate=function(beta.input){
      ZZ=matrix(0,nrow=nrow(y),ncol=1)
       for (i in 1:nrow(ZZ)){
         mu.z.gen=t(X[i,])%*%beta.input
         if (y[i,1]==0){ZZ[i,1]=Trunc_Normal(1,mu.z.gen,1,-Inf,0,10000)}
         else {ZZ[i,1]=Trunc_Normal(1,mu.z.gen,1,0,Inf,10000)}
       } #end of for loop     
      return(ZZ)                 
      } #end of z.generate

######beta.post is posterior likelihood function

      beta.post=function(beta.input,Z.input){
      B.beta.post=solve(solve(Sigma.0.inv)+t(X)%*%X)
      mu.beta=B.beta.post%*%(solve(Sigma.0.inv)%*%beta.0+t(X)%*%Z.input)
      part1=-nrow(y)*log(2*pi)/2-nrow(y)*det(B.beta.post)/2
      part2=t(beta.input-mu.beta)%*%solve(B.beta.post)%*%(beta.input-mu.beta)
      return(part1-0.5*part2)}


######par.update() is a function that will update the estimated beta & z
######input: current beta,z and covariate matrix
######output: updated estimate

    par.update=function(beta.input,z.input,sigma.input){
        Z.par.update=z.generate(beta.input)
        beta.par.cand=beta.generate(Z.par.update,sigma.input)
        compare=beta.post(beta.par.cand,Z.par.update)-beta.post(beta.input,
                Z.par.update)
        logU=log(runif(1))
        if (logU<=compare){
         result=list(beta.par.cand,Z.par.update,1)
        }
        else{
         result=list(beta.input,Z.par.update,0)}
        return(result)
    } #end of par.update

######start the draws from mean of covariates 
######and use beta.update() to update the estimates
######check acceptance rate and tune the covariance matrix for proposal density
######decrease the covariance matrix if acceptance is greater than 0.3
######increase the covariance matrix if acceptance is greater than 0.5

      
      GLM=glm(y~X_1+X_2+X_3+X_4+X_5+X_6+X_7+X_8-1,
          data=Data,family=binomial(link = "probit"))
      beta.cur=summary(GLM)$coefficients[,1]
      z.cur=z.generate(beta.cur)
      draws[[1]]=list(beta.cur,z.cur,0)
      sigma.cur=summary(GLM)$cov.unscaled
      for (x in 2:burnin){
      beta.cur=draws[[x-1]][[1]]
      z.cur=draws[[x-1]][[2]]
      if (x %% retune == 0){
       num.a=lapply((x-retune+1):(x-1),function(tt){draws[[tt]][[3]]})
       accept.rate=mean(unlist(num.a))
       if (accept.rate>=0.5){
         sigma.cur=(1-accept.rate)*t(X)%*%X+accept.rate*sigma.cur}
       else if (accept.rate<0.3 &accept.rate!=0){
         sigma.cur=(1-accept.rate)*summary(GLM)$cov.unscaled+accept.rate*sigma.cur}
       print(paste("accept-ratio=",round(accept.rate,3),sep=""))
       }
      draws[[x]]=par.update(beta.cur,z.cur,sigma.cur)
      }
######We can abandon estimates in draws()

######After burnin process in draws() we fix covariance matrix for proposal density
######and we continue to sample from posterior distribution
######Moredraws are the estimate we use.

      Moredraws=list()
      Moredraws[[1]]=draws[[burnin]]
      for (T in 2:niter){
       beta.cur=Moredraws[[T-1]][[1]]
       z.cur=Moredraws[[T-1]][[2]]
       Moredraws[[T]]=par.update(beta.cur,z.cur,sigma.cur) 
       if(T  %% print.every == 0){
       M.num.a=lapply((T-print.every+1):(T-1),
               function(TT){Moredraws[[TT]][[3]]})
       M.accept.rate=mean(unlist(M.num.a))
       print(paste("update-accept-ratio=",round(M.accept.rate,3)))}
      }

      return(Moredraws)


} #end of probit_mcmc_cpu


##use mini data to test probit_mcmc_cpu()

Data=read.table("mini_data.txt",fill=TRUE,header=TRUE,sep = "")
Par=read.table("mini_pars.txt",fill=TRUE,header=TRUE,sep = "")

y=as.matrix(Data[,1])
X=as.matrix(Data[,2:9])
Sigma.0.inv=diag(8)
beta.0=matrix(0,nrow=8)

probit_result=probit_mcmc_cpu(y,X,beta.0,Sigma.0.inv,
                   niter=2000,burnin=500,
                   print.every=2000,retune=100)

system.time(
    probit_mcmc_cpu(y,X,beta.0,Sigma.0.inv,
                   niter=2000,burnin=500,
                   print.every=2000,retune=100)
)


