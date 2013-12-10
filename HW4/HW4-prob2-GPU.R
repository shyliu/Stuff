library(RCUDA)
library(MASS)

m = loadModule("rtruncnorm.ptx")
k = m$rtruncnorm_kernel

##compute grid&block size
compute_grid <- function(N,sqrt_threads_per_block=16L,grid_nd=1)
{
    # if...
    # N = 1,000,000
    # => 1954 blocks of 512 threads will suffice
    # => (62 x 32) grid, (512 x 1 x 1) blocks
    # Fix block dims:
    block_dims <- c(as.integer(sqrt_threads_per_block), as.integer(sqrt_threads_per_block), 1L)
    threads_per_block <- prod(block_dims)
    if (grid_nd==1){
      grid_d1 <- as.integer(max(1L,ceiling(N/threads_per_block)))
      grid_d2 <- 1L
    } else {
      grid_d1 <- as.integer(max(1L, floor(sqrt(N/threads_per_block))))
      grid_d2 <- as.integer(ceiling(N/(grid_d1*threads_per_block)))
    }
    grid_dims <- c(grid_d1, grid_d2, 1L)
    return(list("grid_dims"=grid_dims,"block_dims"=block_dims))
}

###### start probit_mcmc_gpu()

probit_mcmc_gpu <- function (y,X,beta.0,Sigma.0.inv,
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

      z.generate=function(beta.input,A,B,C){
       N=as.integer(nrow(y))
       vals = matrix(0,N,1)
       mu = X%*%beta.input
       sigma = matrix(1,nrow=N,ncol=1)
       rng_a = as.integer(A)
       rng_b = as.integer(B)
       rng_c = as.integer(C)
       bg <- compute_grid(N)
       grid_dims=bg$grid_dims
       block_dims=bg$block_dims
       maxtries = 1000L
       lo = matrix(0,nrow=N,ncol=1) 
       hi = matrix(0,nrow=N,ncol=1)
       for (i in 1:N){
         if (y[i,1]==0){lo[i,1]=-100000L
            hi[i,1]=0L}
         else {lo[i,1]=0L
               hi[i,1]=100000L}
       } #end of for loop
      Z.sample=.cuda(k,"vals"=vals,N,mu,sigma,rng_a,rng_b,rng_c,lo,hi,maxtries,gridDim=grid_dims,blockDim=block_dims,outputs="vals")     
      return(Z.sample)                 
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

    par.update=function(beta.input,z.input,sigma.input, A.input, B.input, C.input){
        Z.par.update=z.generate(beta.input,A.input, B.input, C.input)
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
      z.cur=z.generate(beta.cur, burnin+1, burnin+2, burnin+3)
      draws[[1]]=list(beta.cur,z.cur,0)
      sigma.cur=summary(GLM)$cov.unscaled
      for (x in 2:burnin){
      beta.cur=draws[[x-1]][[1]]
      z.cur=draws[[x-1]][[2]]
      if (x %% retune == 0){
       num.a=lapply((x-retune+1):(x-1),function(tt){draws[[tt]][[3]]})
       accept.rate=mean(unlist(num.a))
       Temp_cov=lapply((x-retune+1):(x-1),function(tt){draws[[tt]][[1]]})
       Cov=cov(t(matrix(unlist(Temp_cov),nrow=ncol(X))))
       if (accept.rate>=0.5){
         sigma.cur=accept.rate*t(X)%*%X+(1-accept.rate)*sigma.cur}
       else if (accept.rate<0.3 &accept.rate!=0){
         sigma.cur=(1-accept.rate)*Cov%*%t(X)%*%X+accept.rate*sigma.cur}
       else if (accept.rate==0){sigma.cur=diag(ncol(X))}

       print(paste("accept-ratio=",round(accept.rate,3),sep=""))
       }
      A.update=2*x
      B.update=3*x
      C.update=x
      draws[[x]]=par.update(beta.cur,z.cur,sigma.cur, A.update, B.update, C.update)
      }

######After burnin process in draws() we fix covariance matrix for proposal density
######and we continue to sample from posterior distribution
######Moredraws are the estimate we use.

      Moredraws=list()
      Moredraws[[1]]=draws[[burnin]]
      for (T in 2:niter){
       beta.cur=Moredraws[[T-1]][[1]]
       z.cur=Moredraws[[T-1]][[2]]
       A.Mupdate=2*T
       B.Mupdate=3*T
       C.Mupdate=T
       Moredraws[[T]]=par.update(beta.cur,z.cur,sigma.cur, A.update, B.update, C.update) 
       if(T  %% print.every == 0){
       M.num.a=lapply((T-print.every+1):(T-1),
               function(TT){Moredraws[[TT]][[3]]})
       M.accept.rate=mean(unlist(M.num.a))
       print(paste("update-accept-ratio=",round(M.accept.rate,3)))}
      }

      return(Moredraws)


} #probit_mcmc_gpu

Data=read.table("mini_data.txt",fill=TRUE,header=TRUE,sep = "")
Par=read.table("mini_pars.txt",fill=TRUE,header=TRUE,sep = "")

y=as.matrix(Data[,1])
X=as.matrix(Data[,2:9])
Sigma.0.inv <- diag(8)
beta.0=matrix(0,nrow=8)

probit_result<- probit_mcmc_gpu(y,X,beta.0,Sigma.0.inv,niter=2000,burnin=500,print.every=2000,retune=100)

GPU_Time= system.time(
    probit_mcmc_gpu(y,X,beta.0,Sigma.0.inv,niter=2000,burnin=500,print.every=2000,retune=100)
)
save(GPU_Time,file="Probit_MCMC_GPU_Time_mini.RData")


