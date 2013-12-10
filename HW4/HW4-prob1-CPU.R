##Trunc_Normal() is a function to generate samples form 
##truncated normal distribution

Trunc_Normal=function(N,mu,sigma,a,b,Maxtries){
     x=matrix(0,nrow=N,ncol=1)
     for (i in 1:N){
      numtries=0
      accept=0    
      while(accept==0 & numtries<Maxtries){
      numtries = numtries+1
      cand=rnorm(1, mean = mu, sd = sigma)
         if (cand >=a & cand <=b){
         x[i,1]=cand
         accept=1
         } ##end of if (cand >=a & cand <=b)
      } ##end of while loop (start Robert sampling if necessary)

      while (accept == 0){ ## start rejection sampling
        sta = (a-mu)/sigma           
        stb = (b-mu)/sigma   
        if (sta>0 & stb>0){imu = sta}   
        else {imu = -stb}    
        alpha = (imu + sqrt(imu^2+4))/2
        testnum = 0      
        maxtestnum = 1000
        testaccept = 0
        while (testaccept==0 & testnum<maxtestnum){
        testnum = testnum + 1
        z = imu - log(runif(1)/alpha)               
        if (imu < alpha){phiz=exp(-((alpha-z)^2/2))}    
	    else {phiz=exp(-((alpha-z)^2/2)+(imu-alpha)^2/2)}         
        testmu=runif(1)   
         if (testmu<phiz){
           testaccept=1                 
           accept = 1
           if (sta>0 & stb>0){
              x[i,1]=mu+sigma*z} 
           else {
              x[i,1]=mu-sigma*z}
         } #end of if (testmu<phiz)	                  
        } #end of while (testaccept==0 & testnum<maxtestnum)

      } ##end of Robert sampling
     } ##end of for loop
     return(x)
}

val=Trunc_Normal(10000,0,1,-Inf,-10,2000)
system.time(Trunc_Normal(10000,0,1,-Inf,-10,2000))

