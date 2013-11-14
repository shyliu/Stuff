
mini <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:
s=5           # number of sampled subset 
r=50          # number of Monte Carlo iterations

#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

# I/O specifications:
datapath <- "/home/pdbaines/data"
outpath <- "output/"

# mini or full?
if (mini){
	rootfilename <- "blb_lin_reg_mini"
} else {
	rootfilename <- "blb_lin_reg_data"
}

# Filenames:
myfile=paste(rootfilename,".desc",sep="")

# Set up I/O stuff:
mypath=file.path("/home/pdbaines/data",myfile)

# Attach big.matrix :
Data=attach.big.matrix(mypath)

# Remaining BLB specs:
Num=sim_num-1000
n=nrow(Data)
Gamma=0.7
b=round(n^Gamma, digits = 0)  # subset size

FindRS=function(number){
  if (number %% 50 ==0){
  S=number  %/% 50 
  R=50
  }
  if (number %% 50 !=0){
  S=number  %/% 50 + 1
  R=number  %% 50 
  }
  return(c(S,R))
}

s_index=FindRS(Num)[1]
r_index=FindRS(Num)[2]

# Extract the subset:
set.seed(s_index)
Sub=sample(n,b,replace = FALSE)
SubSample=Data[sort(Sub),]
SubSample=data.frame(SubSample)

# Reset simulation seed:
set.seed(sim_num)

# Bootstrap dataset:
W=rmultinom(1, size = n, prob = rep(1/b,b))

# Fit lm:
#D = dim(SubSample)
#NAME = names(SubSample)

#fit=lm(SubSample[,ncol(SubSample)]~SubSample[,1:(ncol(SubSample)-1)]-1,
#      weights=W)
FIT=lm(X1001~.-1,data=SubSample,weights=W)
BETA=matrix(FIT$coefficients)

# Output file:
outfile=paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")

# Save estimates to file:

write.table(BETA,file=outfile,quote=FALSE,
            row.names=FALSE,col.names=TRUE) 


