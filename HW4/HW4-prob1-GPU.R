library(RCUDA)
library(MASS)

m = loadModule("test.ptx")
k = m$rtruncnorm_kernel
cat("done. Setting up miscellaneous stuff...\n")

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


# Normal parameters:
N = 1000000L # 1e6L fails on my mac... :/
verbose = FALSE
vals = matrix(0,N,1)
mu = matrix(0,nrow=N,ncol=1)
sigma = matrix(1,nrow=N,ncol=1)
rng_a = 33L
rng_b = 7L
rng_c = 28L
lo = matrix(-1000000000,nrow=N,ncol=1)
hi = matrix(-10,nrow=N,ncol=1)
maxtries = 2000L

sub=.cuda(k, "vals"=vals, N, mu, sigma,rng_a,rng_b,rng_c,lo,hi,maxtries, gridDim=grid_dims, blockDim=block_dims,outputs="vals")
#save(sub,file="Probit_sub.RData")

TN_GPU_Time=list()
for (t in 1:8){
N=as.integer(10^t)
vals = matrix(0,N,1)
mu = matrix(2,nrow=N,ncol=1)
sigma = matrix(1,nrow=N,ncol=1)
rng_a = 10L
rng_b = 1L
rng_c = 5L
bg <- compute_grid(N)
grid_dims=bg$grid_dims
block_dims=bg$block_dims
lo = matrix(0,nrow=N,ncol=1)
hi = matrix(1.5,nrow=N,ncol=1)
maxtries = 1000L
TN_GPU_Time[[t]] <- system.time({
    .cuda(k, "vals"=vals, N, mu, sigma,rng_a,rng_b,rng_c,lo,hi,maxtries, gridDim=grid_dims, blockDim=block_dims,outputs="vals")
    cudaDeviceSynchronize()
})
}



