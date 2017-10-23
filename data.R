source('lassoshooting.R')

### helper function to generate data matrices X and Y
get.ARX<-function(n, p)
{
    X = matrix(0, nrow = n, ncol = p)
    for (j in 1:p)
    {
        if (j==1)
            X[,j] = rnorm(n)
        else{
            X[,j] = X[,j-1]*0.8 + 0.6*rnorm(n)
        }
    }
    return(X)
}

get.linearY<-function(X, n, beta)
{
    p = length(beta)
    eps = rnorm(n, sd = 2)
    
    Y = X%*%beta + eps
    return(as.vector(Y))
}

# helper function to calculate MSE of beta
mse = function(beta,bhat) {
    sum((beta-bhat)^2)
}

# Data generation
set.seed(1024)

n = 100
p = 1000
beta = c(3, 1.5, 0, 0, 2, rep(0, p-5))
X = get.ARX(n, p)
Y = get.linearY(X, n, beta)

# Select lambda with HBIC
res = HBIC(X, Y, c(3,3.5,4,4.5,5))
ind = which.min(res[[2]])
bhat_best = res[[1]][[ind]]
mse(beta, bhat_best)

# Calculate MSE with different lambda settings
res2 = MSE(X, Y, c(1,2,3,4,5), beta)