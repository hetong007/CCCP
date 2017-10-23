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

scale.X = function(X) {
    mns = colMeans(X)
    std = sqrt(colMeans(X^2) - mns^2)
    res = t((t(X) - mns)/std)
    return(res)
}

# Data generation
set.seed(1024)

n = 100
p = 3000
beta = c(3, 1.5, 0, 0, 2, rep(0, p-5))
X = get.ARX(n, p)
X = scale.X(X)
Y = get.linearY(X, n, beta)

# Select lambda with HBIC for cccp, 2013
res = HBIC(X, Y, c(1,2,3,4,5), method="cccp")
ind = which.min(res$hbic.val)
bhat_best = res$bhat[[ind]]
mse(beta, bhat_best)

# Select lambda with HBIC for ncv, 2011
res = HBIC(X, Y, c(1,2,3,4,5,6), method="ncv")
ind = which.min(res$hbic.val)
bhat_best = res$bhat[[ind]]
mse(beta, bhat_best)

# Select lambda with cross validation for vanilla lasso
res = lasso.cv(X, Y, c(0, 0.1, 0.2, 0.3, 0.4, 0.8), 5, seed = 1024)
ind = which.min(res$rmse)
bhat_lasso = lasso.cd(X, Y, res$lambda[ind])
mse(beta, bhat_lasso)
