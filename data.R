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
            X[,j] = X[,j-1]*0.6 + 0.8*rnorm(n)
        }
    }
    return(X)
}

get.linearY<-function(X, n, beta)
{
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

# Store Result 
cccp.res = matrix(0, 100, 10)
colnames(cccp.res) = c("Lambda.hbic","HBIC.hbic","MSE.hbic","TP.hbic","FP.hbic",
                       "Lambda.mse","HBIC.mse","MSE.mse","TP.mse","FP.mse")
ncv.res = cccp.res
lasso.res = cccp.res
colnames(lasso.res) = c("Lambda.rmse","HBIC.rmse","MSE.rmse","TP.rmse","FP.rmse",
                        "Lambda.mse","HBIC.mse","MSE.mse","TP.mse","FP.mse")


lambda.vec.cccp = c(0, 0.1, 0.2, 0.3, 0.4, 0.8, 1,2,3,4,5)
lambda.vec.ncv = c(0, 0.1, 0.2, 0.3, 0.4, 0.8, 1,2,3,4,5)
lambda.vec.lasso = c(0, 0.1, 0.2, 0.3, 0.4, 0.8, 1, 2)

for (i in 1:10) {
    # Data generation
    cat('Dataset No.', i, '\n')
    set.seed(i)
    
    n = 100
    p = 3000
    beta = c(3, 1.5, 0, 0, 2, rep(0, p-5))
    X = get.ARX(n, p)
    X = scale.X(X)
    Y = get.linearY(X, n, beta)
    
    # Select lambda with HBIC for cccp, 2013
    res_cccp = HBIC(X, Y, lambda.vec.cccp, method="cccp", beta)
    ind_hbic = which.min(res_cccp$hbic.val)
    ind_mse = which.min(res_cccp$mse.val)
    cccp.res[i,] = c(lambda.vec.cccp[ind_hbic], res_cccp$hbic.val[ind_hbic],
                     res_cccp$mse.val[ind_hbic],
                     res_cccp$tp.val[ind_hbic], res_cccp$fp.val[ind_hbic],
                     lambda.vec.cccp[ind_mse], res_cccp$hbic.val[ind_mse],
                     res_cccp$mse.val[ind_mse],
                     res_cccp$tp.val[ind_mse], res_cccp$fp.val[ind_mse])
    
    # Select lambda with HBIC for ncv, 2011
    res_ncv = HBIC(X, Y, lambda.vec.ncv, method="ncv", beta)
    ind_hbic = which.min(res_ncv$hbic.val)
    ind_mse = which.min(res_ncv$mse.val)
    ncv.res[i,] = c(lambda.vec.ncv[ind_hbic], res_ncv$hbic.val[ind_hbic],
                    res_ncv$mse.val[ind_hbic],
                    res_ncv$tp.val[ind_hbic], res_ncv$fp.val[ind_hbic],
                    lambda.vec.ncv[ind_mse], res_ncv$hbic.val[ind_mse],
                    res_ncv$mse.val[ind_mse],
                    res_ncv$tp.val[ind_mse], res_ncv$fp.val[ind_mse])
    
    # Select lambda with cross validation for vanilla lasso
    res_lasso = lasso.cv(X, Y, lambda.vec.lasso, beta, 5)
    ind_rmse = which.min(res_lasso$rmse)
    ind_mse = which.min(res_lasso$mse)
    lasso.res[i,] = c(res_lasso$lambda[ind_rmse], res_lasso$rmse[ind_rmse], 
                      res_lasso$mse[ind_rmse],
                      res_lasso$tp.val[ind_rmse], res_lasso$fp.val[ind_rmse],
                      res_lasso$lambda[ind_mse], res_lasso$rmse[ind_mse], 
                      res_lasso$mse[ind_mse],
                      res_lasso$tp.val[ind_rmse], res_lasso$fp.val[ind_rmse])
}

simulation.res.10 = list(CCCP = cccp.res[1:10,],
                         NCV = ncv.res[1:10,],
                         Lasso = lasso.res[1:10,])

save(simulation.res.10, file = 'simulation.10.rda')

# simulation.res.100 = list(CCCP = cccp.res,
#                          NCV = ncv.res,
#                          Lasso = lasso.res)
# save(simulation.res.100, file = 'simulation.100.rda')
