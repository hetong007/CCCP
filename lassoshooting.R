soft.thred = function(x, tau) {
    sign(x)*max(abs(x)-tau, 0)
}

# Normal Lasso with Coordinate Descent
lasso.cd = function(x, y, lambda, max.step = 10000, thred = 0.001, verbose = FALSE) {
    n = nrow(x)
    p = ncol(x)
    if (length(y) != n)
        stop("Wrong length of y and x.")
    
    beta = rep(0, p)
    txx = colSums(x^2)
    xni.sum = rep(0, n)
    
    step = 1
    current.diff = Inf
    while (step <= max.step && thred < current.diff) {
        current.diff = 0
        for (i in 1:p) {
            up = t(x[,i,drop = FALSE]) %*% (y - xni.sum + x[,i]*beta[i])
            new.val = soft.thred(up/txx[i], lambda*n/txx[i])[1,1]
            
            val.diff = new.val - beta[i]
            current.diff = current.diff + abs(val.diff)
            
            beta[i] = new.val
            xni.sum = xni.sum + x[,i]*val.diff
        }
        if (verbose)
            cat(step, current.diff, '\n')
        step = step + 1
    }
    
    return(beta)
}

scad.loss = function(beta, lambda, gamma) {
    theta = abs(beta)
    if (theta <= lambda) {
        return(lambda*theta)
    }
    if (theta > gamma*lambda) {
        res = lambda^2*(gamma^2-1)/2/(gamma-1)
        return(res)
    }
    res = (lambda*theta*gamma-0.5*(theta^2+lambda^2))/(gamma-1)
    return(res)
}

lasso.cd.cccp = function(x, y, lambda, beta.init,
                         max.step = 10000, thred = 0.001, verbose = FALSE) {
    n = nrow(x)
    p = ncol(x)
    if (length(y) != n)
        stop("Wrong length of y and x.")
    
    beta = beta.init
    txx = colSums(x^2)
    xni.sum = rep(0, n)
    param.a = 3.7
    
    step = 1
    current.loss = Inf
    last.loss = Inf
    loss.diff = Inf
    while (step <= max.step && thred < loss.diff) {
        #if (step > 100)
            #browser()
        last.loss = current.loss
        current.loss = 0
        for (i in 1:p) {
            up = t(x[,i,drop = FALSE]) %*% (y - xni.sum + x[,i]*beta[i])
            if (abs(beta[i]) <= lambda) {
                lambda.spec = lambda
            } else {
                lambda.spec = max(param.a*lambda-abs(beta[i]), 0)/(param.a-1)
            }
            
            new.val = soft.thred(up/txx[i], lambda.spec*n/txx[i])[1,1]
            
            val.diff = new.val - beta[i]
            beta[i] = new.val
            xni.sum = xni.sum + x[,i]*val.diff
            
            val.loss = scad.loss(beta[i], lambda, param.a)
            current.loss = current.loss + val.loss
        }
        
        current.loss = current.loss + mean((y-x %*% beta)^2)/2
        if (verbose) {
            cat(step, last.loss - current.loss, '\n')
        }
        
        loss.diff = abs(current.loss - last.loss)
        step = step + 1
    }
    
    return(beta)
}

cccp.solve = function(x, y, lambda, max.step = 10000, thred = 0.001) {
    n = nrow(x)
    if (n <= 0)
        stop("Wrong input X")
    bhat_1 = lasso.cd(x, y, lambda/log(n), max.step, thred, F)
    bhat_res = lasso.cd.cccp(x, y, lambda, bhat_1, max.step, thred, F)
}

hbic.calc = function(x, y, lambda, beta) {
    n = nrow(x)
    p = ncol(x)
    if (n <= 1)
        stop("Wrong Input")
    
    res = log(mean((y-x %*% beta)^2)) + 
        sum(beta != 0) * log(log(n)) * log(p) / n
    return(res)
}

HBIC = function(x, y, lambda.vec) {
    len = length(lambda.vec)
    res = list()
    hbic.val = rep(0, len)
    
    for (i in 1:len) {
        lambda = lambda.vec[i]
        cat(i, lambda, '\t')
        bhat = cccp.solve(x, y, lambda, max.step = 10000, thred = 0.001)
        
        hbic.val[i] = hbic.calc(x, y, lambda, bhat)
        res[[i]] = bhat
        cat(hbic.val[i],'\n')
    }
    
    return(list(bhat = res, hbic.val = hbic.val))
}

mse.calc = function(beta,bhat) {
    sum((beta-bhat)^2)
}

MSE = function(x, y, lambda.vec, beta) {
    len = length(lambda.vec)
    res = list()
    hbic.val = rep(0, len)
    
    for (i in 1:len) {
        lambda = lambda.vec[i]
        cat(i, lambda, '\t')
        bhat = cccp.solve(x, y, lambda, max.step = 10000, thred = 0.001)
        
        hbic.val[i] = mse.calc(beta, bhat)
        res[[i]] = bhat
        cat(hbic.val[i],'\n')
    }
    
    return(list(bhat = res, hbic.val = hbic.val))
}
