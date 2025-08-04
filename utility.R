#######################################################################################
###################################data generating new ################################
## Mean change model
# n: sample size
# p: dimension
# cps: change points location
# a: amplitude
# dis: error distribution, including three modes: c('normal(0,1)', 't(3)/\sqrt{3}', '(chi(3)-3)/\sqrt{6}') 
# sparsity : number of changes in coordinates 
# rho : correlation of COV(X) (chi squre does not have this setting!)

#n = 1000
#p = 10
#cps = seq.int(100, 900, by = 100)
#a = 1
#sparsity = p
#dis = 'chi'
#rho = 0
#dense = FALSE

MCP_model <- function(n, p, cps, a, sparsity, dis = c('normal', 't', 'chi'), rho = 0, dense = FALSE, df = 0){
  if(sparsity > p){
    warning('sparsity is larger than dimension')
  }
  
  
  if(dense){
    Sigma = matrix(rep(rho, p*p), ncol = p)
    diag(Sigma) <- rep(1, p)
  }else{
    Sigma = toeplitz(rho^(0:(p-1)))
  }
  
  cps1 <- c(0, cps, n) 
  
  m_c <- matrix(nrow = p, ncol = n)
  res_mean <- matrix(rep(0,p), ncol = 1)
  
  for(i in 1:(length(cps1)- 1)){
    start <- cps1[i] + 1
    end <- cps1[i+1]
    res_change <- sample(p, sparsity)
    res_mean[res_change,] <- res_mean[res_change,] + a
    m_c[,start:end] <- res_mean
    a <- -a 
  }
  
  if(dis == 'normal'){
    Eps <- mvtnorm::rmvnorm(n, sigma = Sigma)
  }else if(dis == 't'){
    Eps <- mvtnorm::rmvt(n, sigma = Sigma, df = 10)
  }else if(dis == 'chi'){
    Eps <- (rchisq(n*p, df = 3) - 3)/sqrt(6)
    Eps <- matrix(Eps, nrow = n, ncol = p)
  }else if(dis == 'df'){
    Eps <- matrix(rt(n*p, df = df), nrow = n, ncol = p) 
  }else{
    warning('distribution setting is wrong')
  }
  d_c<- m_c + t(Eps)
  return(d_c)
}


## Mean change model
# n: sample size
# p: dimension
# cps: change points location
# a: amplitude
# dis: error distribution, including three modes: c('normal(0,1)', 't(3)/\sqrt{3}', '(chi(3)-3)/\sqrt{6}') 
# s : number of nonzero elements in regression coefficients 
# rho : correlation of COR(X)
# dense : the covariance matrix is dense or not

#n = 1000
#p = 10
#cps = seq.int(100, 900, by = 100)
#a = 1
#s = p
#dis = 'chi'
#rho = 0
#dense = FALSE



Regression_model <- function(n, p, cps, a, s, dis = c('normal', 't', 'chi'), rho = 0, dense = FALSE){
  if(s > p){
    warning('sparsity is larger than dimension')
  }
  
  if(dense){
    Sigma = matrix(rep(rho, p*p), ncol = p)
    diag(Sigma) <- rep(1, p)
  }else{
    Sigma = toeplitz(rho^(0:(p-1)))
  }
  
  X = matrix(rnorm(n*p),ncol = p) %*% chol(Sigma)
  
  y <- matrix(ncol = 1, nrow = n)
  idx_c <- c(0, cps, n)
  for(i in 1:(length(idx_c) - 1)){
    star <- idx_c[i] + 1
    end <- idx_c[i+1]
    nonzero = sample(p, s)
    beta = matrix(a * (1:p %in% nonzero), ncol = 1)
    y[star:end] <- X[star:end,] %*% beta
    a <- -a
  }
  
  if(dis == 'normal'){
    eps <- rnorm(n)
  }else if(dis == 't'){
    eps <- rt(n, df = 5)/sqrt(5)
  }else if(dis == 'chi'){
    eps <- (rchisq(n, df = 3)-3)/sqrt(6)
  }else{
    warning('distribution setting is wrong')
  }
  y = y + matrix(eps, ncol = 1)
  
  return(list(X = X, y = y))
}


#######################################################################################
################################### Evaluation function ################################

### bandwidth construction function
## change_points : change_points estimated by upstream processes
## n : number of time points
## h : user-specified bandwidth, if NULL, then maximum bandwidth will be applied

## return a matrix: each row indicates the interval for each change-point
bandwidth_fun <- function(change_points, n, h = NULL){
  K <- length(change_points)
  tau <- c(0, change_points, n)
  if (!is.null(h) && h > min(floor(diff(tau)/2))){
    warning('h is too large, set h = NULL automatically, maximum bandwidth will be applied')
    h = NULL
  }
  # maximum bandwidth
  if(is.null(h)){
    cutoff <- sapply(1:(K+1), function(t){ceiling((tau[t] + tau[t+1])/2)})
    interval_m <- matrix(ncol = 2, nrow = K)
    for( i in 1:K){
      interval_m[i,1] <-  cutoff[i]
      interval_m[i,2] <-  cutoff[i+1] - 1
    }
  }
  else{
    interval_m <- matrix(ncol = 2, nrow = K)
    for( i in 1:K){
      interval_m[i,1] <-  change_points[i] - h
      interval_m[i,2] <-  change_points[i] + h - 1
    }
  }
  return(interval_m)
}


#sel : selected change points by fdr control methods
#ge : true change points
#change_points:  selected change points by upstream model


Performance.fun <- function(sel, ge, h = NULL, change_points){
  
  num <- length(sel)
  
  if(num == 0){
    fdr <- 0
    power <- 0
    fwer <- 0
  }else{
    tab <- bandwidth_fun(change_points, n = n, h = h)[change_points %in% sel, ]
    tab <- matrix(tab, ncol = 2)
    sel_c <- NULL
    for(i in 1:dim(tab)[1]){
      sel_c <- c(sel_c, tab[i,1] : tab[i,2])
    }
    
    TP <- sum(ge %in% sel_c)
    fdr <- (num -  TP)/ num
    fwer <- (num -  TP)
    power <- TP / length(ge) 
  }
  return(list(fdr = fdr, power = power, fwer = fwer))
}


fwer_thre_fun <- function(W, k){
  if(k < 0){
    warning('k should be positive')
  }
  
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t) (sum(W <= -t)))
  ok = which(ratio <= k)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

#W <- c(1,2,2,1,-1,-1,-1,3)

################################### Bootstrap FDR ##########################
# According bootstrap method to get mirror statistic T^\star
# cutoff: removal parameter
# alpha: l_alpha norm

lp <- function(x, alpha){
  if(alpha == 'inf'){
    res <- max(abs(x))
  }else{
    res <- (sum(abs(x)^alpha))^{1/alpha}
  }
  return(res)
}

Bootstrap <- function(X, cutoff, alpha = alpha){
  p <- dim(X)[1]
  n <- dim(X)[2]
  e <- rnorm(n)
  Z_c <- matrix(nrow = p, ncol = n - 2*cutoff)
  for(s in (cutoff+1):(n-cutoff)){
    #bar_X_neg <- rowSums(X[,1:s, drop = FALSE])/s
    bar_X_neg <- rowMeans(X[,1:s, drop = FALSE])
    #bar_X_pos <- rowSums(X[,(s+1):n, drop = FALSE])/(n-s)
    bar_X_pos <- rowMeans(X[,(s+1):n, drop = FALSE])
    
    
    E1 <- X[,1:s, drop = FALSE] - matrix(rep(bar_X_neg, s), ncol = s, byrow = FALSE)
    #E1 <- sweep(X[,1:s, drop = FALSE], MARGIN = 1, bar_X_neg, FUN = '-')
    E1 <- matrix(rowSums(E1 %*% diag(e[1:s])), ncol = 1)
    left <- sqrt((n-s)/(n*s)) * E1
    
    
    E2 <- X[,(s+1):n, drop = FALSE] - matrix(rep(bar_X_pos, n-s), ncol = n-s, byrow = FALSE)
    #E2 <- sweep(X[,(s+1):n, drop = FALSE], MARGIN = 1, bar_X_pos, FUN = '-')
    E2 <- matrix(rowSums(E2 %*% diag(e[(s+1):n])), ncol = 1)
    right <- sqrt(s/(n*(n-s))) * E2
    
    res_Z <- left - right
    Z_c[,s-cutoff] <- matrix(res_Z, ncol = 1)
  }
  
  back <- apply(Z_c, 2, lp, alpha)
  return(max(back))
}


#sweep(matrix(c(1,1,2,1,1,1), ncol = 2), MARGIN = 1, c(1,2,1), FUN = '-')

#matrix(c(1,1,2,1,1,1), ncol = 2) %*% diag(c(1,2))
#TRUE cusum statistic, T


cusum <- function(X, cutoff, alpha){
  p <- dim(X)[1]
  n <- dim(X)[2]
  
  Z_c <- matrix(ncol = n - 2*cutoff, nrow = p)
  for(s in (cutoff+1):(n-cutoff)){
    #bar_X_neg <- matrix(rowSums(X[,1:s, drop = FALSE]), ncol = 1)/s
    bar_X_neg <- matrix(rowMeans(X[,1:s, drop = FALSE]), ncol = 1)
    #bar_X_pos <- matrix(rowSums(X[,(s+1):n, drop = FALSE]), ncol = 1)/(n-s)
    bar_X_pos <- matrix(rowMeans(X[,(s+1):n, drop = FALSE]), ncol = 1)
    res_Z <-  sqrt((s*(n-s))/(n)) * (bar_X_neg  - bar_X_pos)
    Z_c[,s-cutoff] <- matrix(res_Z, ncol = 1)
  }
  back <- apply(Z_c, 2, lp, alpha)
  return(max(back))
}



Bootstrap_fdr <- function(X, change_points, fdr, cutoff, alpha = 'inf', h = NULL, side_info = TRUE){
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  K <- length(change_points)
  n <- dim(X)[2] 
  p <- dim(X)[1]
  
  if(n%%2 != 0){n = n - 1} 
  
  ## data splitting
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  O <- X[,index_O, drop = FALSE]
  E <- X[,index_E, drop = FALSE]
  
  n_O <- length(index_O)
  
  interval_fdr <- bandwidth_fun(floor(change_points/2), n = n_O, h = h)
  
  W_c <- NULL
  if(side_info == TRUE){
    for(i in 1:K){
      #print(i)
      res_interval <- interval_fdr[i,]
      res_O <- O[,res_interval[1] : res_interval[2], drop = FALSE]
      res_E <- E[,res_interval[1] : res_interval[2], drop = FALSE]
      #side information
      T3 <- cusum(res_E, cutoff = cutoff, alpha = alpha)
      T1 <- cusum(res_O, cutoff = cutoff, alpha = alpha)
      #T2 <- Bootstrap(res_E, cutoff = cutoff, alpha = alpha)
      T2 <- Bootstrap(res_O, cutoff = cutoff, alpha = alpha)
      W_c[i] <- (T1 - T2)*T3
    }
  }else{
    #print(i)
    for(i in 1:K){
    res_interval <- interval_fdr[i,]
    res_O <- O[,res_interval[1] : res_interval[2], drop = FALSE]
    res_E <- E[,res_interval[1] : res_interval[2], drop = FALSE]
    T1 <- cusum(res_O, cutoff = cutoff, alpha = alpha)
    #T2 <- Bootstrap(res_E, cutoff = cutoff, alpha = alpha)
    T2 <- Bootstrap(res_O, cutoff = cutoff, alpha = alpha)
    W_c[i] <- (T1 - T2)
    }
  }
 
  # choose thresholds
  thre <- knockoff::knockoff.threshold(W_c, fdr = fdr)
  return(change_points[which(W_c >= thre)])
}




Bootstrap_fwer <- function(X, change_points, fwer, cutoff, alpha = 'inf', h = NULL, side_info = TRUE){
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  K <- length(change_points)
  n <- dim(X)[2] 
  p <- dim(X)[1]
  
  if(n%%2 != 0){n = n - 1} 
  
  ## data splitting
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  O <- X[,index_O, drop = FALSE]
  E <- X[,index_E, drop = FALSE]
  
  n_O <- length(index_O)
  
  interval_fdr <- bandwidth_fun(floor(change_points/2), n = n_O, h = h)
  
  W_c <- NULL
  if(side_info == TRUE){
    for(i in 1:K){
      #print(i)
      res_interval <- interval_fdr[i,]
      res_O <- O[,res_interval[1] : res_interval[2], drop = FALSE]
      res_E <- E[,res_interval[1] : res_interval[2], drop = FALSE]
      #side information
      T3 <- cusum(res_E, cutoff = cutoff, alpha = alpha)
      T1 <- cusum(res_O, cutoff = cutoff, alpha = alpha)
      #T2 <- Bootstrap(res_E, cutoff = cutoff, alpha = alpha)
      T2 <- Bootstrap(res_O, cutoff = cutoff, alpha = alpha)
      W_c[i] <- (T1 - T2)*T3
    }
  }else{
    #print(i)
    for(i in 1:K){
      res_interval <- interval_fdr[i,]
      res_O <- O[,res_interval[1] : res_interval[2], drop = FALSE]
      res_E <- E[,res_interval[1] : res_interval[2], drop = FALSE]
      T1 <- cusum(res_O, cutoff = cutoff, alpha = alpha)
      #T2 <- Bootstrap(res_E, cutoff = cutoff, alpha = alpha)
      T2 <- Bootstrap(res_O, cutoff = cutoff, alpha = alpha)
      W_c[i] <- (T1 - T2)
    }
  }
  
  # choose thresholds
  thre <- fwer_thre_fun(W_c, fwer)
  return(change_points[which(W_c >= thre)])
}



################### MOPS FDR ##########################
R_MOPS <- function(X, change_points, fdr, h = NULL){
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  K <- length(change_points)
  n <- dim(X)[2] 
  p <- dim(X)[1]
  
  if(n%%2 != 0){n = n - 1} 
  
  ## data splitting
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  O <- X[,index_O, drop = FALSE]
  E <- X[,index_E, drop = FALSE]
  
  n_O <- length(index_O)
  interval_fdr <- bandwidth_fun(floor(change_points/2), n = n_O, h = h)
  
  W_c <- NULL
  
  for(i in 1:K){
    res_interval <- interval_fdr[i,]
    res_mid <- floor((res_interval[2] + res_interval[1])/2)
    idx_neg <- res_interval[1] : res_mid 
    idx_pos <- (res_mid+1) : res_interval[2]
    n_neg <- length(idx_neg)
    n_pos <- length(idx_pos)
    
    SE_neg <- matrix(rowSums(E[,idx_neg, drop = FALSE]), ncol = 1)
    SO_neg <- matrix(rowSums(O[,idx_neg, drop = FALSE]), ncol = 1)
    SE_pos <- matrix(rowSums(E[,idx_pos, drop = FALSE]), ncol = 1)
    SO_pos <- matrix(rowSums(O[,idx_pos, drop = FALSE]), ncol = 1)
      
    #resE <- E[,(tau[i]+1): tau[i+2]]
      
    #S_diff <- resE[,1:(dim(resE)[2] - 1)] - resE[,2:(dim(resE)[2])]
    #Omega_inv <- solve((1/2*(dim(resE)[2] - 1)) * S_diff %*% t(S_diff))
    Omega_inv <- diag(p)
    W_c[i] <-  as.numeric(((n_neg * n_pos)/(n_neg + n_pos)) * t(SE_neg - SE_pos) %*% Omega_inv %*% (SO_neg - SO_pos))
  }
  thre <- knockoff::knockoff.threshold(W_c, fdr = fdr)
  return(change_points[which(W_c >= thre)])
}


R_MOPS_fwer <- function(X, change_points, fwer, h = NULL){
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  K <- length(change_points)
  n <- dim(X)[2] 
  p <- dim(X)[1]
  
  if(n%%2 != 0){n = n - 1} 
  
  ## data splitting
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  O <- X[,index_O, drop = FALSE]
  E <- X[,index_E, drop = FALSE]
  
  n_O <- length(index_O)
  interval_fdr <- bandwidth_fun(floor(change_points/2), n = n_O, h = h)
  
  W_c <- NULL
  
  for(i in 1:K){
    res_interval <- interval_fdr[i,]
    res_mid <- floor((res_interval[2] + res_interval[1])/2)
    idx_neg <- res_interval[1] : res_mid 
    idx_pos <- (res_mid+1) : res_interval[2]
    n_neg <- length(idx_neg)
    n_pos <- length(idx_pos)
    
    SE_neg <- matrix(rowSums(E[,idx_neg, drop = FALSE]), ncol = 1)
    SO_neg <- matrix(rowSums(O[,idx_neg, drop = FALSE]), ncol = 1)
    SE_pos <- matrix(rowSums(E[,idx_pos, drop = FALSE]), ncol = 1)
    SO_pos <- matrix(rowSums(O[,idx_pos, drop = FALSE]), ncol = 1)
    
    #resE <- E[,(tau[i]+1): tau[i+2]]
    
    #S_diff <- resE[,1:(dim(resE)[2] - 1)] - resE[,2:(dim(resE)[2])]
    #Omega_inv <- solve((1/2*(dim(resE)[2] - 1)) * S_diff %*% t(S_diff))
    Omega_inv <- diag(p)
    W_c[i] <-  as.numeric(((n_neg * n_pos)/(n_neg + n_pos)) * t(SE_neg - SE_pos) %*% Omega_inv %*% (SO_neg - SO_pos))
  }
  thre <- fwer_thre_fun(W_c, fwer)
  return(change_points[which(W_c >= thre)])
}



MOPS <- function(X, change_points, fdr){
  #change_points <- changepos
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  changepos <- floor(change_points/2)
  
  K <- length(changepos)
  n <- dim(X)[2] 
  p <- dim(X)[1]
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  O <- X[,index_O]
  E <- X[,index_E]
  if(p == 1){
    n_O <- length(O)
  }else{
    n_O <- dim(O)[2]
    
  }
  tau <- c(0, changepos, n_O)
  
  # if p == 1, all are vector
  W_c <- NULL
  
  
  for(i in 1:K){
    idx_neg <- (tau[i]+1) : tau[i+1] 
    idx_pos <- (tau[i+1]+1) : tau[i+2]
    n_neg <- length(idx_neg)
    n_pos <- length(idx_pos)
      
    SE_neg <- matrix(rowSums(E[,idx_neg, drop = FALSE]), ncol = 1)
    SO_neg <- matrix(rowSums(O[,idx_neg, drop = FALSE]), ncol = 1)
    SE_pos <- matrix(rowSums(E[,idx_pos, drop = FALSE]), ncol = 1)
    SO_pos <- matrix(rowSums(O[,idx_pos, drop = FALSE]), ncol = 1)
      
    #resE <- E[,(tau[i]+1): tau[i+2]]
      
    #S_diff <- resE[,1:(dim(resE)[2] - 1)] - resE[,2:(dim(resE)[2])]
    #Omega_inv <- solve((1/2*(dim(resE)[2] - 1)) * S_diff %*% t(S_diff))
    Omega_inv <- diag(p)
    W_c[i] <-  as.numeric(((n_neg * n_pos)/(n_neg + n_pos)) * t(SE_neg - SE_pos) %*% Omega_inv %*% (SO_neg - SO_pos))
  }
  
  thre <- knockoff::knockoff.threshold(W_c, fdr = fdr)
  
  return(change_points[which(W_c >= thre)])
}

MOPS_fwer <- function(X, change_points, fwer){
  #change_points <- changepos
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  changepos <- floor(change_points/2)
  
  K <- length(changepos)
  n <- dim(X)[2] 
  p <- dim(X)[1]
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  O <- X[,index_O]
  E <- X[,index_E]
  if(p == 1){
    n_O <- length(O)
  }else{
    n_O <- dim(O)[2]
    
  }
  tau <- c(0, changepos, n_O)
  
  # if p == 1, all are vector
  W_c <- NULL
  
  
  for(i in 1:K){
    idx_neg <- (tau[i]+1) : tau[i+1] 
    idx_pos <- (tau[i+1]+1) : tau[i+2]
    n_neg <- length(idx_neg)
    n_pos <- length(idx_pos)
    
    SE_neg <- matrix(rowSums(E[,idx_neg, drop = FALSE]), ncol = 1)
    SO_neg <- matrix(rowSums(O[,idx_neg, drop = FALSE]), ncol = 1)
    SE_pos <- matrix(rowSums(E[,idx_pos, drop = FALSE]), ncol = 1)
    SO_pos <- matrix(rowSums(O[,idx_pos, drop = FALSE]), ncol = 1)
    
    #resE <- E[,(tau[i]+1): tau[i+2]]
    
    #S_diff <- resE[,1:(dim(resE)[2] - 1)] - resE[,2:(dim(resE)[2])]
    #Omega_inv <- solve((1/2*(dim(resE)[2] - 1)) * S_diff %*% t(S_diff))
    Omega_inv <- diag(p)
    W_c[i] <-  as.numeric(((n_neg * n_pos)/(n_neg + n_pos)) * t(SE_neg - SE_pos) %*% Omega_inv %*% (SO_neg - SO_pos))
  }
  
  thre <- fwer_thre_fun(W_c, k = fwer)
  
  return(change_points[which(W_c >= thre)])
}



################################## DERANDOM ######################

DeRandom_Bootstrap_fdr <- function(X, change_points, fdr, cutoff, B = 100, alpha = 'inf', h = NULL){
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  K <- length(change_points)
  n <- dim(X)[2] 
  p <- dim(X)[1]
  
  if(n%%2 != 0){n = n - 1} 
  
  ## data splitting
  tt <- 1:n
  index_O <- tt[1:n %% 2 == 1]
  index_E <- tt[1:n %% 2 == 0]
  O <- X[,index_O, drop = FALSE]
  E <- X[,index_E, drop = FALSE]
  
  n_O <- length(index_O)
  
  interval_fdr <- bandwidth_fun(floor(change_points/2), n = n_O, h = h)
  
  W_c <- NULL
  for(i in 1:K){
    res_interval <- interval_fdr[i,]
    res_O <- O[,res_interval[1] : res_interval[2], drop = FALSE]
    res_E <- E[,res_interval[1] : res_interval[2], drop = FALSE]
    T1 <- cusum(res_O, cutoff = cutoff, alpha = alpha)
    #T2 <- Bootstrap(res_E, cutoff = cutoff, alpha = alpha)
    bootstrap_value <- replicate(B, Bootstrap(res_E, cutoff = cutoff, alpha = alpha))
    Pi <- mean(bootstrap_value >= T1)
    #T2 <- Bootstrap(res_O, cutoff = cutoff, alpha = alpha)
    W_c[i] <- Pi
    #W_c[i] <- (T1 - T2)
  }
  
  # choose thresholds
  adj <- p.adjust(W_c, method = 'BH')
  #thre <- knockoff::knockoff.threshold(W_c, fdr = fdr)
  #return(change_points[which(W_c >= thre)])
  return(change_points[which(adj <= fdr)])
}




