############### FDR CONTROL FOR CHANGE-POINTS VIA MIRROR #########################

#######################################################################################
######################### DEFINE SOME USEFUL FUNCTION #################################
#######################################################################################


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


Performance.fun <- function(sel, ge, h = NULL, change_points){
  
  num <- length(sel)
  
  if(num == 0){
    fdr <- 0
    power <- 0
  }else{
    tab <- bandwidth_fun(change_points, n = n, h = h)[change_points %in% sel, ]
    tab <- matrix(tab, ncol = 2)
    sel_c <- NULL
    for(i in 1:dim(tab)[1]){
      sel_c <- c(sel_c, tab[i,1] : tab[i,2])
    }
    
    TP <- sum(ge %in% sel_c)
    fdr <- (num -  TP)/ num
    power <- TP / length(ge) 
  }
  return(list(fdr = fdr, power = power))
}

# summation statistic function
sum_statistic <- function(x){
  n <- length(x)
  if(n%%2 !=0){
    n = n -1
  }
  res <- sum(x[1:(n/2)]) - sum(x[(n/2 +1) : n])
  res <- res/sqrt(n/2)
  return(as.numeric(abs(res)))
}


#######################################################################################
########################### DATA GENERATING FUNCTION ##################################
#######################################################################################

# n: sample size
# p: feature size
# bw: the distance between two change-points
# a: maximum change magnitude
# df: distribution
# sparsity: how many position changes 
Data_t <- function(n, p, bw, a, sparsity, df){
  
  true_position <- seq.int(from = 0, to = n, by = bw)
  m_c <- matrix(nrow = p, ncol = n)
  res_mean <- matrix(rep(0,p), ncol = 1)
  
  for(i in 1:(length(true_position)- 1)){
    start <- true_position[i] + 1
    end <- true_position[i+1]
    res_change <- sample(p, sparsity)
    res_mean[res_change,] <- res_mean[res_change,] + a
    m_c[,start:end] <- res_mean
    a <- -a 
  }
  
  d_c<- m_c + matrix(rt(n*p, df = df), nrow = p, ncol = n)
  return(list(x = d_c, true_position = true_position[-c(1, length(true_position))]))
}


Data_mvt <- function(n, p, bw, a, sparsity, df, phi){
  
  true_position <- seq.int(from = 0, to = n, by = bw)
  m_c <- matrix(nrow = p, ncol = n)
  res_mean <- matrix(rep(0,p), ncol = 1)
  
  for(i in 1:(length(true_position)- 1)){
    start <- true_position[i] + 1
    end <- true_position[i+1]
    res_change <- sample(p, sparsity)
    res_mean[res_change,] <- res_mean[res_change,] + a
    m_c[,start:end] <- res_mean
    a <- -a 
  }
  
  t <- seq.int(1,p, by = 1)
  Sigma <- phi^abs(outer(t,t,"-"))
  Eps <- mvtnorm::rmvt(n, sigma = Sigma, df = df)
  d_c<- m_c + t(Eps)
  return(list(x = d_c, true_position = true_position[-c(1, length(true_position))]))
}

Data_norm <- function(n, p, bw, a, sparsity, sd = 1){
  
  true_position <- seq.int(from = 0, to = n, by = bw)
  m_c <- matrix(nrow = p, ncol = n)
  res_mean <- matrix(rep(0,p), ncol = 1)
  
  for(i in 1:(length(true_position)- 1)){
    start <- true_position[i] + 1
    end <- true_position[i+1]
    res_change <- sample(p, sparsity)
    res_mean[res_change,] <- res_mean[res_change,] + a
    m_c[,start:end] <- res_mean
    a <- -a 
  }
  
  d_c<- m_c + matrix(rnorm(n*p, sd = sd), nrow = p, ncol = n)
  return(list(x = d_c, true_position = true_position[-c(1, length(true_position))]))
}


#######################################################################################
################################ FDR CONTROL FUNCTIONS ################################
#######################################################################################


################################### Bootstrap FDR ##########################
# According bootstrap method to get mirror statistic T^\star
# cutoff: removal parameter
Bootstrap <- function(X, cutoff){
  p <- dim(X)[1]
  n <- dim(X)[2]
  #p = 1, use drop = FALSE to aviod matrix becoming vector
  if(p == 1){
    Z_c <- c()
    for(s in (cutoff+1):(n-cutoff)){
      bar_X_neg <- sum(X[1:s])/s
      bar_X_pos <- sum(X[,(s+1):n])/(n-s)
      
      e1 <- rnorm(s)
      E1 <- X[1:s] - bar_X_neg 
      E1 <- sum(E1*e1)
      left <- sqrt((n-s)/(n*s)) * E1
      
      e2 <- rnorm(n-s)
      E2 <- X[(s+1):n] - bar_X_pos
      E2 <- sum(E2*e2)
      right <- sqrt(s/(n*(n-s))) * E2
      
      res_Z <- left - right
      Z_c <-  c(Z_c, res_Z)
    }
  }else{
    Z_c <- matrix(nrow = p, ncol = n - 2*cutoff)
    for(s in (cutoff+1):(n-cutoff)){
      bar_X_neg <- rowSums(X[,1:s])/s
      bar_X_pos <- rowSums(X[,(s+1):n])/(n-s)
      e1 <- rnorm(s)
      E1 <- sweep(X[,1:s], MARGIN = 1, bar_X_neg, FUN = '-')
      E1 <- matrix(rowSums(E1 %*% diag(e1)), ncol = 1)
      left <- sqrt((n-s)/(n*s)) * E1
      
      e2 <- rnorm(n-s)
      E2 <- sweep(X[,(s+1):n], MARGIN = 1, bar_X_pos, FUN = '-')
      E2 <- matrix(rowSums(E2 %*% diag(e2)), ncol = 1)
      right <- sqrt(s/(n*(n-s))) * E2
      
      res_Z <- left - right
      Z_c[,s-cutoff] <- matrix(res_Z, ncol = 1)
    }
  }
  return(max(abs(Z_c)))
}

#sweep(matrix(c(1,1,2,1,1,1), ncol = 2), MARGIN = 1, c(1,2,1), FUN = '-')

#matrix(c(1,1,2,1,1,1), ncol = 2) %*% diag(c(1,2))
#TRUE cusum statistic, T


cusum <- function(X, cutoff){
  p <- dim(X)[1]
  n <- dim(X)[2]
  
  Z_c <- matrix(ncol = n - 2*cutoff, nrow = p)
  for(s in (cutoff+1):(n-cutoff)){
    bar_X_neg <- matrix(rowSums(X[,1:s, drop = FALSE]), ncol = 1)/s
    bar_X_pos <- matrix(rowSums(X[,(s+1):n, drop = FALSE]), ncol = 1)/(n-s)
    res_Z <-  sqrt((s*(n-s))/(n)) * (bar_X_neg  - bar_X_pos)
    Z_c[,s-cutoff] <- matrix(res_Z, ncol = 1)
  }
  return(max(abs(Z_c)))
}


#Bootstrap(X, cutoff = 10)
#cusum(X, cutoff = 10)


Bootstrap_fdr <- function(X, change_points, fdr, cutoff, h = NULL){
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
    #print(i)
    res_interval <- interval_fdr[i,]
    res_O <- O[,res_interval[1] : res_interval[2], drop = FALSE]
    res_E <- E[,res_interval[1] : res_interval[2], drop = FALSE]
    T1 <- cusum(res_O, cutoff = cutoff)
    T2 <- Bootstrap(res_E, cutoff = cutoff)
    W_c[i] <- T1 - T2
  }
  # choose thresholds
  thre <- knockoff::knockoff.threshold(W_c, fdr = fdr)
  return(change_points[which(W_c >= thre)])
}

################### MOPS FDR ##########################
MOPS <- function(X, changepos, fdr){
  change_points <- changepos
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  changepos <- floor(changepos/2)
  
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
  
  if(p == 1){
    for(i in 1:K){
      idx_neg <- (tau[i]+1) : tau[i+1] 
      idx_pos <- (tau[i+1]+1) : tau[i+2]
      n_neg <- length(idx_neg)
      n_pos <- length(idx_pos)
      
      SE_neg <- sum(E[idx_neg])
      SO_neg <- sum(O[idx_neg])
      SE_pos <- sum(E[idx_pos])
      SO_pos <- sum(O[idx_pos])
      
      #resE <- E[(tau[i]+1): tau[i+2]]
      #Omega_inv <- 1 / ((1/2*(length(resE) - 1)) * sum((diff(resE))^2))
      Omega_inv <- 1
      W_c[i] <-  as.numeric(((n_neg * n_pos)/(n_neg + n_pos)) * (SE_neg - SE_pos) * Omega_inv %*% (SO_neg - SO_pos))
    }
  }
  else{
    for(i in 1:K){
      idx_neg <- (tau[i]+1) : tau[i+1] 
      idx_pos <- (tau[i+1]+1) : tau[i+2]
      n_neg <- length(idx_neg)
      n_pos <- length(idx_pos)
      
      SE_neg <- matrix(rowSums(E[,idx_neg]), ncol = 1)
      SO_neg <- matrix(rowSums(O[,idx_neg]), ncol = 1)
      SE_pos <- matrix(rowSums(E[,idx_pos]), ncol = 1)
      SO_pos <- matrix(rowSums(O[,idx_pos]), ncol = 1)
      
      #resE <- E[,(tau[i]+1): tau[i+2]]
      
      #S_diff <- resE[,1:(dim(resE)[2] - 1)] - resE[,2:(dim(resE)[2])]
      #Omega_inv <- solve((1/2*(dim(resE)[2] - 1)) * S_diff %*% t(S_diff))
      Omega_inv <- diag(p)
      W_c[i] <-  as.numeric(((n_neg * n_pos)/(n_neg + n_pos)) * t(SE_neg - SE_pos) %*% Omega_inv %*% (SO_neg - SO_pos))
    }
  }
  
  thre <- knockoff::knockoff.threshold(W_c, fdr = fdr)
  
  return(change_points[W_c >= thre])
}



###################### SARA FDR ###########################
SaRa <- function(X, change_points, fdr, h = NULL, sigma){
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  K <- length(change_points)
  n <- dim(X)[2] 
  p <- dim(X)[1]
  if(p == 1){
    # p = 1, no projection and no data-splitting
    x_proj <- as.vector(X)
    interval_fdr <- bandwidth_fun(change_points, n = n, h = h)
    
  }else{
    
    if(n%%2 != 0){n = n - 1} 
    
    ## data splitting
    tt <- 1:n
    index_O <- tt[1:n %% 2 == 1]
    index_E <- tt[1:n %% 2 == 0]
    O <- X[,index_O]
    E <- X[,index_E]
    
    n_O <- length(index_O)
    
    ## projection in Odd data
    interval_m <- bandwidth_fun(floor(change_points/2), n = n_O, h = NULL)
    
    x_proj <- rep(0, n_O)
    for(i in 1:K){
      res_interval <- interval_m[i,]
      res_X <- O[,res_interval[1] : res_interval[2]]
      res_proj <- InspectChangepoint::locate.change(res_X)$vector.proj
      res_proj <- matrix(res_proj, ncol = 1)
      res_pv <- t(E[,res_interval[1] : res_interval[2]]) %*% res_proj
      x_proj[res_interval[1] : res_interval[2]] <- as.vector(res_pv)
    }
    
    interval_fdr <- bandwidth_fun(floor(change_points/2), n = n_O, h = h)
    
  }
  
  ### FDR control procedure
  
  W_c <- NULL
  for(i in 1:K){
    #print(i)
    res_interval <- interval_fdr[i,]
    res_x <- x_proj[res_interval[1] : res_interval[2]]/sigma
    
    T1 <- sum_statistic(res_x)
    W_c[i] <- 2*pnorm(T1, lower.tail = FALSE)
  }
  # choose thresholds
  sel <- which(p.adjust(W_c, method = 'BY') <= 0.2)
  return(change_points[sel])
}


##########################################################################################
######################################### Oracle #########################################
##########################################################################################
my_rnorm <- function(n,p, df = NULL, sigma = 1, phi = NULL){
  return(matrix(rnorm(n*p, sd = sigma), nrow = p, ncol = n))
}

my_rt <- function(n,p,df, sigma = NULL, phi=NULL){
  return(matrix(rt(n*p, df = df), nrow = p, ncol = n))
}

#n = 10; p = 3; df = 7


my_elliptical <- function(n, p, df, sigma = NULL, phi){
  t <- seq.int(1,p, by = 1)
  Sigma <- phi^abs(outer(t,t,"-"))
  return(t(mvtnorm::rmvt(n, sigma = Sigma, df = df)))
}
#my_elliptical(n,p,df)


Oracle_fdr <- function(X, change_points, fdr, cutoff, distribution = c('normal', 't', 'elliptical'), h = NULL, 
                       sigma = NULL, df = NULL, phi = NULL){
  if(!is.matrix(X)){
    stop('X should be a matrix')
  }
  
  distribution <- switch (distribution,
                          normal = my_rnorm,
                          t = my_rt,
                          elliptical = my_elliptical
  )
  
  
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
    #print(i)
    res_interval <- interval_fdr[i,]
    res_O <- O[,res_interval[1] : res_interval[2], drop = FALSE]
    #res_E <- E[,res_interval[1] : res_interval[2], drop = FALSE]
    res_E <- distribution(dim(res_O)[2], dim(res_O)[1], df = df, sigma = sigma, phi = phi)
    T1 <- cusum(res_O, cutoff = cutoff)
    T2 <- cusum(res_E, cutoff = cutoff)
    W_c[i] <- T1 - T2
  }
  # choose thresholds
  thre <- knockoff::knockoff.threshold(W_c, fdr = fdr)
  return(change_points[which(W_c >= thre)])
}

##########################################################################################
######################################### multiple #########################################
##########################################################################################
# X : data matrix
# cut : removal parameter
# B: multiple times
multiple_bootstrap <- function(X, change_points, fdr, cutoff, h = NULL, B = 500){
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
    #print(i)
    res_interval <- interval_fdr[i,]
    res_O <- O[,res_interval[1] : res_interval[2], drop = FALSE]
    res_E <- E[,res_interval[1] : res_interval[2], drop = FALSE]
    T1 <- cusum(res_O, cutoff = cutoff)
    bootstrap_value <- replicate(B, Bootstrap(res_E, cutoff = cutoff))
    Pi <- mean(bootstrap_value >= T1)
    #T2 <- Bootstrap(res_E, cutoff = cutoff)
    W_c[i] <- Pi
  }
  # choose thresholds
  #thre <- knockoff::knockoff.threshold(W_c, fdr = fdr)
  adj <- p.adjust(W_c, method = 'BH')
  
  return(change_points[which(adj <= (1:length(adj))/length(adj))])
}

