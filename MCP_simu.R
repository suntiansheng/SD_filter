setwd("~/bootstrap_fdr/new20221007/")
source('./utility.R')
#library(foreach)
#library(doParallel)

#cps: true change points
#h: upper error bound
#num_fds: number of false discoveries
#candidate_set_fun <- function(cps, h, num_fds){
  n <- length(cps)
  if(num_fds > n-1){
    warnings('numbers of FD are too large')
  }
  sign <- ifelse(rbinom(n, 1, 0.5) == 0, -1, 1)
  shifts = floor(runif(n) * h) * sign
  can <- cps + shifts
  fd <- sample(cps[-n] + diff(cps)/2, num_fds)
  re <- c(can, fd)
  return(sort(re))
}

#############################################################################
################################# Signal magnitude ######################
#############################################################################

set.seed(1234)
n = 4000
p = 50
cps = seq.int(200, 3800, by = 200)
rho = 0
sparsity = 1
dis = 't'
dense = FALSE
fdr = 0.2
h = NULL
cutoff = 10

#rho_c = seq(0.3, 0.7, by = 0.1)
a_c <- seq(from = 1.2, to = 2.2, by = 0.2)


change_points <- seq.int(from = 0, to = n, by = 150)
change_points <- change_points[-1]

#Performance.fun(sel = change_points, ge = cps, h = NULL, change_points = change_points)
#change_points <- candidate_set_fun(cps = cps, h = 30, num_fds = 5)

#change_points <- cps

simulation_time <- 20

iter_time <- length(a_c)

RESULT <- matrix(ncol = 7, nrow = iter_time)

for(iter in 1:iter_time){
  
  a <- a_c[iter]
  
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  RMO_fdr <- NULL
  RMO_power <- NULL
  
  for(i in 1:simulation_time){
    
    X <- MCP_model(n = n, p = p, cps = cps, dis = dis, a = a, sparsity = sparsity, rho = rho, dense = dense) 
    
    #Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = 0.2, cutoff = 5, h = NULL)
    Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 'inf', side_info = TRUE)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = cps, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    
    #Derandom_sel <- DeRandom_Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 'inf', B = 100)
    #Derandom_per <- Performance.fun(Derandom_sel, ge = cps, change_points = change_points, h = h)
    #DE_fdr[i] <- Derandom_per$fdr
    #DE_power[i] <- Derandom_per$power
    
    MOPS_sel <- MOPS(X = X, change_points = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = cps, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    R_MOPS_sel <- R_MOPS(X = X, change_points =  change_points, fdr =  fdr)
    R_MOPS_per <- Performance.fun(R_MOPS_sel, ge = cps, change_points = change_points, h = h)
    RMO_fdr[i] <- R_MOPS_per$fdr
    RMO_power[i] <- R_MOPS_per$power
    
    
  }
  print(c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), a))
  RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), a)
}

colnames(RESULT) <- c('BO_fdr', 'BO_power', 'MO_fdr', 'MO_power', 'RMO_fdr', 'RMO_power', 'a')
RESULT
write.csv(RESULT, file = './results20201109/ts_vary_a.csv', row.names = FALSE)


###### normal distribution with different correlation ########

set.seed(1234)
n = 4000
p = 50
cps = seq.int(200, 3800, by = 200)
a = 2
sparsity = 3
dis = 'normal'
dense = FALSE
fdr = 0.1
h = NULL
cutoff = 10

rho_c = seq(0.1, 0.9, by = 0.2)
#a_c <- seq(from = 3, to = 4, by = 0.2)


change_points <- seq.int(from = 0, to = n, by = 150)
change_points <- change_points[-1]
#change_points <- candidate_set_fun(cps = cps, h = 30, num_fds = 5)

#change_points <- cps

simulation_time <- 20

iter_time <- length(rho_c)

RESULT <- matrix(ncol = 7, nrow = iter_time)

for(iter in 1:iter_time){
  
  rho <- rho_c[iter]
  
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  RMO_fdr <- NULL
  RMO_power <- NULL
  
  for(i in 1:simulation_time){
    
    X <- MCP_model(n = n, p = p, cps = cps, dis = dis, a = a, sparsity = sparsity, rho = rho, dense = dense) 
    
    #Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = 0.2, cutoff = 5, h = NULL)
    Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 'inf', side_info = TRUE)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = cps, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    MOPS_sel <- MOPS(X = X, change_points = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = cps, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    R_MOPS_sel <- R_MOPS(X = X, change_points =  change_points, fdr =  fdr)
    R_MOPS_per <- Performance.fun(R_MOPS_sel, ge = cps, change_points = change_points, h = h)
    RMO_fdr[i] <- R_MOPS_per$fdr
    RMO_power[i] <- R_MOPS_per$power
  }
  print(c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), rho))
  RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), rho)
}

colnames(RESULT) <- c('BO_fdr', 'BO_power', 'MO_fdr', 'MO_power', 'RMO_fdr', 'RMO_power', 'rho')

RESULT
#write.csv(RESULT, file = './results20201109/ts_normal_vary_rho.csv', row.names = FALSE)

######  t distribution with different df ####

#n = 1000
#p = 10
#cps = seq.int(100, 900, by = 100)
#a = 1
#sparsity = p
#dis = 'chi'
#rho = 0
#dense = FALSE


n = 4000
p = 50
cps = seq.int(400, 3600, by = 400)
a = 1.5
sparsity = 3
dis = 'df'
dense = FALSE
fdr = 0.2
h = NULL
cutoff = 10

#rho_c = seq(0, 0.5, by = 0.1)
df_c <- seq(from = 6, to = 10, by = 1)

change_points <- seq.int(from = 0, to = n, by = 150)
change_points <- change_points[-1]
#change_points <- cps

simulation_time <- 20

iter_time <- length(df_c)

RESULT <- matrix(ncol = 7, nrow = iter_time)

for(iter in 1:iter_time){
  
  df <- df_c[iter]
  
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  RMO_fdr <- NULL
  RMO_power <- NULL
  
  for(i in 1:simulation_time){
    
    X <- MCP_model(n = n, p = p, cps = cps, dis = dis, a = a, df = df, sparsity = sparsity, rho = rho, dense = dense) 
    
    #Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = 0.2, cutoff = 5, h = NULL)
    Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 'inf', side_info = FALSE)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = cps, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    MOPS_sel <- MOPS(X = X, change_points = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = cps, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    R_MOPS_sel <- R_MOPS(X = X, change_points =  change_points, fdr =  fdr)
    R_MOPS_per <- Performance.fun(R_MOPS_sel, ge = cps, change_points = change_points, h = h)
    RMO_fdr[i] <- R_MOPS_per$fdr
    RMO_power[i] <- R_MOPS_per$power
  }
  RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), df)
}


RESULT

colnames(RESULT) <- c('BO_fdr', 'BO_power', 'MO_fdr', 'MO_power', 'RMO_fdr', 'RMO_power', 'df')

write.csv(RESULT, file = './results20201109/t_vary_df.csv', row.names = FALSE)





############### ############### ############### ############### ############### 
############### K-fwer   ###############################################
############### ############### ############### ############### ############### 
simulation_time = 25
fwer = 3
RESULT <- matrix(ncol = 7, nrow = iter_time)

for(iter in 1:iter_time){
  
  rho <- rho_c[iter]
  
  
  Bo_fwer <- NULL
  Bo_power <- NULL
  MO_fwer <- NULL
  MO_power <- NULL
  RMO_fwer <- NULL
  RMO_power <- NULL
  
  for(i in 1:simulation_time){
    
    X <- MCP_model(n = n, p = p, cps = cps, dis = dis, a = a, sparsity = sparsity, rho = rho, dense = dense) 
    
    #Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = 0.2, cutoff = 5, h = NULL)
    Bootstrap_sel <- Bootstrap_fwer(X = X, change_points = change_points, fwer = fwer, cutoff = cutoff, h = h, alpha = 'inf', side_info = FALSE)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = cps, change_points = change_points, h = h)
    Bo_fwer[i] <- Bootstrap_per$fwer
    Bo_power[i] <- Bootstrap_per$power
    
    MOPS_sel <- MOPS_fwer(X = X, change_points = change_points, fwer = fwer)
    MOPS_per <- Performance.fun(MOPS_sel, ge = cps, change_points = change_points, h = h)
    MO_fwer[i] <- MOPS_per$fwer
    MO_power[i] <- MOPS_per$power
    
    R_MOPS_sel <- R_MOPS_fwer(X = X, change_points =  change_points, fwer =  fwer)
    R_MOPS_per <- Performance.fun(R_MOPS_sel, ge = cps, change_points = change_points, h = h)
    RMO_fwer[i] <- R_MOPS_per$fwer
    RMO_power[i] <- R_MOPS_per$power
  }
  RESULT[iter,] <- c(mean(Bo_fwer), mean(Bo_power), mean(MO_fwer), mean(MO_power), mean(RMO_fwer), mean(RMO_power), rho)
}

RESULT



#################################################
#################### REGRESSION ##########################
#################################################

##### vary A
n = 8000
p = 10
cps = seq.int(1000, 7000, by = 1000)
#rho = 0
rho = 0.2
s = 10
dis = 'normal'
dense = FALSE
fdr = 0.2
h = NULL
cutoff = 30

#rho_c = seq(0.0, 0.5, by = 0.1)
#a_c <- seq(from = 0.15, to = 0.25, by = 0.02)
a_c <- seq(from = 0.25, to = 0.35, by = 0.02)

change_points <- seq.int(from = 0, to = n, by = 450)
change_points <- change_points[-c(1, length(change_points))]
#change_points <- seq.int(from = 1000, to = 7000, by = 1000)
#change_points <- seq.int(from = 0, to = n, by = 300)
#change_points <- change_points[-1]

#Performance.fun(sel = change_points, ge = cps, h = NULL, change_points = change_points)
#change_points <- cps 

#set.seed(1)
#change_points <- candidate_set_fun(cps = cps, h = 200, num_fds = 5)
simulation_time <- 20

iter_time <- length(a_c)

RESULT <- matrix(ncol = 7, nrow = iter_time)
#RESULT <- matrix(ncol = 9, nrow = iter_time)

res_time <- Sys.time()
for(iter in 1:iter_time){
  
  a <- a_c[iter]
  
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  DE_fdr <- NULL
  DE_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  RMO_fdr <- NULL
  RMO_power <- NULL
  
  for(i in 1:simulation_time){
    Tuple <- Regression_model(n = n, p = p, cps = cps, a = a, s = s, dis = dis, rho = rho, dense = dense)
    XX <-Tuple$X
    yy <- Tuple$y
    
    X <- t(XX) %*% diag(as.numeric(yy))
    #Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = 0.2, cutoff = 5, h = NULL)
    Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 2, side_info = TRUE)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = cps, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    
    #Derandom_sel <- DeRandom_Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 2, B = 500)
    #Derandom_per <- Performance.fun(Derandom_sel, ge = cps, change_points = change_points, h = h)
    #DE_fdr[i] <- Derandom_per$fdr
    #DE_power[i] <- Derandom_per$power
    
    MOPS_sel <- MOPS(X = X, change_points = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = cps, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    R_MOPS_sel <- R_MOPS(X = X, change_points =  change_points, fdr =  fdr)
    R_MOPS_per <- Performance.fun(R_MOPS_sel, ge = cps, change_points = change_points, h = h)
    RMO_fdr[i] <- R_MOPS_per$fdr
    RMO_power[i] <- R_MOPS_per$power
  }
  RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), a)
  #RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(DE_fdr), mean(DE_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), a)
}
Sys.time() - res_time

RESULT
colnames(RESULT) <- c('BO_fdr', 'BO_power', 'MO_fdr', 'MO_power', 'RMO_fdr', 'RMO_power', 'a')

write.csv(RESULT, file = './results20201109//regression_normalvaryA.csv', row.names = F)
####为什么H0没有信号的没有对称的感觉，怎么就全负了？
####### K-FWER 



##### t-distribuion
n = 8000
p = 10
cps = seq.int(1000, 7000, by = 1000)
rho = 0
#rho = 0.1
s = 5
dis = 't'
dense = FALSE
fdr = 0.2
h = NULL
cutoff = 30

#rho_c = seq(0.0, 0.5, by = 0.1)
a_c <- seq(from = 0.25, to = 0.35, by = 0.02)


change_points <- seq.int(from = 0, to = n, by = 450)
change_points <- change_points[-c(1, length(change_points))]
#change_points <- seq.int(from = 1000, to = 7000, by = 1000)
#change_points <- seq.int(from = 0, to = n, by = 300)
#change_points <- change_points[-1]

#Performance.fun(sel = change_points, ge = cps, h = NULL, change_points = change_points)
#change_points <- cps 

#set.seed(1)
#change_points <- candidate_set_fun(cps = cps, h = 200, num_fds = 5)
simulation_time <- 20

iter_time <- length(a_c)

RESULT <- matrix(ncol = 7, nrow = iter_time)
#RESULT <- matrix(ncol = 9, nrow = iter_time)

res_time <- Sys.time()
for(iter in 1:iter_time){
  
  a <- a_c[iter]
  
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  #DE_fdr <- NULL
  #DE_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  RMO_fdr <- NULL
  RMO_power <- NULL
  
  for(i in 1:simulation_time){
    Tuple <- Regression_model(n = n, p = p, cps = cps, a = a, s = s, dis = dis, rho = rho, dense = dense)
    XX <-Tuple$X
    yy <- Tuple$y
    
    X <- t(XX) %*% diag(as.numeric(yy))
    #Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = 0.2, cutoff = 5, h = NULL)
    Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 2, side_info = TRUE)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = cps, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    
    #Derandom_sel <- DeRandom_Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 2, B = 500)
    #Derandom_per <- Performance.fun(Derandom_sel, ge = cps, change_points = change_points, h = h)
    #DE_fdr[i] <- Derandom_per$fdr
    #DE_power[i] <- Derandom_per$power
    
    MOPS_sel <- MOPS(X = X, change_points = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = cps, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    R_MOPS_sel <- R_MOPS(X = X, change_points =  change_points, fdr =  fdr)
    R_MOPS_per <- Performance.fun(R_MOPS_sel, ge = cps, change_points = change_points, h = h)
    RMO_fdr[i] <- R_MOPS_per$fdr
    RMO_power[i] <- R_MOPS_per$power
  }
  RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), a)
  #RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(DE_fdr), mean(DE_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), a)
}
Sys.time() - res_time
colnames(RESULT) <- c('BO_fdr', 'BO_power', 'MO_fdr', 'MO_power', 'RMO_fdr', 'RMO_power', 'a')
RESULT

write.csv(RESULT, file = './results20201109//regression_varyA_disT.csv', row.names = F)


############## rho #######
n = 8000
p = 10
cps = seq.int(1000, 7000, by = 1000)
a = 0.25
s = 5
dis = 'normal'
dense = FALSE
fdr = 0.2
h = NULL
cutoff = 30

#rho_c = seq(0.0, 0.5, by = 0.1)
rho_c <- seq(from = 0, to = 0.5, by = 0.1)

change_points <- seq.int(from = 0, to = n, by = 450)
change_points <- change_points[-c(1, length(change_points))]
#change_points <- seq.int(from = 1000, to = 7000, by = 1000)
#change_points <- seq.int(from = 0, to = n, by = 300)
#change_points <- change_points[-1]

#Performance.fun(sel = change_points, ge = cps, h = NULL, change_points = change_points)
#change_points <- cps 

#set.seed(1)
#change_points <- candidate_set_fun(cps = cps, h = 200, num_fds = 5)
simulation_time <- 20

iter_time <- length(rho_c)

RESULT <- matrix(ncol = 7, nrow = iter_time)
#RESULT <- matrix(ncol = 9, nrow = iter_time)

res_time <- Sys.time()
for(iter in 1:iter_time){
  
  rho <- rho_c[iter]
  
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  #DE_fdr <- NULL
  #DE_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  RMO_fdr <- NULL
  RMO_power <- NULL
  
  for(i in 1:simulation_time){
    Tuple <- Regression_model(n = n, p = p, cps = cps, a = a, s = s, dis = dis, rho = rho, dense = dense)
    XX <-Tuple$X
    yy <- Tuple$y
    
    X <- t(XX) %*% diag(as.numeric(yy))
    #Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = 0.2, cutoff = 5, h = NULL)
    Bootstrap_sel <- Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 2, side_info = TRUE)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = cps, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    
    #Derandom_sel <- DeRandom_Bootstrap_fdr(X = X, change_points = change_points, fdr = fdr, cutoff = cutoff, h = h, alpha = 2, B = 500)
    #Derandom_per <- Performance.fun(Derandom_sel, ge = cps, change_points = change_points, h = h)
    #DE_fdr[i] <- Derandom_per$fdr
    #DE_power[i] <- Derandom_per$power
    
    MOPS_sel <- MOPS(X = X, change_points = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = cps, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    R_MOPS_sel <- R_MOPS(X = X, change_points =  change_points, fdr =  fdr)
    R_MOPS_per <- Performance.fun(R_MOPS_sel, ge = cps, change_points = change_points, h = h)
    RMO_fdr[i] <- R_MOPS_per$fdr
    RMO_power[i] <- R_MOPS_per$power
  }
  RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), rho)
  #RESULT[iter,] <- c(mean(Bo_fdr), mean(Bo_power), mean(DE_fdr), mean(DE_power), mean(MO_fdr), mean(MO_power), mean(RMO_fdr), mean(RMO_power), a)
}
Sys.time() - res_time
colnames(RESULT) <- c('BO_fdr', 'BO_power', 'MO_fdr', 'MO_power', 'RMO_fdr', 'RMO_power', 'a')
RESULT

write.csv(RESULT, file = './results20201109//regression_varyrho.csv', row.names = F)


