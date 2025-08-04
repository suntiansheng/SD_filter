setwd("~/bootstrap_fdr")
source('./utility_bootstrap.R')
library(foreach)
library(doParallel)

######  t distribution with different df ####
n = 2000
p = 10
bw = 200
fdr = 0.2
sigma = 1
h = NULL
#a_c <- seq(from = 3, to = 4, by = 0.2)
df_c <- seq.int(from = 5, to = 10, by = 1)
a = 3
sparsity = 3






cl <- parallel::makeCluster(length(df_c))
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/suntiansheng/R/x86_64-pc-linux-gnu-library/3.5/"))

### initial
iter_time <- length(df_c)
simulation_time = 20
change_points <- seq.int(from = 0, to = n, by = 75)
change_points <- change_points[-1]

res_time <- Sys.time()
result_matrix <- foreach(iter = 1:iter_time, .combine = 'rbind') %dopar% {
  df <- df_c[iter]
  
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  Sa_fdr <- NULL
  Sa_power <- NULL
  Or_fdr <- NULL
  Or_power <- NULL
  multi_fdr <- NULL
  multi_power <- NULL
  
  for(i in 1:simulation_time){
    
    DATA_tuple <- Data_t(n, p, bw, a, sparsity, df)
    X <- DATA_tuple$x
    gen <- DATA_tuple$true_position
    
    
    Bootstrap_sel <- Bootstrap_fdr(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = gen, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    MOPS_sel <- MOPS(X = X, changepos = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = gen, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    SaRa_sel <- SaRa(X = X, change_points = change_points, fdr = fdr, sigma = sigma, h = h)
    SaRa_per <- Performance.fun(SaRa_sel, ge = gen, change_points = change_points, h = h)
    Sa_fdr[i] <- SaRa_per$fdr
    Sa_power[i] <- SaRa_per$power
    
    multi_sel <- multiple_bootstrap(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    multi_per <- Performance.fun(multi_sel, ge = gen, change_points = change_points, h = h)
    multi_fdr[i] <- multi_per$fdr
    multi_power[i] <- multi_per$power
    
  }
  
  c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(Sa_fdr), mean(Sa_power), mean(multi_fdr), mean(multi_power), df)
  
}
Sys.time() - res_time


colnames(result_matrix) <- c('Bo_fdr', 'Bo_power', 'Mo_fdr', 'Mo_power', 'Sa_fdr', 'Sa_power', 'multi_fdr','multi_power','df')
result_matrix

#write.csv(result_matrix, file = './result/vary_df.csv', row.names = FALSE)




######  normal distribution with different sigma ####
n = 2000
p = 10
bw = 200
fdr = 0.2
sigma = 1
h = NULL
#a_c <- seq(from = 3, to = 4, by = 0.2)
sd_c <- seq.int(from = 1, to = 2, by = 0.2)
a = 3
sparsity = 3

### initial change position
change_points <- seq.int(from = 0, to = n, by = 75)
change_points <- change_points[-1]
#change_points <- c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800)


cl <- parallel::makeCluster(length(sd_c))
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/suntiansheng/R/x86_64-pc-linux-gnu-library/3.5/"))

### initial
iter_time <- length(sd_c)
simulation_time = 20
change_points <- seq.int(from = 0, to = n, by = 75)
change_points <- change_points[-1]

res_time <- Sys.time()
result_matrix <- foreach(iter = 1:iter_time, .combine = 'rbind') %dopar% {
  sd <- sd_c[iter]  
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  Sa_fdr <- NULL
  Sa_power <- NULL
  Or_fdr <- NULL
  Or_power <- NULL
  multi_fdr <- NULL
  multi_power <- NULL
  
  for(i in 1:simulation_time){
    
    DATA_tuple <- Data_norm(n, p, bw, a, sparsity, sd)
    X <- DATA_tuple$x
    gen <- DATA_tuple$true_position
    
    
    
    Bootstrap_sel <- Bootstrap_fdr(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = gen, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    MOPS_sel <- MOPS(X = X, changepos = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = gen, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    SaRa_sel <- SaRa(X = X, change_points = change_points, fdr = fdr, sigma = sigma, h = h)
    SaRa_per <- Performance.fun(SaRa_sel, ge = gen, change_points = change_points, h = h)
    Sa_fdr[i] <- SaRa_per$fdr
    Sa_power[i] <- SaRa_per$power
    
    multi_sel <- multiple_bootstrap(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    multi_per <- Performance.fun(multi_sel, ge = gen, change_points = change_points, h = h)
    multi_fdr[i] <- multi_per$fdr
    multi_power[i] <- multi_per$power
    
  }
  
  c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(Sa_fdr), mean(Sa_power), mean(multi_fdr), mean(multi_power), sd)
  
}
Sys.time() - res_time


colnames(result_matrix) <- c('Bo_fdr', 'Bo_power', 'Mo_fdr', 'Mo_power', 'Sa_fdr', 'Sa_power', 'multi_fdr','multi_power','df')
result_matrix
#write.csv(result_matrix, file = './result/vary_sd.csv', row.names = FALSE)



########### vary a in t distribution ################
n = 2000
p = 10
bw = 200
fdr = 0.2
sigma = 1
h = NULL
df = 7
a_c <- seq(from = 2, to = 6, by = 1)
#df_c <- seq.int(from = 5, to = 10, by = 1)
sparsity = 3

### initial
change_points <- seq.int(from = 0, to = n, by = 75)
change_points <- change_points[-1]
#change_points <- c(100, 200, 300, 400, 500, 600, 700, 800, 900)

cl <- parallel::makeCluster(length(a_c))
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/suntiansheng/R/x86_64-pc-linux-gnu-library/3.5/"))

### initial
iter_time <- length(a_c)
simulation_time = 20
change_points <- seq.int(from = 0, to = n, by = 75)
change_points <- change_points[-1]

res_time <- Sys.time()
result_matrix <- foreach(iter = 1:iter_time, .combine = 'rbind') %dopar% {
  a <- a_c[iter]
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  Sa_fdr <- NULL
  Sa_power <- NULL
  Or_fdr <- NULL
  Or_power <- NULL
  multi_fdr <- NULL
  multi_power <- NULL
  
  for(i in 1:simulation_time){
    
    DATA_tuple <- Data_t(n, p, bw, a, sparsity, df)
    X <- DATA_tuple$x
    gen <- DATA_tuple$true_position
    
    
    
    Bootstrap_sel <- Bootstrap_fdr(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = gen, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    MOPS_sel <- MOPS(X = X, changepos = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = gen, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    SaRa_sel <- SaRa(X = X, change_points = change_points, fdr = fdr, sigma = sigma, h = h)
    SaRa_per <- Performance.fun(SaRa_sel, ge = gen, change_points = change_points, h = h)
    Sa_fdr[i] <- SaRa_per$fdr
    Sa_power[i] <- SaRa_per$power
    
    multi_sel <- multiple_bootstrap(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    multi_per <- Performance.fun(multi_sel, ge = gen, change_points = change_points, h = h)
    multi_fdr[i] <- multi_per$fdr
    multi_power[i] <- multi_per$power
    
  }
  
  c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(Sa_fdr), mean(Sa_power), mean(multi_fdr), mean(multi_power), a)
  
}
Sys.time() - res_time
colnames(result_matrix) <- c('Bo_fdr', 'Bo_power', 'Mo_fdr', 'Mo_power', 'Sa_fdr', 'Sa_power', 'multi_fdr','multi_power','a')
result_matrix
#write.csv(result_matrix, file = './result/vary_a.csv', row.names = FALSE)







########### vary change number ################
n = 2000
p = 10
bw = 200
fdr = 0.2
sigma = 1
h = NULL
df = 7
a = 3
by_c <- seq(from = 75, to = 200, by = 20)
#df_c <- seq.int(from = 5, to = 10, by = 1)
sparsity = 3


cl <- parallel::makeCluster(length(by_c))
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/suntiansheng/R/x86_64-pc-linux-gnu-library/3.5/"))


### initial
iter_time <-length(by_c)
simulation_time = 20
#change_points <- seq.int(from = 0, to = n, by = 75)
#change_points <- change_points[-1]

res_time <- Sys.time()
result_matrix <- foreach(iter = 1:iter_time, .combine = 'rbind') %dopar% {
  by <- by_c[iter]
  change_points <- seq.int(from = 0, to = n, by = by)
  change_points <- change_points[-1]
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  Sa_fdr <- NULL
  Sa_power <- NULL
  multi_fdr <- NULL
  multi_power <- NULL
  
  for(i in 1:simulation_time){
    
    DATA_tuple <- Data_t(n, p, bw, a, sparsity, df)
    X <- DATA_tuple$x
    gen <- DATA_tuple$true_position
    
    
    
    Bootstrap_sel <- Bootstrap_fdr(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = gen, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    MOPS_sel <- MOPS(X = X, changepos = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = gen, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    SaRa_sel <- SaRa(X = X, change_points = change_points, fdr = fdr, sigma = sigma, h = h)
    SaRa_per <- Performance.fun(SaRa_sel, ge = gen, change_points = change_points, h = h)
    Sa_fdr[i] <- SaRa_per$fdr
    Sa_power[i] <- SaRa_per$power
    
    multi_sel <- multiple_bootstrap(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    multi_per <- Performance.fun(multi_sel, ge = gen, change_points = change_points, h = h)
    multi_fdr[i] <- multi_per$fdr
    multi_power[i] <- multi_per$power
    
  }
  
  c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(Sa_fdr), mean(Sa_power), mean(multi_fdr), mean(multi_power), by)
  
}
Sys.time() - res_time


colnames(result_matrix) <- c('Bo_fdr', 'Bo_power', 'Mo_fdr', 'Mo_power', 'Sa_fdr', 'Sa_power', 'multi_fdr','multi_power','by')
result_matrix[,9] <- by_c
#write.csv(result_matrix, file = './result/vary_dis.csv', row.names = FALSE)




######  MVA t distribution with different covariance matrix ####
n = 2000
p = 10
bw = 200
fdr = 0.2
sigma = 1
h = NULL
#a_c <- seq(from = 3, to = 4, by = 0.2)
phi_c <- seq(from = 0, to = 0.5, length.out = 6)
df = 7
a = 3
sparsity = 3

### initial
change_points <- seq.int(from = 0, to = n, by = 75)
change_points <- change_points[-1]
#change_points <- c(100, 200, 300, 400, 500, 600, 700, 800, 900)




cl <- parallel::makeCluster(length(phi_c))
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/suntiansheng/R/x86_64-pc-linux-gnu-library/3.5/"))


### initial
iter_time <-length(phi_c)
simulation_time = 20
#change_points <- seq.int(from = 0, to = n, by = 75)
#change_points <- change_points[-1]

res_time <- Sys.time()
result_matrix <- foreach(iter = 1:iter_time, .combine = 'rbind') %dopar% {
  phi <- phi_c[iter]
  
  Bo_fdr <- NULL
  Bo_power <- NULL
  MO_fdr <- NULL
  MO_power <- NULL
  Sa_fdr <- NULL
  Sa_power <- NULL
  multi_fdr <- NULL
  multi_power <- NULL
  
  for(i in 1:simulation_time){
    
    
    DATA_tuple <- Data_mvt(n, p, bw, a, sparsity, df, phi = phi)
    X <- DATA_tuple$x
    gen <- DATA_tuple$true_position
    
    
    
    Bootstrap_sel <- Bootstrap_fdr(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    Bootstrap_per <- Performance.fun(Bootstrap_sel, ge = gen, change_points = change_points, h = h)
    Bo_fdr[i] <- Bootstrap_per$fdr
    Bo_power[i] <- Bootstrap_per$power
    
    MOPS_sel <- MOPS(X = X, changepos = change_points, fdr = fdr)
    MOPS_per <- Performance.fun(MOPS_sel, ge = gen, change_points = change_points, h = h)
    MO_fdr[i] <- MOPS_per$fdr
    MO_power[i] <- MOPS_per$power
    
    SaRa_sel <- SaRa(X = X, change_points = change_points, fdr = fdr, sigma = sigma, h = h)
    SaRa_per <- Performance.fun(SaRa_sel, ge = gen, change_points = change_points, h = h)
    Sa_fdr[i] <- SaRa_per$fdr
    Sa_power[i] <- SaRa_per$power
    
    multi_sel <- multiple_bootstrap(X, change_points = change_points, fdr = fdr, cutoff = 5, h = h)
    multi_per <- Performance.fun(multi_sel, ge = gen, change_points = change_points, h = h)
    multi_fdr[i] <- multi_per$fdr
    multi_power[i] <- multi_per$power
    
  }
  
  c(mean(Bo_fdr), mean(Bo_power), mean(MO_fdr), mean(MO_power), mean(Sa_fdr), mean(Sa_power), mean(multi_fdr), mean(multi_power), phi)
  
}
Sys.time() - res_time


#write.csv(result_matrix, file = './result/vary_phi.csv', row.names = FALSE)

