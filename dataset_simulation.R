rm(list=ls())
setwd('~/Documenti/compoteam/paper_stattest/compca')
source("utils.R")

#install.packages("rockchalk")
library(rockchalk)
library(compositions)
library(utils)

# Step 1. Definition of the parameters
machine_toll = 5*.Machine$double.eps

num_sim = 1000
D = 6
n_y = 500
n_z = n_y
K = 2
Q = 2

#mu_y = rep(0, D-Q-1)
#mu_z= rep(0, D-1)
mu_y = c(20, 10, 3)
mu_z = c(2, 1, 0, 10, 10)

# Case 1. Diagonal matrices with structural zeros
#omega_y = diag(c(10, 5, 3, 0, 0))
#omega_z = diag(c(8, 7, 3, 2, 1))

# Case 2. Non-diagonal matrices with structural zeros
temp_covs = rand_covmat(D, K, Q, eigenval_y=c(10, 5, 3), 
                        eigenval_z=c(8, 7, 3, 2, 1))

# Case 3. Matrices that do not satisfy H_0
#temp_covs = rand_indep_covmat(D, K, Q, eigenval_y=c(10, 5, 3), 
#                        eigenval_z=c(8, 7, 3, 2, 1))

omega_y = temp_covs$omega_y
omega_z = temp_covs$omega_z
rm(temp_covs)

stat_values = matrix(data=NA, nrow=num_sim, ncol=1)

for (isim in 1:num_sim){
  
  if (isim %% 10 == 1){print(paste0('Simulation num = ', isim))}

  # Step 2. Random samples from multivariate Gaussian in R^(D-1)
  Y_or = mvrnorm(n=n_y, mu=mu_y, Sigma=omega_y)
  Z_or = mvrnorm(n=n_z, mu=mu_z, Sigma=omega_z)
  
  # Step 3. Derive the compositional dataset in S_D 
  Y_comp = matrix(data=NA, nrow=n_y, ncol=D-Q)
  Z_comp = matrix(data=NA, nrow=n_z, ncol=D)
  for (id in 1:n_y){
    Y_comp[id,] <- ilrInv(Y_or[id,])
  }
  for (id in 1:n_z){
    Z_comp[id,] <- ilrInv(Z_or[id,])
  }
  
  if (!all(abs(rowSums(Y_comp)-1)<machine_toll) || 
      !all(abs(rowSums(Z_comp)-1)<machine_toll)){
    stop("Data are note 1-closed compositions")
  }
  
  # Step 4. Compute ilr transformed data
  Y_transf = matrix(data=NA, nrow=n_y, ncol=D-Q-1)
  Z_transf = matrix(data=NA, nrow=n_z, ncol=D-1)
  for (id in 1:n_y){
    Y_transf[id,] <- ilr(Y_comp[id,])
  }
  for (id in 1:n_z){
    Z_transf[id,] <- ilr(Z_comp[id,])
  }
  
  #print("Difference between original and transformed dataset")
  #print(paste0("Y -> ", max(abs(Y_transf-Y_or))))
  #print(paste0("Z -> ", max(abs(Z_transf-Z_or))))
  
  # !!!! Note: Probably Step 3 and 4 are not really necessary !!!!!
  
  # Step 5. Estimate sample covariance matrix
  est_cov_y = crossprod(sweep(Y_transf, 2, colMeans(Y_transf)))
  est_cov_z = crossprod(sweep(Z_transf, 2, colMeans(Z_transf)))
  
  # Step 6. Estimate sample value of the test statistic.
  eigen_cov_y = eigen(est_cov_y)
  eigen_cov_z = eigen(est_cov_z)
  eigen_sum = eigen(adiag(est_cov_y, matrix(0, Q, Q)) + est_cov_z)
  
  stat_values[isim] = sum(eigen_cov_y$values[1:K] + eigen_cov_z$values[1:K] - 
                   eigen_sum$values[1:K])
}

hist(stat_values, 200)

result = list('param' = list('n_y' = n_y, 'n_z' = n_z, 'K'=K, 'Q'=Q,
                          'omega_y' = omega_y, 'omega_z' = omega_z),
                    'stat_values'=stat_values)
save(result, file="prova_zeros_means.Rdata")

#mu_stat_values = mean(stat_values)
#var_stat_values = var(stat_values)
#print(2*mu_stat_values/var_stat_values)

# Todo:
# 1. Considerare altre medie e matrici di covarianza.
#    Eventualmente creare un dataset dedicato per gestire la loro definizione.

# Commenti:
# 1. ilrinv da composizioni che sommano a 1 a meno della precisione di macchina