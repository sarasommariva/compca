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
n_y = 100
n_z = n_y
K = 2
Q = 2

# Bootstrap parameters
num_rip_boot = 100

mu_y = rep(0, D-Q-1)
mu_z= rep(0, D-1)
#mu_y = c(20, 10, 3)
#mu_z = c(2, 1, 0, 10, 10)

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

est_mu_T = matrix(data=NA, nrow=num_sim, ncol=1)
est_sigma2_T = matrix(data=NA, nrow=num_sim, ncol=1)
est_pvalues = matrix(data=NA, nrow=num_sim, ncol=1)

boot_mu_T = matrix(data=NA, nrow=num_sim, ncol=1)
boot_sigma2_T = matrix(data=NA, nrow=num_sim, ncol=1)
boot_pvalues = matrix(data=NA, nrow=num_sim, ncol=1)


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
  
  # Step 7. Estimate the p-value
  #   7.a. Estimate eigenvalues
  est_alpha = c(eigen_cov_y$values/(n_y-1), matrix(0, Q, 1))
  est_beta = eigen_cov_z$values/(n_z -1)
  est_psi =  eigen_sum$values / (n_y + n_z - 2)
  #   7.b. Estimate eigenvectors
  est_mat_K = eigen_sum$vectors
  temp1 = eigen(t(est_mat_K[,1:K]) %*% 
                  adiag(est_cov_y, matrix(0, Q, Q)) %*% est_mat_K[,1:K])
  temp2 = eigen(t(est_mat_K[,(K+1):(D-1)]) %*% 
                  adiag(est_cov_y, matrix(0, Q, Q)) %*% est_mat_K[,(K+1):(D-1)])
  est_mat_U = adiag(temp1$vectors, temp2$vectors); rm(temp1, temp2)
  temp1 = eigen(t(est_mat_K[,1:K]) %*% est_cov_z %*% est_mat_K[,1:K])
  temp2 = eigen(t(est_mat_K[,(K+1):(D-1)]) %*% est_cov_z %*% est_mat_K[,(K+1):(D-1)])
  est_mat_V = adiag(temp1$vectors, temp2$vectors); rm(temp1, temp2)
  #   7.c. Estimate parameters of the null distribution of the test statistic
  est_paramsT = compute_param_nulldistr(est_alpha, est_beta, est_psi, 
                                        est_mat_U, est_mat_V, n_y, n_z)
  est_mu_T[isim] = est_paramsT$mu_T
  est_sigma2_T[isim] = est_paramsT$sigma2_T
  #   7.d. Compute p-value
  est_const_chi = 0.5*est_sigma2_T[isim] / est_mu_T[isim]
  est_df_chi = max(round(2*est_mu_T[isim]^2/est_sigma2_T[isim]), 1)
  
  est_pvalues[isim] = pchisq(stat_values[isim]/est_const_chi,
                              est_df_chi, lower.tail = FALSE)
  
  # Step 8. Estimate p-value via bootstrap
  boot_stat_values = matrix(data=NA, nrow=num_rip_boot, ncol=1)

  #   8.1. Randomly rotate the eigenvector of Y.
  Rb_1 = randortho(K, type="orthonormal")
  Rb_2 = randortho(D-K-1, type="orthonormal")
  est_mat_Ub = adiag(eigen_cov_y$vectors, diag(Q)) %*% adiag(Rb_1, Rb_2)
  #   8.2. Rotate the dataset Z so that its PCs coincide with those just definied
  Rb = est_mat_Ub %*% t(eigen_cov_z$vectors)
  Zb_transf = Z_transf %*% t(Rb)
  for (iboot in 1:num_rip_boot){
  #   8.3. Define a bootstrap dataset
       Y_boot = Y_transf[randi(n_y, n_y, 1), ]
       Z_boot = Zb_transf[randi(n_z, n_z, 1), ]
  #   8.4. Compute value of the test statistic for the bootstrap samplee
       boot_cov_y = crossprod(sweep(Y_boot, 2, colMeans(Y_boot)))
       boot_cov_z = crossprod(sweep(Z_boot, 2, colMeans(Z_boot)))
       boot_eigen_cov_y = eigen(boot_cov_y)
       boot_eigen_cov_z = eigen(boot_cov_z)
       boot_eigen_sum = eigen(adiag(boot_cov_y, matrix(0, Q, Q)) + boot_cov_z)
   
       boot_stat_values[iboot] = sum(boot_eigen_cov_y$values[1:K]
                   + boot_eigen_cov_z$values[1:K] - boot_eigen_sum$values[1:K])
  }
  #   8.5. Estimate the bootstrap p-value
  boot_pvalues[isim] = sum(boot_stat_values > stat_values[isim])/num_rip_boot
  
  #rm(boot_stat_values)
}

hist(stat_values, 200)

result = list('param' = list('n_y' = n_y, 'n_z' = n_z, 'K'=K, 'Q'=Q,
                          'omega_y' = omega_y, 'omega_z' = omega_z),
                    'stat_values'=stat_values, 'est_mu_T'=est_mu_T,
                    'est_sigma2_T'=est_sigma2_T, 'est_pvalues'=est_pvalues,
                    'num_rip_boot'=num_rip_boot, 'boot_pvalues'=boot_pvalues)
save(result, file="prova_h0_nynz100.Rdata")

#mu_stat_values = mean(stat_values)
#var_stat_values = var(stat_values)
#print(2*mu_stat_values/var_stat_values)

# Todo:
# 1. Considerare altre medie e matrici di covarianza.
#    Eventualmente creare un dataset dedicato per gestire la loro definizione.

# Commenti:
# 1. ilrinv da composizioni che sommano a 1 a meno della precisione di macchina