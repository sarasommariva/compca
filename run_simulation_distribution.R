rm(list=ls())
setwd(getSrcDirectory(function(){})[1])

source("utils.R")
library(compositions)

# Step 1. Parameters definition
K_sim = 2 # Value of K used to simulate the data (H0: K=2)
eigenval_y=c(10, 9, 1, 1, 0.5)
eigenval_z=c(6, 5, 1, 0.9, 0.3, 0.1, 0.02)
n_y = 100; n_z = 100

D = 8
Q = 2
K_H0 = 2
num_sim = 1000

# - Simplex basis for the ilr
V_ilr_y = def_base_ilr(D-Q)
V_ilr_z = def_base_ilr(D)

machine_toll = 5*.Machine$double.eps

folder_res = file.path('.', 'results')
dir.create(folder_res, showWarnings = FALSE)

# Step 2. Randomly generate covariance matrix
temp_covs = rand_covmat(D, K_sim, Q, eigenval_y, eigenval_z)
omega_y = temp_covs$omega_y; omega_z = temp_covs$omega_z; rm(temp_covs)

# Step 3. Compute Schott's parameters by using the true covariance matrices
fake_stat = NA
schott = compute_pvalue_schott_theo(K_H0, fake_stat, omega_y, omega_z, n_y, n_z)

# Step 4. Simulate data and compute test statistics
distributions = c('Gaussian', 'Uniform', 'Mult4', 'Mult8', 'Mult40')

stat_values_gaussian = matrix(data=NA, nrow=num_sim)
stat_values_uniform = matrix(data=NA, nrow=num_sim)
stat_values_mult4 = matrix(data=NA, nrow=num_sim)
stat_values_mult8 = matrix(data=NA, nrow=num_sim)
stat_values_mult40 = matrix(data=NA, nrow=num_sim)

for (idist in 1:length(distributions)){
  
  distribution = distributions[idist]
  if (distribution == 'Gaussian'){df=NULL; aux_distr='Gaussian'
  }else if(distribution == 'Mult4'){df = 4; aux_distr='Mult'
  }else if(distribution == 'Mult8'){df = 8; aux_distr='Mult'
  }else if(distribution == 'Mult40'){df = 40; aux_distr='Mult'
  }else if(distribution == 'Uniform'){df = NULL; aux_distr='Uniform'}
  
  for (isim in 1:num_sim){
    
    if (isim %% 10 == 1){print(paste0(distribution, ' simulation num = ', isim))}
    
    # 4.a. Generate in R^(D-1)
    temp_dataset = rand_sample(n_y, n_z, omega_y, omega_z, aux_distr, df)
    Y_or = temp_dataset$Y_or; Z_or = temp_dataset$Z_or; rm(temp_dataset)
    # 4.b. Derive the compositional dataset in S_D 
    Y_comp = matrix(data=NA, nrow=n_y, ncol=D-Q)
    Z_comp = matrix(data=NA, nrow=n_z, ncol=D)
    for (id in 1:n_y){
      Y_comp[id,] <- ilrInv(Y_or[id,], V_ilr_y)
    }
    for (id in 1:n_z){
      Z_comp[id,] <- ilrInv(Z_or[id,], V_ilr_z)
    }
    if (!all(abs(rowSums(Y_comp)-1)<machine_toll) || 
        !all(abs(rowSums(Z_comp)-1)<machine_toll)){
      stop("Data are note 1-closed compositions")
    }
    # 4.c. Compute ilr transformed data
    Y_transf = matrix(data=NA, nrow=n_y, ncol=D-Q-1)
    Z_transf = matrix(data=NA, nrow=n_z, ncol=D-1)
    for (id in 1:n_y){
      Y_transf[id,] <- ilr(Y_comp[id,], V=V_ilr_y)
    }
    for (id in 1:n_z){
      Z_transf[id,] <- ilr(Z_comp[id,], V=V_ilr_z)
    }
    aux = compute_stat_value(K_H0, Y_transf, Z_transf)
  
    if (distribution=='Gaussian'){
      stat_values_gaussian[isim] = aux$stat_value
    }else if(distribution=='Mult4'){
      stat_values_mult4[isim]= aux$stat_value
    }else if(distribution=='Mult8'){
      stat_values_mult8[isim]= aux$stat_value
    }else if(distribution=='Mult40'){
      stat_values_mult40[isim]= aux$stat_value
    }else if(distribution=='Uniform'){
      stat_values_uniform[isim]= aux$stat_value
    }
  }
  rm(df)
}

# Step 5. Save
result = list('param' = list('D'=D, 'Q'=Q, 'K_H0'=K_H0, 'K_sim'=K_sim,
                             "n_y" = n_y, "n_z" = n_z,
                             'eigenval_y'=eigenval_y, 'eigenval_z'=eigenval_z), 
              'schott.mu_T'=schott$mu_T, 'schott.sigma2_T'=schott$sigma2_T,
              'stat_values_gaussian'=stat_values_gaussian, 
              'stat_values_mult4'=stat_values_mult4,
              'stat_values_mult8'=stat_values_mult8,
              'stat_values_mult40'=stat_values_mult40,
              'stat_values_uniform'=stat_values_uniform)
save(result, file=file.path(folder_res, 
                             paste0('empirical_distr_ny', n_y, '_nz', n_z,'.Rdata')))
 
# # Plot
# const_chi = 0.5*result$schott.sigma2_T / result$schott.mu_T
# df_chi = max(round(2*result$schott.mu_T^2/result$schott.sigma2_T ), 1)
# delta_bin = 0.1
# 
# points = seq(min(result$stat_values_gaussian)-delta_bin, 
#              max(result$stat_values_gaussian)+delta_bin, by=delta_bin)
# dist_chi = const_chi^-1 * dchisq(points/const_chi, df_chi)
# h = hist(result$stat_values_gaussian, breaks=300, freq=FALSE)
# lines(points, dist_chi)
# 
# points = seq(min(result$stat_values_gaussian)-delta_bin, 
#              max(result$stat_values_gaussian)+delta_bin, by=delta_bin)
# dist_chi = const_chi^-1 * dchisq(points/const_chi, df_chi)
# h = hist(result$stat_values_mult100, breaks=300, freq=FALSE)
# lines(points, dist_chi)
# 
# points = seq(min(result$stat_values_gaussian)-delta_bin, 
#              max(result$stat_values_gaussian)+delta_bin, by=delta_bin)
# dist_chi = const_chi^-1 * dchisq(points/const_chi, df_chi)
# h = hist(result$stat_values_mult4, breaks=300, freq=FALSE)
# lines(points, dist_chi)

