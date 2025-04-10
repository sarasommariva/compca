rm(list=ls())
setwd(getSrcDirectory(function(){})[1])

source("utils.R")
library(rockchalk)
library(compositions)
library(utils)


# Step 1. Parameters definition

#   1.1 Parameters to be explored
K_sim = 1 # Value of K used to simulate the data (H0: K=2)
distribution = 'Uniform'

#   1.2. Pre-determined parameters
D = 8
Q = 2
K_H0 = 2  # Value of K used with the null hypothesis
ny_all= c(20, 60, 100)
nz_all = c(20, 60, 100)
num_sim = 1000
num_boot = 100

# - Simplex basis for the ilr
V_ilr_y = def_base_ilr(D-Q)
V_ilr_z = def_base_ilr(D)

machine_toll = 5*.Machine$double.eps

# Step 3. Initialization
folder_res = file.path('.', 'results')
dir.create(folder_res, showWarnings = FALSE)

results_list = list()
results_list_name = list()
for (iy in 1:length(ny_all)){
  for (iz in 1:length(nz_all)){
    results_list[[length(nz_all)*(iy-1)+iz]] = data.frame(
      stat_values = matrix(data=NA, nrow=num_sim),
      pvalues_schott_est = matrix(data=NA, nrow=num_sim),
      pvalues_schott_theo = matrix(data=NA, nrow=num_sim),
      pvalues_boot = matrix(data=NA, nrow=num_sim),

      err_ilr_Y = matrix(data=NA, nrow=num_sim),
      err_ilr_Z = matrix(data=NA, nrow=num_sim),
      est_mu_T = matrix(data=NA, nrow=num_sim),
      est_sigma2_T = matrix(data=NA, nrow=num_sim),
      theo_mu_T = matrix(data=NA, nrow=num_sim),
      theo_sigma2_T = matrix(data=NA, nrow=num_sim)
    )
    results_list_name[[length(nz_all)*(iy-1)+iz]] = paste0('res_ny', ny_all[iy], 
                                                           '_nz', nz_all[iz])
  }
}

for (isim in 1:num_sim){

  if (isim %% 10 == 1){print(paste0('Simulation num = ', isim))}
  
# Step 4. Randomly generate covariance matrix
  temp_covs = rand_covmat(D, K_sim, Q, eigenval_y=c(10, 5, 3, 1, 0.5), 
                          eigenval_z=c(8, 7, 3, 2, 0.3, 0.1, 0.02))
  omega_y = temp_covs$omega_y; omega_z = temp_covs$omega_z; rm(temp_covs)

  for (iy in 1:length(ny_all)){
    for (iz in 1:length(nz_all)){

      n_y = ny_all[iy]; n_z = nz_all[iz]; aux_idx = length(nz_all)*(iy-1)+iz
        
# Step 5. Randomly generate data
      # 5.a. Generate in R^(D-1)
      temp_dataset = rand_sample(n_y, n_z, omega_y, omega_z, distribution)
      Y_or = temp_dataset$Y_or; Z_or = temp_dataset$Z_or; rm(temp_dataset)
      # 5.b. Derive the compositional dataset in S_D 
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
      
      # 5.c. Compute ilr transformed data
      Y_transf = matrix(data=NA, nrow=n_y, ncol=D-Q-1)
      Z_transf = matrix(data=NA, nrow=n_z, ncol=D-1)
      for (id in 1:n_y){
        Y_transf[id,] <- ilr(Y_comp[id,], V=V_ilr_y)
      }
      for (id in 1:n_z){
        Z_transf[id,] <- ilr(Z_comp[id,], V=V_ilr_z)
      }
      
      results_list[[length(nz_all)*(iy-1)+iz]]$err_ilr_Y[isim] = max(abs(Y_transf-Y_or))
      results_list[[length(nz_all)*(iy-1)+iz]]$err_ilr_Z[isim] = max(abs(Z_transf-Z_or))

      # Step 6. Estimate sample value of the test statistic.
      aux = compute_stat_value(K_H0, Y_transf, Z_transf)
      temp_stat_value = aux$stat_value; 
      est_cov_y = aux$est_cov_y; est_cov_z = aux$est_cov_z;
      eigen_cov_y = aux$eigen_cov_y; eigen_cov_z = aux$eigen_cov_z;
      eigen_sum = aux$eigen_sum; rm(aux)
      
      # Step 7. Estimate p-value by using Schott's formula on the sample
      #         covariance matrix
      aux = compute_pvalue_schott(K_H0, temp_stat_value, est_cov_z, est_cov_y, 
                                eigen_cov_y, eigen_cov_z, eigen_sum, n_y, n_z)
      temp_pvalue_schott = aux$pvalue
      results_list[[aux_idx]]$est_mu_T[isim] = aux$est_mu_T
      results_list[[aux_idx]]$est_sigma2_T[isim] = aux$est_sigma2_T
      rm(aux)
      
      # Step 8. Estimate p-value by using Schott's formula on the true 
      #         covariance matrices
      aux = compute_pvalue_schott_theo(K_H0, temp_stat_value, omega_y, omega_z, 
                                             n_y, n_z)
      temp_pvalue_schott_theo = aux$pvalue
      results_list[[aux_idx]]$theo_mu_T[isim] = aux$mu_T
      results_list[[aux_idx]]$theo_sigma2_T[isim] = aux$sigma2_T
      rm(aux)
      
      
      # Step 9. Estimate p-value by using a bootstrap procedure
      aux = compute_pvalue_boot(K_H0, num_boot, temp_stat_value, 
                          eigen_cov_y, eigen_cov_z, Y_transf, Z_transf)
      temp_pvalue_boot = aux$pvalue; rm(aux)
      
      # Step ???. Store results
      results_list[[aux_idx]]$err_ilr_Y[isim] = max(abs(Y_transf-Y_or))
      results_list[[aux_idx]]$err_ilr_Z[isim] = max(abs(Z_transf-Z_or))
      
      results_list[[aux_idx]]$stat_values[isim] = temp_stat_value
      results_list[[aux_idx]]$pvalues_schott_est[isim] = temp_pvalue_schott
      results_list[[aux_idx]]$pvalues_schott_theo[isim] = temp_pvalue_schott_theo
      results_list[[aux_idx]]$pvalues_boot[isim] = temp_pvalue_boot
      
      rm(Y_or, Y_comp, Y_transf, Z_or, Z_comp, Z_transf, temp_stat_value)
      
      
    }
  }
}

# Step ??. Save result
result = list('param' = list('D'=D, 'Q'=Q, 'K_H0'=K_H0, 'K_sim'=K_sim,
                             "ny_all" = ny_all, "nz_all" = nz_all,
                             'distribution'=distribution, 'num_boot'=num_boot), 
              'results' = results_list, 'results_name' = results_list_name)
save(result, file=file.path(folder_res, 
                    paste0('res_K_', K_sim, '_distr_', distribution, '_newilr_corrpSchott.Rdata')))

