# Todo:
# Provare a creare le tabelle direttamente via R

rm(list=ls())
setwd(getSrcDirectory(function(){})[1])

# Step 1. Parameters definition
#  1.a. Parameters for loading
folder_res = file.path('.', 'results')
K_sim = 2
distribution = 'Gaussian'
#  1.b. Parameters for analysing
alpha_sign = 0.05

# Step 1. Load data

#   1.b. Load
load(file.path(folder_res, 
                paste0('res_K_', K_sim, '_distr_', distribution, '_2.Rdata')))

ny_all= result$param$ny_all  # Read from file!!!!
nz_all = result$param$nz_all

num_sim = length(result$results[[1]]$pvalues_schott_est)

# Step 2. Compute power
# 2.a. Initialize 
power_schott_theo = matrix(data=NA, nrow=length(nz_all), ncol=length(ny_all))
power_schott = matrix(data=NA, nrow=length(nz_all), ncol=length(ny_all))
power_boot = matrix(data=NA, nrow=length(nz_all), ncol=length(ny_all))
power_fortex = matrix(data=NA, nrow=3*length(nz_all), ncol=length(ny_all))

for (iy in 1:length(ny_all)){
  for (iz in 1:length(nz_all)){
    aux_idx = length(nz_all)*(iy-1)+iz;
    power_schott_theo[iz, iy] = sum(
      result$results[[aux_idx]]$pvalues_schott_theo < alpha_sign)/num_sim
    power_schott[iz, iy] = sum(
      result$results[[aux_idx]]$pvalues_schott_est < alpha_sign)/num_sim
    power_boot[iz, iy] = sum(
      result$results[[aux_idx]]$pvalues_boot < alpha_sign)/num_sim
    power_fortex[3*(iz-1)+1, iy] = power_schott_theo[iz, iy]
    power_fortex[3*(iz-1)+2, iy] = power_schott[iz, iy]
    power_fortex[3*(iz-1)+3, iy] = power_boot[iz, iy]
  }
}

print(power_fortex)

# Some plots for start checking
aux_ny = 100; aux_nz = 100;
iy = match(aux_ny, ny_all); iz = match(aux_nz, nz_all)
aux_ris = result$results[[length(nz_all)*(iy-1)+iz]]

plot(aux_ris$theo_mu_T, aux_ris$pvalues_schott_theo)
title(result$results_name[[length(nz_all)*(iy-1)+iz]])

aux_ny = 100; aux_nz = 100;
iy = match(aux_ny, ny_all); iz = match(aux_nz, nz_all)
hist(result$results[[length(nz_all)*(iy-1)+iz]]$pvalues_schott_est, ylim=c(0, 120))

aux_ny = 60; aux_nz = 60;
iy = match(aux_ny, ny_all); iz = match(aux_nz, nz_all)
hist(result$results[[length(nz_all)*(iy-1)+iz]]$pvalues_schott_est, ylim=c(0, 120))

aux_ny = 20; aux_nz = 20;
iy = match(aux_ny, ny_all); iz = match(aux_nz, nz_all)
hist(result$results[[length(nz_all)*(iy-1)+iz]]$pvalues_schott_est, ylim=c(0, 120))