# TODO:
# 1. Rendere questo codice una funzione.
#    O cmq essere sicura che i parametri che passo in input ai diversi codici
#    siano gli stessi
# 2. Aggiungere un check che le matrici U e V siano diagonali a blocchi.
#    Altrimenti va cambiata un po' la definizione di mat_UlamU e mat_VlamV

rm(list=ls())
setwd('~/Documenti/compoteam/paper_stattest/compca')

load('prova.Rdata')

# Step 1. Definition of the parameters
D = dim(result$param$omega_z)[1] + 1
n_y = result$param$n_y
n_z = result$param$n_z
omega_y = result$param$omega_y
omega_z = result$param$omega_z

K = 2 # <--- READ FROM PARAM!!!!!

# Step 2. Diagonalization of the covariance matrices
pooled_omega = ((n_y-1)*omega_y + (n_z-1)*omega_z) / (n_y+n_z-2)

eigen_omega_y = eigen(omega_y)
eigen_omega_z = eigen(omega_z)
eigen_pooled = eigen(pooled_omega)

alpha = eigen_omega_y$values
beta = eigen_omega_z$values
psi = eigen_pooled$values
mat_K = eigen_pooled$vectors
mat_U = t(mat_K) %*% eigen_omega_y$vectors
mat_V = t(mat_K) %*% eigen_omega_z$vectors
mat_UlamU = mat_U %*% diag(alpha) %*% t(mat_U)
mat_VlamV = mat_V %*% diag(beta) %*% t(mat_V)

# Step 3. Compute mean of the test statistic
mu_T = 0
for (ii in 1:K){
  for (jj in (K+1):(D-1)){
      mu_T = mu_T + 
        (alpha[ii] * alpha[jj])/(alpha[ii] - alpha[jj]) +
        (beta[ii] * beta[jj])/(beta[ii] - beta[jj]) -
        ((n_y+n_z-2)*(psi[ii]-psi[jj]))^-1 * 
          ( (n_y-1)*mat_UlamU[ii,ii]*mat_UlamU[jj,jj] +  (n_z-1)*mat_VlamV[ii,ii]*mat_VlamV[jj,jj])
  }
}

# Step 4. Compute variance of the test statistic
sigma2_T = 0
for (ii in 1:K){
  for (jj in (K+1):(D-1)){
    sum1 = 0
    sum2 = 0
    for (hh in 1:K){
      for (ll in (K+1):(D-1)){
        sum1 = sum1 + 
          (n_y-1)*(mat_U[ii,hh]*mat_U[jj,ll]*alpha[hh]*alpha[ll])^2/(alpha[hh]-alpha[ll]) +
          (n_z-1)*(mat_V[ii,hh]*mat_V[jj,ll]*beta[hh]*beta[ll])^2/(beta[hh]-beta[ll])
        sum2 = sum2 + 
          ((n_y-1)*mat_UlamU[ii,hh]*mat_UlamU[jj,ll]+(n_z-1)*mat_VlamV[ii,hh]*mat_VlamV[jj,ll])/(psi[hh]-psi[ll])
      }
    }
    sigma2_T = sigma2_T + 
      2 * ( (alpha[ii] * alpha[jj])/(alpha[ii] - alpha[jj]) )^2 +
      2 * ( (beta[ii] * beta[jj])/(beta[ii] - beta[jj]) )^2 -
      4 * ((n_y+n_z-2)*(psi[ii]-psi[jj]))^-1 * sum1 + 
      2 * ((n_y+n_z-2)^2*(psi[ii]-psi[jj]))^-1 * sum2
  }
}

# Step 5. Define parameters of the Chi
const_chi = 0.5*sigma2_T / mu_T
df_chi = floor(2*mu_T^2/sigma2_T)

delta_bin = 0.1
points = seq(min(result$stat_values)-delta_bin, max(result$stat_values)+delta_bin, by=delta_bin)

dist_chi = const_chi^-1 * dchisq(points/const_chi, df_chi)

h = hist(result$stat_values, breaks=300, freq=FALSE)
lines(points, dist_chi)





