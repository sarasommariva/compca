# TODO:
# 1. Rendere questo codice una funzione.
#    O cmq essere sicura che i parametri che passo in input ai diversi codici
#    siano gli stessi
# 2. Aggiungere un check che le matrici U e V siano diagonali a blocchi.
#    Altrimenti va cambiata un po' la definizione di mat_UlamU e mat_VlamV

rm(list=ls())
setwd('~/Documenti/compoteam/paper_stattest/compca')

source("utils.R")

load('prova_h0.Rdata')

# Step 1. Definition of the parameters
D = dim(result$param$omega_z)[1] + 1
K = result$param$K
#Q = result$param$Q
Q = 2
n_y = result$param$n_y
n_z = result$param$n_z
omega_y = result$param$omega_y
omega_z = result$param$omega_z

# Step 2. Diagonalization of the covariance matrices
eigen_omega_y = eigen(omega_y)
eigen_omega_z = eigen(omega_z)

#   - Pooled covariance matrix
pooled_omega = ((n_y-1)*adiag(omega_y, matrix(0, Q, Q)) + (n_z-1)*omega_z) / (n_y+n_z-2)
eigen_pooled = eigen(pooled_omega)

alpha = c(eigen_omega_y$values, matrix(0, Q, 1))
beta = eigen_omega_z$values
psi = eigen_pooled$values
mat_K = eigen_pooled$vectors
mat_U = t(mat_K) %*% adiag(eigen_omega_y$vectors, diag(Q))
mat_V = t(mat_K) %*% eigen_omega_z$vectors

# Step 3. Compute mean and variance of the test statistics
paramsT = compute_param_nulldistr(alpha, beta, psi, mat_U, mat_V, n_y, n_z)
mu_T = paramsT$mu_T
sigma2_T = paramsT$sigma2_T
rm(paramsT)

# Compute also sample mean and variance for comparison
mu_T_camp = mean(result$stat_values)
sigma2_T_camp = var(result$stat_values)

# Step 5. Define parameters of the Chi
const_chi = 0.5*sigma2_T / mu_T
df_chi = 2*mu_T^2/sigma2_T

#   Plots
# P1. Histogram of the sample values against the theoretical distribution
delta_bin = 0.1
points = seq(min(result$stat_values)-delta_bin, max(result$stat_values)+delta_bin, by=delta_bin)
dist_chi = const_chi^-1 * dchisq(points/const_chi, df_chi)
h = hist(result$stat_values, breaks=300, freq=FALSE)
lines(points, dist_chi)

# P2. 





