# TODO:
# 1. Pulisci e commenta
# 2. Studia cosa fa la funzione randortho


#install.packages("pracma")
#install.packages("magic")

library(pracma)
library(magic)

rand_covmat <- function(D, K, Q, eigenval_y, eigenval_z) {
  
  # Input:
  # D
  # K
  # Q
  # eigenval_y
  # eigenval_z
  
  # Step 1. Covariance matrix of the r.v. with structural zeros
  mat_U = randortho(D-Q-1, type="orthonormal")
  omega_y = mat_U %*% diag(eigenval_y) %*% t(mat_U)
  
  # Step 2. Covariance matrix of the complete r.v.
  #     2.1. Complete U to a basis of R^(D-1)
  mat_U_comp = adiag(mat_U, diag(Q))
  
  #     2.2. Define the two rotation matrix R1 and R2
  det_R1 = -1
  while (det_R1 != 1){
    mat_R1 = randortho(K, type="orthonormal")
    det_R1 = round(det(mat_R1)) # Round avoide machine errors in while condition.
  }
  det_R2 = -1
  while (det_R2 != 1){
    mat_R2 = randortho(D-1-K, type="orthonormal")
    det_R2 = round(det(mat_R2))
  }
  
  #     2.3. Define the covariance matrix
  mat_V = cbind(mat_U_comp[, 1:K] %*% mat_R1, mat_U_comp[, (K+1):(D-1)] %*% mat_R2)
  omega_z = mat_V %*% diag(eigenval_z) %*% t(mat_V)
  
  # Output
  return(list("omega_y"=omega_y, "omega_z"=omega_z))
}


rand_indep_covmat <- function(D, K, Q, eigenval_y, eigenval_z) {

  # Step 1. Covariance matrix of the r.v. with structural zeros
  mat_U = randortho(D-Q-1, type="orthonormal")
  omega_y = mat_U %*% diag(eigenval_y) %*% t(mat_U)
  
  # Step 1. Covariance matrix of the complete r.v.
  mat_V = randortho(D-1, type="orthonormal")
  omega_z = mat_V %*% diag(eigenval_z) %*% t(mat_V)
  
  # Output
  return(list("omega_y"=omega_y, "omega_z"=omega_z))
}

compute_param_nulldistr <- function(alpha, beta, psi, mat_U, mat_V, n_y, n_z){
  
  mat_UlamU = mat_U %*% diag(alpha) %*% t(mat_U)
  mat_VlamV = mat_V %*% diag(beta) %*% t(mat_V)
  
  # Step 1. Compute mean of the test statistic
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
  
  # Step 2. Compute variance of the test statistic
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
            ((n_y-1)*mat_UlamU[ii,hh]*mat_UlamU[jj,ll]+(n_z-1)*mat_VlamV[ii,hh]*mat_VlamV[jj,ll])^2/(psi[hh]-psi[ll])
        }
      }
      sigma2_T = sigma2_T + 
        2 * ( (alpha[ii] * alpha[jj])/(alpha[ii] - alpha[jj]) )^2 +
        2 * ( (beta[ii] * beta[jj])/(beta[ii] - beta[jj]) )^2 -
        4 * ((n_y+n_z-2)*(psi[ii]-psi[jj]))^-1 * sum1 + 
        2 * ((n_y+n_z-2)^2*(psi[ii]-psi[jj]))^-1 * sum2
    }
  }

  return(list("mu_T"=mu_T, "sigma2_T"=sigma2_T))  
}
  