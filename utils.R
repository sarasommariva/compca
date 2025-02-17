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
  