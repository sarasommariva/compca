# TODO:
# 1. Pulisci e commenta
# 2. Studia cosa fa la funzione randortho
# 3. Usare mvtnorm per estrarre da Gaussiane ed eliminare rockchalk?


#install.packages("pracma")
#install.packages("magic")

library(pracma)
library(magic)
library(mvtnorm)
library(rockchalk)

my_ilr <- function(comp){
  
  D = size(comp)[2]
  comp_transf = matrix(data=NA, nrow=1, ncol=D-1)
  for (i in 1:D-1){
    comp_transf[i] <- sqrt((i/(i+1))) * log( comp[i+1] / exp(mean(log(comp[1:i]))) )
  }
  
  return(comp_transf)
  
}

def_base_ilr <- function(D, transf='backpivot'){

  V_ilr<-matrix(ncol=D-1,nrow=D,0)

  if (transf=='pivot'){

    for(j in 1:(D-1)){
      V_ilr[j,j]<-(D-j)/(sqrt((D-j+1)*(D-j)))
      for(i in (j+1):D)
        V_ilr[i,j]<- -1/(sqrt((D-j+1)*(D-j)))
    }

    }else if(transf=='backpivot'){

      for(i in 1:(D-1)){
        aux <- sqrt(i/(i+1))
        for(j in 1:i){
          V_ilr[j,i] <- -1/i*aux}
        V_ilr[i+1, i] <- aux}

    }
  return(V_ilr)
}

def_projection <- function(M, Q1, Q2){
  
  EY = rbind(cbind(diag(M-1),          matrix(0, M-1, Q2)),
             cbind(matrix(0, Q2, M-1), diag(Q2)),
             matrix(0, Q1, M-1+Q2)
  )
  EZ = rbind(cbind(diag(M-1),          matrix(0, M-1, Q1)),
             matrix(0, Q2, M-1+Q1),
             cbind(matrix(0, Q1, M-1), diag(Q1))
  )
  
  return(list("EY"=EY, "EZ"=EZ))
  
}

rand_covmat <- function(M, Q1, Q2, K, eigenval_y, eigenval_z) {
  
  # Input:
  # D
  # K
  # Q
  # eigenval_y
  # eigenval_z
  
  # Step 1. Covariance matrix of Y
  mat_U_1 = randortho(M-1, type="orthonormal")
  if (Q2 > 0){
    mat_U_2 = randortho(Q2, type="orthonormal")}
  else{
    mat_U_2 = diag(0)  
    }
  mat_U = adiag(mat_U_1, mat_U_2)
  omega_y = mat_U %*% diag(eigenval_y) %*% t(mat_U)
  
  if (K == 0){
    # Step 2a. Independent covariance of Z
    mat_V_1 = randortho(M-1, type="orthonormal")
    mat_V_2 = randortho(Q1, type="orthonormal")
    mat_V = adiag(mat_V_1, mat_V_2)
  } else {
    # Step 2b. Covariance matrix of Z under H0(K)
    #     2b.1. Complete U to a basis of R^(D2-1)
    mat_U_comp = adiag(mat_U_1, diag(Q1))
    #   2b.2. Define the two rotation matrix R1 and R2
    det_R1 = -1
    while (det_R1 != 1){
      mat_R1 = randortho(K, type="orthonormal")
      det_R1 = round(det(mat_R1)) # Round avoids machine errors in while condition.
    }
    det_R2 = -1
    while (det_R2 != 1){
      mat_R2 = randortho(M+Q1-1-K, type="orthonormal")
      det_R2 = round(det(mat_R2))
    }
    #   2b.3. Define the covariance matrix
    mat_V = cbind(mat_U_comp[, 1:K] %*% mat_R1, mat_U_comp[, (K+1):(M+Q1-1)] %*% mat_R2)
    #mat_V_bis = mat_U_comp %*% adiag(mat_R1, mat_R2) # Alternative computation
  }
  omega_z = mat_V %*% diag(eigenval_z) %*% t(mat_V)
  # Output
  return(list("omega_y"=omega_y, "omega_z"=omega_z))
}


rand_covmat_single <- function(D, K, Q, eigenval_y, eigenval_z) {
  
  # Input:
  # D
  # K
  # Q
  # eigenval_y
  # eigenval_z
  
  # Step 1. Covariance matrix of the r.v. with structural zeros
  mat_U = randortho(D-Q-1, type="orthonormal")
  omega_y = mat_U %*% diag(eigenval_y) %*% t(mat_U)
  
  if (K == 0){
  # Step 2a. Independent covariance matrix of the complete r.v.
    mat_V = randortho(D-1, type="orthonormal")
  } else {
  # Step 2b. Covariance matrix of the complete r.v. under H0(K)
    #     2b.1. Complete U to a basis of R^(D-1)
    mat_U_comp = adiag(mat_U, diag(Q))
    #   2b.2. Define the two rotation matrix R1 and R2
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
    #   2b.3. Define the covariance matrix
    mat_V = cbind(mat_U_comp[, 1:K] %*% mat_R1, mat_U_comp[, (K+1):(D-1)] %*% mat_R2)
  }
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

rand_sample <- function(n_y, n_z, omega_y, omega_z, distribution, df=NULL){
  
  if (distribution=='Gaussian'){
    # Case 1. Gaussian variables (with zero mean)
    mu_y = rep(0, size(omega_y)[1])
    mu_z = rep(0, size(omega_z)[1])
    Y_or = mvrnorm(n=n_y, mu=mu_y, Sigma=omega_y)
    Z_or = mvrnorm(n=n_z, mu=mu_z, Sigma=omega_z)
  } else if (distribution=='Uniform'){
    # Case 2. Uniform distribution with ????
    L_y = chol(omega_y) # I use the upper triangular matrix 
    L_z = chol(omega_z)
    Y_or = t(replicate(n_y, runif(size(omega_y)[1], -sqrt(3), sqrt(3)))) %*% L_y
    Z_or = t(replicate(n_z, runif(size(omega_z)[1], -sqrt(3), sqrt(3)))) %*% L_z
  }else if (distribution=='Mult'){
    # Case 3. Multivariate Student t.
    scale = (df-2)/df;
    mu_y = rep(0, size(omega_y)[1])
    mu_z = rep(0, size(omega_z)[1])
    Y_or = rmvt(n=n_y, delta=mu_y, sigma=scale*omega_y, df = df)
    Z_or = rmvt(n=n_z, delta=mu_z, sigma=scale*omega_z, df = df)
  }else if (distribution=='Gaussian2'){
    # Case 4. Gaussian distribution using Cholesky decomposition
    L_y = chol(omega_y) # I use the upper triangular matrix 
    L_z = chol(omega_z)
    mu_y = rep(0, size(omega_y)[1])
    mu_z = rep(0, size(omega_z)[1])
    Y_or = mvrnorm(n=n_y, mu=mu_y, Sigma=eye(size(omega_y)[1])) %*% L_y
    Z_or = mvrnorm(n=n_z, mu=mu_z, Sigma=eye(size(omega_z)[1])) %*% L_z
  }else(print(paste('Cannot sample from', distribution)))
  
  # Output
  return(list("Y_or"=Y_or, "Z_or"=Z_or))
  
}

compute_stat_value <- function(M, Q1, Q2, K, Y_transf, Z_transf){

  projs = def_projection(M, Q1, Q2)
  EY = projs$EY; EZ = projs$EZ

  # Step 1. Estimate (unnormalized) sample covariance matrix
  est_cov_y = crossprod(sweep(Y_transf, 2, colMeans(Y_transf)))
  est_cov_z = crossprod(sweep(Z_transf, 2, colMeans(Z_transf)))
  
  # Step 2. Diagonalize sample covariance matrix 
  eigen_cov_y = eigen(est_cov_y)
  eigen_cov_z = eigen(est_cov_z)
  eigen_sum = eigen(EY %*% est_cov_y %*% t(EY) + EZ %*% est_cov_z %*% t(EZ))
  
  # Step 3. Compute test statistic 
  stat_value = sum(eigen_cov_y$values[1:K] + eigen_cov_z$values[1:K] 
                   - eigen_sum$values[1:K])
  # Output
  return(list('stat_value'=stat_value, 
              'est_cov_y'=est_cov_y, 'est_cov_z'=est_cov_z,
              'eigen_cov_y'=eigen_cov_y, 'eigen_cov_z'=eigen_cov_z, 
              'eigen_sum'=eigen_sum))
  
}

compute_stat_value_single <- function(K, Y_transf, Z_transf){
  
  Q = size(Z_transf)[2] - size(Y_transf)[2]
  
  # Step 1. Estimate (unnormalized) sample covariance matrix
  est_cov_y = crossprod(sweep(Y_transf, 2, colMeans(Y_transf)))
  est_cov_z = crossprod(sweep(Z_transf, 2, colMeans(Z_transf)))
  
  # Step 2. Diagonalize sample covariance matrix 
  eigen_cov_y = eigen(est_cov_y)
  eigen_cov_z = eigen(est_cov_z)
  eigen_sum = eigen(adiag(est_cov_y, matrix(0, Q, Q)) + est_cov_z)
  
  # Step 3. Compute test statistic 
  stat_value = sum(eigen_cov_y$values[1:K] + eigen_cov_z$values[1:K] 
                   - eigen_sum$values[1:K])
  # Output
  return(list('stat_value'=stat_value, 
              'est_cov_y'=est_cov_y, 'est_cov_z'=est_cov_z,
              'eigen_cov_y'=eigen_cov_y, 'eigen_cov_z'=eigen_cov_z, 
              'eigen_sum'=eigen_sum))
  
}

compute_param_nulldistr <- function(K, alpha, beta, psi, mat_U, mat_V, n_y, n_z){
  
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

compute_pvalue_schott_theo <- function(M, Q1, Q2, K, stat_value, omega_y, omega_z, 
                                       n_y, n_z){
  
  D = M + Q1 + Q2
  projs = def_projection(M, Q1, Q2)
  EY = projs$EY; EZ = projs$EZ
  
  # Step 1. Diagonalize covariance matrices
  eigen_omega_y = eigen(omega_y)
  eigen_omega_z = eigen(omega_z)
  pooled_omega = ((n_y-1)* EY %*% omega_y %*% t(EY) + 
                  (n_z-1)* EZ %*% omega_z %*% t(EZ)) / (n_y+n_z-2)
  eigen_pooled = eigen(pooled_omega)
  
  # Step 2. Collect eigen values
  alpha = as.vector(EY %*% eigen_omega_y$values)
  beta  = as.vector(EZ %*% eigen_omega_z$values)
  psi = eigen_pooled$values
  
  # Step 3. Rotate eigen vectors
  mat_K = eigen_pooled$vectors
  mat_U = t(mat_K) %*% adiag(eigen_omega_y$vectors, diag(Q1))
  temp_mat = EZ %*%  eigen_omega_z$vectors %*% t(EZ) + 
             adiag(diag(0, nrow=M-1), diag(1, nrow=Q2), diag(0, Q1)) 
  mat_V = t(mat_K) %*% temp_mat
  # Transform into a block diagonal matrix (this property is already 
  # satisfied only under H_0)
  mat_U = adiag(mat_U[1:K, 1:K], mat_U[(K+1):(D-1), (K+1):(D-1)])
  mat_V = adiag(mat_V[1:K, 1:K], mat_V[(K+1):(D-1), (K+1):(D-1)])
  
  # Step 4. Estimate parameters of the null distribution
  paramsT = compute_param_nulldistr(K, alpha, beta, psi, mat_U, mat_V, n_y, n_z)
  
  # Step 5. Compute p-value
  const_chi = 0.5*paramsT$sigma2_T / paramsT$mu_T
  df_chi = max(round(2*paramsT$mu_T^2/paramsT$sigma2_T), 1)

  if (const_chi > 0){lower_tail=FALSE}else{lower_tail=TRUE}
  pvalues = pchisq(stat_value/const_chi, df_chi, lower.tail=lower_tail)

  # Output
  return(list("pvalue"=pvalues, 
              "mu_T"=paramsT$mu_T, "sigma2_T"=paramsT$sigma2_T))
}

compute_pvalue_schott_theo_single <- function(K, stat_value, omega_y, omega_z, 
                                       n_y, n_z){
  
  Q = size(omega_z)[1] - size(omega_y)[1]
  D = size(omega_z)[1] + 1
  
  # Step 1. Diagonalize covariance matrices
  eigen_omega_y = eigen(omega_y)
  eigen_omega_z = eigen(omega_z)
  pooled_omega = ((n_y-1)*adiag(omega_y, matrix(0, Q, Q)) + (n_z-1)*omega_z) / (n_y+n_z-2)
  eigen_pooled = eigen(pooled_omega)
  
  # Step 2. Collect eigen values
  alpha = c(eigen_omega_y$values, matrix(0, Q, 1))
  beta = eigen_omega_z$values
  psi = eigen_pooled$values
  
  # Step 3. Rotate eigen vectors
  mat_K = eigen_pooled$vectors
  mat_U = t(mat_K) %*% adiag(eigen_omega_y$vectors, diag(Q))
  mat_V = t(mat_K) %*% eigen_omega_z$vectors
  # Transform into a block diagonal matrix (this property is already 
  # satisfied only under H_0)
  mat_U = adiag(mat_U[1:K, 1:K], mat_U[(K+1):(D-1), (K+1):(D-1)])
  mat_V = adiag(mat_V[1:K, 1:K], mat_V[(K+1):(D-1), (K+1):(D-1)])
  
  # Step 4. Estimate parameters of the null distribution
  paramsT = compute_param_nulldistr(K, alpha, beta, psi, mat_U, mat_V, n_y, n_z)
  
  # Step 5. Compute p-value
  const_chi = 0.5*paramsT$sigma2_T / paramsT$mu_T
  df_chi = max(round(2*paramsT$mu_T^2/paramsT$sigma2_T), 1)
  
  if (const_chi > 0){lower_tail=FALSE}else{lower_tail=TRUE}
  pvalues = pchisq(stat_value/const_chi, df_chi, lower.tail=lower_tail)
  
  # Output
  return(list("pvalue"=pvalues, 
              "mu_T"=paramsT$mu_T, "sigma2_T"=paramsT$sigma2_T))
}


compute_pvalue_schott <- function(M, Q1, Q2, K, stat_value, cov_z, cov_y, 
                                  eigen_cov_y, eigen_cov_z, eigen_sum, n_y, n_z){
  
  D = M + Q1 + Q2
  projs = def_projection(M, Q1, Q2)
  EY = projs$EY; EZ = projs$EZ
  
  #   Step 1. Normalize eigenvalue
  est_alpha = as.vector(EY %*% eigen_cov_y$values/(n_y-1))
  est_beta = as.vector(EZ %*% eigen_cov_z$values/(n_z -1))
  est_psi =  as.vector(eigen_sum$values / (n_y + n_z - 2))
  
  #  Step 2. Rotate eigenvectors
  # !!!! Controlla !!!!
  proj_cov_y = EY %*% cov_y %*% t(EY)
  proj_cov_z = EZ %*% cov_z %*% t(EZ)
  est_mat_K = eigen_sum$vectors
  temp1 = eigen(t(est_mat_K[,1:K]) %*% 
                  proj_cov_y %*% est_mat_K[,1:K])
  temp2 = eigen(t(est_mat_K[,(K+1):(D-1)]) %*% 
                  proj_cov_y %*% est_mat_K[,(K+1):(D-1)])
  est_mat_U = adiag(temp1$vectors, temp2$vectors); rm(temp1, temp2)
  temp1 = eigen(t(est_mat_K[,1:K]) %*% 
                  proj_cov_z %*% est_mat_K[,1:K])
  temp2 = eigen(t(est_mat_K[,(K+1):(D-1)]) %*% 
                  proj_cov_z %*% est_mat_K[,(K+1):(D-1)])
  est_mat_V = adiag(temp1$vectors, temp2$vectors); rm(temp1, temp2)
  
  #   Step 3. Estimate parameters of the null distribution of the test statistic
  est_paramsT = compute_param_nulldistr(K, est_alpha, est_beta, est_psi, 
                                        est_mat_U, est_mat_V, n_y, n_z)
  
  #   Step 4. Compute p-value
  est_const_chi = 0.5*est_paramsT$sigma2_T / est_paramsT$mu_T
  est_df_chi = max(round(2*est_paramsT$mu_T^2/est_paramsT$sigma2_T), 1)
  
  if (est_const_chi > 0){lower_tail=FALSE}else{lower_tail=TRUE}
  est_pvalues = pchisq(stat_value/est_const_chi, est_df_chi, lower.tail=lower_tail)
  
  # Output
  return(list("pvalue"=est_pvalues, 
              "est_mu_T"=est_paramsT$mu_T, "est_sigma2_T"=est_paramsT$sigma2_T))
  
  
}


compute_pvalue_schott_single <- function(K, stat_value, cov_z, cov_y, 
                                eigen_cov_y, eigen_cov_z, eigen_sum, n_y, n_z){
  
  Q = size(cov_z)[1] - size(cov_y)[1]
  D = size(cov_z)[1] + 1
  
  #   Step 1. Normalize eigenvalue
  est_alpha = c(eigen_cov_y$values/(n_y-1), matrix(0, Q, 1))
  est_beta = eigen_cov_z$values/(n_z -1)
  est_psi =  eigen_sum$values / (n_y + n_z - 2)
  
  #  Step 2. Rotate eigenvectors
  est_mat_K = eigen_sum$vectors
  temp1 = eigen(t(est_mat_K[,1:K]) %*% 
                  adiag(cov_y, matrix(0, Q, Q)) %*% est_mat_K[,1:K])
  temp2 = eigen(t(est_mat_K[,(K+1):(D-1)]) %*% 
                  adiag(cov_y, matrix(0, Q, Q)) %*% est_mat_K[,(K+1):(D-1)])
  est_mat_U = adiag(temp1$vectors, temp2$vectors); rm(temp1, temp2)
  temp1 = eigen(t(est_mat_K[,1:K]) %*% cov_z %*% est_mat_K[,1:K])
  temp2 = eigen(t(est_mat_K[,(K+1):(D-1)]) %*% cov_z %*% est_mat_K[,(K+1):(D-1)])
  est_mat_V = adiag(temp1$vectors, temp2$vectors); rm(temp1, temp2)
  
  #   Step 3. Estimate parameters of the null distribution of the test statistic
  est_paramsT = compute_param_nulldistr(K, est_alpha, est_beta, est_psi, 
                                        est_mat_U, est_mat_V, n_y, n_z)

  #   Step 4. Compute p-value
  est_const_chi = 0.5*est_paramsT$sigma2_T / est_paramsT$mu_T
  est_df_chi = max(round(2*est_paramsT$mu_T^2/est_paramsT$sigma2_T), 1)
  
  if (est_const_chi > 0){lower_tail=FALSE}else{lower_tail=TRUE}
  est_pvalues = pchisq(stat_value/est_const_chi, est_df_chi, lower.tail=lower_tail)
  
  # Output
  return(list("pvalue"=est_pvalues, 
              "est_mu_T"=est_paramsT$mu_T, "est_sigma2_T"=est_paramsT$sigma2_T))
  
  
}

compute_pvalue_boot <- function(M, Q1, Q2, K, num_rip, stat_value, 
                                eigen_cov_y, eigen_cov_z, Y_transf, Z_transf){
  D1 <- M+Q2; D2 <- M+Q1
  n_y = size(Y_transf)[1]; n_z = size(Z_transf)[1];
  
  # Step 1. Initialize
  boot_stat_values = matrix(data=NA, nrow=num_rip, ncol=1)
  
  # Step 2. Qr decomposition of the eigenvector of Y
  qr_firstvect = qr(eigen_cov_y$vectors[,1:K])
  #qr_firstvect = qr(eigen_cov_y$vectors[,1:M-1])   # Alternative
  qmat = qr.Q(qr_firstvect, complete=TRUE)
  rmat = qr.R(qr_firstvect)
  Yb_transf = Y_transf %*% qmat
  
  # Step 2. Randomly rotate the eigenvector of Y.
  Rb_1 = randortho(K, type="orthonormal")
  Rb_2 = randortho(D2-K-1, type="orthonormal")
  est_mat_Ub = adiag(rmat, diag(D2-1-K)) %*% adiag(Rb_1, Rb_2)
  #est_mat_Ub = adiag(rmat, diag(Q1)) %*% adiag(Rb_1, Rb_2) # Alternative
  # Step 3. Rotate the dataset Z so that its PCs coincide with those just definied
  Rb = est_mat_Ub %*% t(eigen_cov_z$vectors)
  Zb_transf = Z_transf %*% t(Rb)
  for (iboot in 1:num_rip){
  # Step 4. Define a bootstrap dataset
    Y_boot = Yb_transf[randi(n_y, n_y, 1), ]
    Z_boot = Zb_transf[randi(n_z, n_z, 1), ]
  # Step 5. Compute value of the test statistic for the bootstrap sample
    aux_boot = compute_stat_value(M, Q1, Q2, K_H0, Y_boot, Z_boot)
    boot_stat_values[iboot] = aux_boot$stat_value
  }
  # Step 6. Estimate the bootstrap p-value
  boot_pvalues = sum(boot_stat_values > stat_value)/num_rip
  
  # Output
  return(list("pvalue"=boot_pvalues, "stat_values"=boot_stat_values))
}

compute_pvalue_boot_single <- function(K, num_rip, stat_value, 
                                eigen_cov_y, eigen_cov_z, Y_transf, Z_transf){
  
  D = size(Z_transf)[2] + 1
  Q = D - (size(Y_transf)[2] + 1)
  n_y = size(Y_transf)[1]; n_z = size(Z_transf)[1];
  
  # Step 1. Initialize
  boot_stat_values = matrix(data=NA, nrow=num_rip, ncol=1)
  
  # Step 2. Randomly rotate the eigenvector of Y.
  Rb_1 = randortho(K, type="orthonormal")
  Rb_2 = randortho(D-K-1, type="orthonormal")
  est_mat_Ub = adiag(eigen_cov_y$vectors, diag(Q)) %*% adiag(Rb_1, Rb_2)
  # Step 3. Rotate the dataset Z so that its PCs coincide with those just definied
  Rb = est_mat_Ub %*% t(eigen_cov_z$vectors)
  Zb_transf = Z_transf %*% t(Rb)
  for (iboot in 1:num_rip){
    # Step 4. Define a bootstrap dataset
    Y_boot = Y_transf[randi(n_y, n_y, 1), ]
    Z_boot = Zb_transf[randi(n_z, n_z, 1), ]
    # Step 5. Compute value of the test statistic for the bootstrap sample
    aux_boot = compute_stat_value_single(K_H0, Y_boot, Z_boot)
    boot_stat_values[iboot] = aux_boot$stat_value
  }
  # Step 6. Estimate the bootstrap p-value
  boot_pvalues = sum(boot_stat_values > stat_value)/num_rip
  
  # Output
  return(list("pvalue"=boot_pvalues, "stat_values"=boot_stat_values))
}
  