rm(list = ls())

library(HMP16SData)
library(compositions)
library(SummarizedExperiment)
library(gridExtra)
source("utils.R")
debugSource("utils_realdata.R")

# Step 1. Load data
se <- V35()
meta <- as.data.frame(colData(se))
subsites_group <- lapply(split(meta$HMP_BODY_SUBSITE, meta$HMP_BODY_SITE), unique) 

# Step 2. Define parameters for the analysis
toll_abundance = 0
toll_prevalence = 0.05

comparisons <- list(c("Right Antecubital Fossa", "Anterior Nares"), 
                    c("Stool", "Saliva"))

# Step 3:  Select specific subsites
ic <- 1
site1 <- comparisons[[ic]][1]
site2 <- comparisons[[ic]][2]

pop1_se <- se[, meta$HMP_BODY_SUBSITE == site1]
pop1_phylum <- make_phylum_table(pop1_se)
pop1_native <- native_support(pop1_phylum, 
                              prevalence = toll_prevalence, mean_abundance = toll_abundance)
pop2_se <- se[, meta$HMP_BODY_SUBSITE == site2]
pop2_phylum <- make_phylum_table(pop2_se)
pop2_native <- native_support(pop2_phylum, 
                              prevalence = toll_prevalence, mean_abundance = toll_abundance)

# Step 3. Reorder data so that the common variable appear first
common_vars = intersect(colnames(pop1_native), colnames(pop2_native))
pop1 <- pop1_native[, c(common_vars, setdiff(colnames(pop1_native), common_vars))]
pop2 <- pop2_native[, c(common_vars, setdiff(colnames(pop2_native), common_vars))]

print(site1)
print(colnames(pop1))
print(site2)
print(colnames(pop2))


# Step 4. Run analysis
#    4.1. Define required parameters
M <- length(common_vars)
Q1 <- ncol(pop2) - M
Q2 <- ncol(pop1) - M

K_H0_tested = seq(1, M-1)
num_boot <- 1000

n_y <- nrow(pop1)
n_z <- nrow(pop2)

D = M + Q1 + Q2
V_ilr_y = def_base_ilr(D-Q1, transf='backpivot')

V_ilr_z = def_base_ilr(D-Q2, transf='backpivot')

#    4.2. Ilr transform data
Y_or_temp = pop1
Z_or_temp = pop2

Y_or_temp[Y_or_temp==0] <- 10^-10*min(min(Y_or_temp[Y_or_temp>0]))
Z_or_temp[Z_or_temp==0] <- 10^-10*min(min(Z_or_temp[Z_or_temp>0]))


Y_transf = matrix(data=NA, nrow=n_y, ncol=D-Q1-1)
Z_transf = matrix(data=NA, nrow=n_z, ncol=D-Q2-1)
for (id in 1:n_y){
  Y_transf[id,] <- ilr(Y_or_temp[id,], V=V_ilr_y)
}
for (id in 1:n_z){
  Z_transf[id,] <- ilr(Z_or_temp[id,], V=V_ilr_z)
}

pvalues.schott <- matrix(data=NA, nrow=length(K_H0_tested))
pvalues.boot <- matrix(data=NA, nrow=length(K_H0_tested))
values_stat <- matrix(data=NA, nrow=length(K_H0_tested)) 

if (Q1 > 0 | Q2 > 0){
  idx <- 1
  for (K_H0 in K_H0_tested){
    
    # Step 4.3. Compute sample value of the test statistic
    aux = compute_stat_value(M, Q1, Q2, K_H0, Y_transf, Z_transf)
    values_stat[idx] = aux$stat_value; 
    est_cov_y = aux$est_cov_y; est_cov_z = aux$est_cov_z;
    eigen_cov_y = aux$eigen_cov_y; eigen_cov_z = aux$eigen_cov_z;
    eigen_sum = aux$eigen_sum; rm(aux)
    
    # Step 4.4. Estimate p-value by using Schott's formula on the sample
    #         covariance matrix
    aux = compute_pvalue_schott(M, Q1, Q2, K_H0, values_stat[idx], est_cov_z, est_cov_y, 
                                eigen_cov_y, eigen_cov_z, eigen_sum, n_y, n_z)
    pvalues.schott[idx] = aux$pvalue;
    
    rm(aux)
    
    # Step 4.5. Estimate p-value by using Bootstrapping
    aux = compute_pvalue_boot(M, Q1, Q2, K_H0, num_boot, values_stat[idx], 
                              eigen_cov_y, eigen_cov_z, Y_transf, Z_transf)
    pvalues.boot[idx] = aux$pvalue; rm(aux)
    
    idx = idx + 1
  }
}

# Test1: relationship with classical compositional PCA
pca_Y <- princomp(acomp(Y_or_temp))
pca_Z <- princomp(acomp(Z_or_temp))

Y_clr <- clr(Y_or_temp)
est_cov_y_clr = crossprod(sweep(Y_clr, 2, colMeans(Y_clr)))/(n_y-1)
eigen_cov_y_clr = eigen(est_cov_y_clr)

print(pca_Y$sdev^2)
print(eigen_cov_y$values/(n_y-1))
print(eigen_cov_y_clr$values)

load_eigen_clr <- eigen_cov_y_clr$vectors
load_pca <- as.matrix(pca_Y$loadings)
load_eigen_ilr <- V_ilr_y %*% eigen_cov_y$vectors
weighted_load_eigen_ilr <- load_eigen_ilr %*% diag(sqrt(eigen_cov_y$values/(n_y-1)))

biplot(pca_Y)

#plot(pca_Y, type = "loadings")
par(mfrow = c(1, 2))
bpy <- biplot_loadings(unclass(pca_Y$loadings)[, 1:2], 
                       paste0('Loadings ', site1))
bpz <- biplot_loadings(unclass(pca_Z$loadings)[, 1:2], 
                       paste0('Loadings ', site2))

# Un po' di prove
#biplot_loadings(unclass(pca_Y$loadings)[, 1:2], 'Loadings Y')
#biplot_loadings(load_eigen_ilr[, 1:2], 'Loadings Y - prod')
#biplot_loadings(weighted_load_eigen_ilr[, 1:2], 'Loadings Y - prod weights')


#library(factoextra)
#p_pcaY <- fviz_pca_var(pca_Y, select.var = list(contrib = 3))
#p_pcaZ <- fviz_pca_var(pca_Z, select.var = list(contrib = 3))

#p_pcaY + xlim(-1.5, 1.5) + ylim(-1.5, 1.5)

#print(p_pcaY)
#print(p_pcaZ)
