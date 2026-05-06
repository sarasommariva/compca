rm(list = ls())

# Step 0: Install required package and data
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("HMP16SData", "SummarizedExperiment"), ask = FALSE)

# Function for preparing table data
make_phylum_table <- function(x, min_reads = 1000) {
  counts <- assay(x, "16SrRNA")
  tax <- as.data.frame(rowData(x))
  
  phylum <- as.character(tax$PHYLUM)
  phylum[is.na(phylum) | phylum == "" | phylum == "NA"] <- "Unclassified"
  
  # aggregazione OTU -> phylum
  phylum_counts <- rowsum(counts, group = phylum, reorder = FALSE)
  
  # rimuovi campioni poveri
  keep_samples <- colSums(phylum_counts) >= min_reads
  phylum_counts <- phylum_counts[, keep_samples, drop = FALSE]
  
  # supporto nativo: tieni solo phyla osservati in quel dataset
  keep_phyla <- rowSums(phylum_counts) > 0
  phylum_counts <- phylum_counts[keep_phyla, , drop = FALSE]
  
  # composizioni campioni x phyla
  comp <- t(t(phylum_counts) / colSums(phylum_counts))
  comp <- t(comp)
  
  return(comp)
}

native_support <- function(x, prevalence = 0.05, mean_abundance = 1e-4) {
  keep <- colMeans(x > 0) >= prevalence & colMeans(x) >= mean_abundance
  x_sub <- x[, keep, drop = FALSE]
  x_sub <- x_sub / rowSums(x_sub)
  return(x_sub)
}


#---------------------------------------------------------------

# Step 1: Load libraries and data
library(HMP16SData)
library(compositions)
library(SummarizedExperiment)
source("utils.R")

# scegli V35; alternativa: V13()
se <- V35()

meta <- as.data.frame(colData(se))
table(meta$HMP_BODY_SUBSITE)

# Step 2: Select sites

# --> All hmp vaginal subsites
vaginal_subsites <- c(
  # "Anterior Nares", # NON usare; solo controllo visivo
  "Mid Vagina",
  "Posterior Fornix",
  "Vaginal Introitus"
)
vaginal_subsites <- intersect(
  c("Mid Vagina", "Posterior Fornix", "Vaginal Introitus"),
  unique(meta$HMP_BODY_SUBSITE)
)
vag_se   <- se[, meta$HMP_BODY_SUBSITE %in% vaginal_subsites]

# --> Stool
stool_se <- se[, meta$HMP_BODY_SUBSITE == "Stool"]

# Step 3. Data preprocessing
stool_phylum <- make_phylum_table(stool_se)
stool_native <- native_support(stool_phylum, 
                               prevalence = 0.05, mean_abundance = 1e-4)
vag_phylum   <- make_phylum_table(vag_se)
vag_native   <- native_support(vag_phylum, 
                               prevalence = 0.05, mean_abundance = 1e-4)


# Some checks
cat("Stool samples:", nrow(stool_phylum), "phyla:", ncol(stool_phylum), "\n", 
    "Vaginal samples:", nrow(vag_phylum), "phyla:", ncol(vag_phylum), "\n")
cat("Stool phyla:\n")
print(colnames(stool_phylum))
cat("Vaginal phyla:\n")
print(colnames(vag_phylum))



dim(stool_native)
dim(vag_native)
colnames(stool_native)
colnames(vag_native)

#write.csv(stool_native, "HMP_V35_stool_phylum_composition.csv")
#write.csv(vag_native, "HMP_V35_vaginal_phylum_composition.csv")


# --------------    PROVO A LANCIARE IL MIO CODICE -----------------------
# TODO:
# 1) Togliere la colonna 'Unclassified'
# 2) Gestire gli zeri non strutturali?

Y_or_temp = vag_native
Z_or_temp = stool_native

n_y = nrow(Y_or_temp)
n_z = nrow(Z_or_temp)
D = ncol(Z_or_temp)
Q = D - ncol(Y_or_temp)

Y_or = clo(Y_or_temp)[,1:(D-Q)]
Z_or = clo(Z_or_temp)

V_ilr_y = def_base_ilr(D-Q)
V_ilr_z = def_base_ilr(D)

Y_transf = matrix(data=NA, nrow=n_y, ncol=D-Q-1)
Z_transf = matrix(data=NA, nrow=n_z, ncol=D-1)
for (id in 1:n_y){
  Y_transf[id,] <- ilr(Y_or[id,], V=V_ilr_y)
}
for (id in 1:n_z){
  Z_transf[id,] <- ilr(Z_or[id,], V=V_ilr_z)
}

K_H0 = 6
num_boot = 3000

# Step 3. Compute sample value of the test statistic
aux = compute_stat_value(K_H0, Y_transf, Z_transf)
stat_value = aux$stat_value; 
est_cov_y = aux$est_cov_y; est_cov_z = aux$est_cov_z;
eigen_cov_y = aux$eigen_cov_y; eigen_cov_z = aux$eigen_cov_z;
eigen_sum = aux$eigen_sum; rm(aux)

# Step 4. Estimate p-value by using Schott's formula on the sample
#         covariance matrix
aux = compute_pvalue_schott(K_H0, stat_value, est_cov_z, est_cov_y, 
                            eigen_cov_y, eigen_cov_z, eigen_sum, n_y, n_z)
pvalue_schott = aux$pvalue;
print('Estimated mu_T')
print(aux$est_mu_T)
print('Estimated sigma2_T')
print(aux$est_sigma2_T)
rm(aux)

# Step 5. Estimate p-value by using Schott's formula on the sample
#         covariance matrix
aux = compute_pvalue_boot(K_H0, num_boot, stat_value, 
                          eigen_cov_y, eigen_cov_z, Y_transf, Z_transf)
pvalue_boot = aux$pvalue; rm(aux)

print('Test statistics')
print(stat_value)
print('Schott:')
print(pvalue_schott)
print('Boot:')
print(pvalue_boot)
