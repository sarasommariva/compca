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

mid_vagina_se <- se[, meta$HMP_BODY_SUBSITE == "Mid Vagina"]
vag_intro_se <- se[, meta$HMP_BODY_SUBSITE == "Vaginal Introitus"]
post_forn_se <- se[, meta$HMP_BODY_SUBSITE == "Posterior Fornix"]


# --> Stool
stool_se <- se[, meta$HMP_BODY_SUBSITE == "Stool"]

# Coppie con cui ottengo buoni risultati:
#pop1 <- mid_vagina_se # Casi con sottospazi comuni
#pop2 <- post_forn_se

#pop1 <- mid_vagina_se
#pop2 <- vag_intro_se

#pop1 <- post_forn_se
#pop2 <- vag_intro_se

pop1 <- post_forn_se  # --> Caso senza sottospazio comune
pop2 <- stool_se


# Step 3. Data preprocessing
pop1_phylum <- make_phylum_table(pop1)
pop1_native <- native_support(pop1_phylum, 
                               prevalence = 0.05, mean_abundance = 1e-4)
pop2_phylum   <- make_phylum_table(pop2)
pop2_native   <- native_support(pop2_phylum, 
                               prevalence = 0.05, mean_abundance = 1e-4)


# Some checks
print("---> Original: \n")
cat("Population 1 samples:", nrow(pop1_phylum), "phyla:", ncol(pop1_phylum), "\n", 
    "Population 2 samples:", nrow(pop2_phylum), "phyla:", ncol(pop2_phylum), "\n")
cat("Population 1 phyla:\n")
print(colnames(pop1_phylum))
cat("Population 2 phyla:\n")
print(colnames(pop2_phylum))

print("---> After preprocessing: \n")
cat("Population 1 samples:", nrow(pop1_native), "phyla:", ncol(pop1_native), "\n", 
    "Population 2 samples:", nrow(pop2_native), "phyla:", ncol(pop2_native), "\n")
cat("Population 1 phyla:\n")
print(colnames(pop1_native))
cat("Population 2 phyla:\n")
print(colnames(pop2_native))

#write.csv(stool_native, "HMP_V35_stool_phylum_composition.csv")
#write.csv(vag_native, "HMP_V35_vaginal_phylum_composition.csv")


# --------------    PROVO A LANCIARE IL MIO CODICE -----------------------
# TODO:
# 1) Togliere la colonna 'Unclassified'
# 2) Gestire gli zeri non strutturali?

Y_or_temp = pop1_native
Z_or_temp = pop2_native

n_y = nrow(Y_or_temp)
n_z = nrow(Z_or_temp)
D = ncol(Z_or_temp)
Q = D - ncol(Y_or_temp)
K_H0_tested = seq(1, D-Q)

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

num_boot = 3000

pvalues.schott <- matrix(data=NA, nrow=length(K_H0_tested))
pvalues.boot <- matrix(data=NA, nrow=length(K_H0_tested))
values_stat <- matrix(data=NA, nrow=length(K_H0_tested)) 

idx <- 1

for (K_H0 in K_H0_tested){
  
  cat("----- Tested K = ", K_H0, " --------- \n")
  
  # Step 3. Compute sample value of the test statistic
  aux = compute_stat_value(K_H0, Y_transf, Z_transf)
  values_stat[idx] = aux$stat_value; 
  est_cov_y = aux$est_cov_y; est_cov_z = aux$est_cov_z;
  eigen_cov_y = aux$eigen_cov_y; eigen_cov_z = aux$eigen_cov_z;
  eigen_sum = aux$eigen_sum; rm(aux)
  
  # Step 4. Estimate p-value by using Schott's formula on the sample
  #         covariance matrix
  aux = compute_pvalue_schott(K_H0, values_stat[idx], est_cov_z, est_cov_y, 
                              eigen_cov_y, eigen_cov_z, eigen_sum, n_y, n_z)
  pvalues.schott[idx] = aux$pvalue;
  
  #print('Estimated mu_T')
  #print(aux$est_mu_T)
  #print('Estimated sigma2_T')
  #print(aux$est_sigma2_T)
  
  rm(aux)
  
  # Step 5. Estimate p-value by using Bootstrapping
  aux = compute_pvalue_boot(K_H0, num_boot, values_stat[idx], 
                            eigen_cov_y, eigen_cov_z, Y_transf, Z_transf)
  pvalues.boot[idx] = aux$pvalue; rm(aux)
  
  print('Test statistics')
  print(values_stat[idx] )
  print('Schott:')
  print(pvalues.schott[idx])
  print('Boot:')
  print(pvalues.boot[idx])
  
  idx = idx + 1
}








