library(zCompositions)

make_phylum_table <- function(x, min_reads = 1000, remove_uncl = TRUE) {
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
  if (remove_uncl){
    phylum_counts <- phylum_counts[rownames(phylum_counts) != "Unclassified", , 
                                   drop = FALSE]
  }
  
  #phylum_counts <- cmultRepl(phylum_counts, method = "CZM")
  
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


biplot_loadings <- function(L, title_text){
  plot(L, type = "n", asp = 1,
     xlim = range(L[,1], 0), ylim = range(L[,2], 0),
     xlab = "Comp.1", ylab = "Comp.2")
abline(h = 0, v = 0, lty = 2, col = "grey")
arrows(0, 0, L[,1], L[,2], length = 0.1, col = "blue")
text(L[,1], L[,2], labels = rownames(L), pos = 3, col = "blue")
title(title_text)
}
