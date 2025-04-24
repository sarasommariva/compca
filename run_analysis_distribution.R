rm(list=ls())
setwd(getSrcDirectory(function(){})[1])

library(ggplot2)
library(RColorBrewer)

# Step 1. Parameters definition.
#  1.a. Parameters for loading
folder_res = file.path('.', 'results')
n_y = 100; n_z = 100;

# Step 2. Load data
load(file=file.path(folder_res, 
                      paste0('empirical_distr_ny', n_y, '_nz', n_z,'.Rdata')))
num_sim = length(result$stat_values_gaussian)

# Plot 1. Comparison of the cdf
df_stat_values = data.frame(
  value = c(result$stat_values_gaussian,
            result$stat_values_uniform,
            result$stat_values_mult4, 
            result$stat_values_mult8, 
            result$stat_values_mult40), 
  distr = factor(rep(c("Gaussian",
                       "Uniform",
                       "Multivariat t (df=4)", 
                       "Multivariat t (df=8)", 
                       "Multivariat t (df=40)"), each = num_sim)))
linestyles = c("Gaussian" = "longdash", 
               "Uniform" = "dotdash",
               "Multivariat t (df=4)" = "dashed", 
               "Multivariat t (df=8)" = "twodash",
               "Multivariat t (df=40)" = "dotted", 
               "Schott" = "solid")
temp_n = length(linestyles)
temp = brewer.pal(n = temp_n-1, name = "Set1")
linecolors = c("Gaussian" = temp[5], 
               "Uniform" = temp[2],
               "Multivariat t (df=4)" = temp[3], 
               "Multivariat t (df=8)" = temp[4],
               "Multivariat t (df=40)" = temp[1], 
               "Schott" = "black")

const_chi = 0.5*result$schott.sigma2_T / result$schott.mu_T
df_chi = 2*result$schott.mu_T^2/result$schott.sigma2_T

p_cdf <- ggplot(df_stat_values, aes(x = value, color = distr, linetype = distr)) +
  labs(x = "Test statistics", y = "Cdf", color = "Distribution", linetype = "Distribution") +
  stat_function(fun = ~ pchisq(.x/const_chi, df_chi), size=0.5,
    aes(color = "Schott", linetype="Schott")) +
  stat_ecdf(geom = "step", size=0.7) +
  scale_linetype_manual(values = linestyles) +
  scale_color_manual(values = linecolors) +
  theme_minimal() + xlim(0, 60)

ggsave('./figure/comparison_cdf.pdf', plot = p_cdf, 
       width = 14, height = 7, units = "cm")
