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
                       "Student t (4 df)", 
                       "Student t (8 df)", 
                       "Student t (40 df)"), each = num_sim)))

desired_order <- c("Gaussian", "Student t (4 df)", "Student t (8 df)", 
                   "Student t (40 df)", "Uniform", "Schott")

linestyles = c("Gaussian" = "longdash", 
               "Uniform" = "dotdash",
               "Student t (4 df)" = "dashed", 
               "Student t (8 df)" = "twodash",
               "Student t (40 df)" = "dotted", 
               "Schott" = "solid")
temp_n = length(linestyles)
temp = brewer.pal(n = temp_n-1, name = "Set1")
linecolors = c("Gaussian" = temp[5], 
               "Uniform" = temp[2],
               "Student t (4 df)" = temp[3], 
               "Student t (8 df)" = temp[4],
               "Student t (40 df)" = temp[1], 
               "Schott" = "black")

const_chi = 0.5*result$schott.sigma2_T / result$schott.mu_T
df_chi = 2*result$schott.mu_T^2/result$schott.sigma2_T

p_cdf <- ggplot(df_stat_values, aes(x = value, color = distr, linetype = distr)) +
  labs(x = "Test statistics", y = "Cdf", color = "Distribution", linetype = "Distribution") +
  stat_function(fun = ~ pchisq(.x/const_chi, df_chi), size=0.5,
    aes(color = "Schott", linetype="Schott")) +
  stat_ecdf(geom = "step", size=0.7) +
  scale_linetype_manual(values = linestyles, breaks=desired_order) +
  scale_color_manual(values = linecolors, breaks=desired_order) +
  theme_minimal() + xlim(0, 60)

print(p_cdf)

ggsave('./figure/comparison_cdf.pdf', plot = p_cdf, 
       width = 14, height = 7, units = "cm")
