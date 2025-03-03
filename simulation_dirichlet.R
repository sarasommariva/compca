#install.packages("Compositional")
library(Compositional)

num_sim = 100
D = 6
Q = 2
n_y = 100
n_z = n_y

alpha_y = matrix( c(3, 30, 45, 45, 20 ,20, 20, 4), 
                   byrow = TRUE, ncol = D-Q) 
prob_y = c(0.5, 0.5) 
alpha_z = matrix(c(3, 30, 45, 45, 12, 3, 20 ,20, 20, 4, 40, 5), 
                  byrow = TRUE, ncol = D)
prob_z = c(0.5, 0.5)

dataset_all = array(0, dim = c(n_y+n_z, D, num_sim))

for (isim in 1:num_sim){
  Y_or = rmixdiri(n_y, alpha_y, prob_y)
  Y_or = Y_or$x
  Z_or = rmixdiri(n_z, alpha_z, prob_z)
  Z_or = Z_or$x
  
  dataset_all[1:n_y,1:(D-Q),isim] = Y_or
  dataset_all[(n_y+1):(n_y+n_z),,isim] = Z_or
}

data = list('param' = list('n_y' = n_y, 'n_z' = n_z, 'Q'=Q, 
                  'alpha_y'=alpha_y, 'alpha_z'=alpha_z, 
                  'prob_y'=prob_y, 'prob_z'=prob_z),
            'dataset_all'=dataset_all)
save(data, file="prova_dataset_dirichlet.Rdata")

# TODO:
# 1. Scegliere K diverso da D-Q
# 2. Scegliere due matrici di covarianze che soddisfino H0, ma siano diverse

