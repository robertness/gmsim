library(gmsim)
# A big tree-like graph --------------------------
# A MVN Network Simulated with --------------------------
num_vars <- 40
num_obs <- 1000
## Mean 0 data from a random MVN covariance matrix
net <- simGaussianNet(num_vars, c(.5, 1.5), "ordered")
## Data simulated directly from the Bayesian network.  
boot <- rbn(net, n = num.obs) %>% 
  boot.strength(R = 500, m = num_obs,
                algorithm = "hc", 
                algorithm.args = list(score = "bic-g"))
ordered_boot <- list(truth = net, boot = boot)
devtools::use_data(ordered_boot, overwrite = TRUE)

