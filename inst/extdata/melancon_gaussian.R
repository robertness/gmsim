library(gmsim)
# A MVN Network Simulated with --------------------------
num_vars <- 15
num_obs <- 1000
## Mean 0 data from a random MVN covariance matrix
net <- simGaussianNet(num_vars, c(.5, 1.5), "melancon")
## Data simulated directly from the Bayesian network.  
boot <- rbn(net, n = num_obs) %>% 
  boot.strength(R = 500, m = num_obs,
                algorithm = "hc", 
                algorithm.args = list(score = "bic-g"))
melancon_boot <- list(truth = net, boot = boot)
devtools::use_data(melancon_boot)

