# -----------------------------------------------------------
# Program: example3.r
# Details:
#           Section 4.2 – Scenario 3
#
#   Estimates are obtained using a Bayesian Hierarchical Model (BHM),
#   resulting in parameters u₁, u₂, r₁, r₂, s₁, and s₂.
#
#   Non-informative priors are assumed for the BHM hyperparameters:
#     μ₀ = 0.1, η = 0.1, ξ = 0.1,
#   with a small amount of borrowing across clinics (τ₀ = 5).
# -----------------------------------------------------------


# please set the current folder as working directory


library(rstan)
library(tidyr)

fullBayesian<-function(x,group,mu0,tau0,eta,xi,seed=123){
  # Prepare data for Stan
  stan_data <- list(
    N = length(x),
    K = length(unique(group)),
    group = group,
    x = x,
    mu0 = mu0,
    tau0 = tau0,
    eta = eta,
    xi = xi
  )
  
  # Compile and run
  fit <- stan(
    file = "fully_bayesian_model.stan",
    data = stan_data,
    iter = 2000,
    warmup = 500,
    chains = 4,
    seed = seed
  )
  return(fit)
}


#==================================================================

# Hyperprior settings
mu0 <- 0.1     # prior mean for u
tau0 <- 5    # prior std for u
eta  <- 0.1  # scale for half-Cauchy on r
xi   <- 0.1  # scale for half-Cauchy on s
# real data
x_ca <- c(0.23, 0.25, 0.28,0.19, 0.22, 0.23, 0.27, 0.21, 0.22, 0.24, 0.28)
group_ca <- c(rep(1,3),rep(2,4),rep(3,4))

x_mo <- c(0.20, 0.29, 0.19, 0.28, 0.34)
group_mo <- c(rep(1,2),rep(2,3))

fullBayesian(x_ca,group_ca,mu0,tau0,eta,xi)->mcmc_ca
fullBayesian(x_mo,group_mo,mu0,tau0,eta,xi)->mcmc_mo

print(mcmc_ca, pars = c("m", "s", "u", "r"))
print(mcmc_mo, pars = c("m", "s", "u", "r"))



# The following code is just for visualization.

#=================================================================
# visualize how the estimated group means m_i are "shrunk" toward 
# the global mean u in the Bayesian hierarchical model.
#=================================================================



library(ggplot2)
library(dplyr)

shrinkage_plot <- function(x,group_x,fit){
  
  
  # Step 1: Extract posterior draws
  posterior_df <- as.data.frame(fit)
  K <- length(unique(group_x))
  m_names <- paste0("m[", 1:K, "]")
  # Extract only m[1], ..., m[K]
  m_draw <- posterior_df[, m_names]
  # Compute posterior summary for m_i
  post_summary <- posterior_df %>%
    select(all_of(m_names)) %>%
    pivot_longer(everything(), names_to = "group", values_to = "m_draw") %>%
    group_by(group) %>%
    summarise(
      m_mean = mean(m_draw),
      m_lower = quantile(m_draw, 0.025),
      m_upper = quantile(m_draw, 0.975)
    ) %>%
    mutate(group_num = as.numeric(gsub("m\\[|\\]", "", group)))
  
  # Step 2: Compute observed group means
  obs_means <- tapply(x, group_x, mean)
  df_obs <- data.frame(
    group_num = 1:K,
    xbar = obs_means
  )
  
  # Step 3: Merge for plotting and prepare lines
  df_plot <- left_join(post_summary, df_obs, by = "group_num")
  
  # Step 4: Global mean (posterior mean of u)
  post_u <- mean(posterior_df$u)
  
  # Step 5: Plot
  g = ggplot(df_plot, aes(x = factor(group_num))) +
    # Dashed line showing shrinkage from observed mean to posterior mean
    geom_segment(aes(x = factor(group_num), xend = factor(group_num),
                     y = xbar, yend = m_mean), linetype = "dashed", color = "gray") +
    
    # Observed means (raw)
    geom_point(aes(y = xbar), color = "blue", size = 2.5) +
    
    # Posterior means with credible intervals
    geom_point(aes(y = m_mean), color = "orange", size = 2.5) +
    geom_errorbar(aes(ymin = m_lower, ymax = m_upper), width = 0.1, color = "orange") +
    
    # Global mean line
    geom_hline(yintercept = post_u, linetype = "dotted", color = "red") +
    
    labs(
      x = "Group",
      y = "Mean",
      title = "Shrinkage of Group Means with 95% Credible Intervals",
      subtitle = "Blue: observed mean, Orange: posterior mean ± 95% CI"
    ) +
    theme_minimal()
  return(g)
}

shrinkage_plot(x_ca,group_ca,mcmc_ca)
shrinkage_plot(x_mo,group_mo,mcmc_mo)










