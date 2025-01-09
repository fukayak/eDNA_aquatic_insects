setwd("<specify your directory>")
load("data_raw.RData")

if (!require("occumb")) install.packages("occumb")
library(occumb)
source("functions.R")



### Build the data object
data <- occumbData(
  y = data_raw$y,
  spec_cov = list(order = data_raw$order,
                  mismatch = data_raw$mismatch),
  repl_cov = list(vol = data_raw$vol)
)
#summary(data)



### Fit the model
fit <- occumb(formula_phi_shared = ~ order + mismatch,
              formula_theta_shared = ~ order + vol,
              n.iter = 250000,
              n.burnin = 50000,
              n.thin = 200,
              data = data,
              parallel = TRUE)
# The above settings are intended to ensure sufficient convergence, and
# it will take several days to execute.
# Preliminary analysis could be run on shorter chains.

#save(fit, file = "fit.RData")
#summary(fit)
#plot(fit)



### Assess goodness-of-fit
gof_fit <- gof(fit)
gof_fit



### Get posterior results
newdata <- occumbData(
  y = data_raw$y[, 1, , drop = FALSE], # selects only one site
  spec_cov = list(order = data_raw$order,
                  mismatch = data_raw$mismatch),
  repl_cov = list(vol = data_raw$vol[1, , drop = FALSE])
)

## phi, theta, and psi estimates
predict_phi   <- predict(fit, newdata, "phi")
predict_theta <- predict(fit, newdata, "theta")[, , , 1] # selects 1000 ml
predict_psi   <- predict(fit, newdata, "psi")

# Figure 1
plot_fig1(
  phi   = predict_phi,
  theta = predict_theta,
  psi   = predict_psi
)

## Correlation coefficients of phi, theta, and psi
post_sample_phi   <- get_post_samples(fit, "phi")
post_sample_theta <- get_post_samples(fit, "theta")[, , 1, 1] # selects only one site and 1000 ml
post_sample_psi   <- get_post_samples(fit, "psi")

cor_mat <- get_cor_mat(
  phi   = post_sample_phi,
  theta = post_sample_theta,
  psi   = post_sample_psi
)
cor_mat # Table 1

## Coefficients
coefficients <- rbind(
  get_post_summary(fit, "beta_shared")[, c("50%", "2.5%", "97.5%")],
  get_post_summary(fit, "alpha_shared")[, c("50%", "2.5%", "97.5%")]
)
coefficients # Table 2



### Compare different study settings
## Unconditional profiles
uncond_1000 <- eval_util_L(
  data.frame(
    K = rep(c(1, 2, 3, 4, 8, 16), each = 7),
    N = rep(c(0.5E4, 1E4, 2E4, 4E4, 8E4, 16E4, 32E4), times = 6)
  ),
  fit,
  theta = predict(fit, newdata, "theta", type = "samples")[, , , 1], # selects 1000 ml
  cores = 8
)
uncond_800 <- eval_util_L(
  data.frame(
    K = rep(c(1, 2, 3, 4, 8, 16), each = 7),
    N = rep(c(0.5E4, 1E4, 2E4, 4E4, 8E4, 16E4, 32E4), times = 6)
  ),
  fit,
  theta = predict(fit, newdata, "theta", type = "samples")[, , , 2], # selects 800 ml
  cores = 8
)
uncond_600 <- eval_util_L(
  data.frame(
    K = rep(c(1, 2, 3, 4, 8, 16), each = 7),
    N = rep(c(0.5E4, 1E4, 2E4, 4E4, 8E4, 16E4, 32E4), times = 6)
  ),
  fit,
  theta = predict(fit, newdata, "theta", type = "samples")[, , , 3], # selects 600 ml
  cores = 8
)
uncond_400 <- eval_util_L(
  data.frame(
    K = rep(c(1, 2, 3, 4, 8, 16), each = 7),
    N = rep(c(0.5E4, 1E4, 2E4, 4E4, 8E4, 16E4, 32E4), times = 6)
  ),
  fit,
  theta = predict(fit, newdata, "theta", type = "samples")[, , , 4], # selects 400 ml
  cores = 8
)

# Figure 2
plot_fig2(
  result1000 = uncond_1000,
  result800  = uncond_800,
  result600  = uncond_600,
  result400  = uncond_400
)

## Conditional profiles
budget <- 540000; lambda1 <- 0.025; lambda2 <- 5750
settings <- list_cond_L(budget, lambda1, lambda2, fit)

cond_1000 <- eval_util_L(
  settings,
  fit,
  theta = predict(fit, newdata, "theta", type = "samples")[, , , 1], # selects 1000 ml
  cores = 8
)
cond_800 <- eval_util_L(
  settings,
  fit,
  theta = predict(fit, newdata, "theta", type = "samples")[, , , 2], # selects 800 ml
  cores = 8
)
cond_600 <- eval_util_L(
  settings,
  fit,
  theta = predict(fit, newdata, "theta", type = "samples")[, , , 3], # selects 600 ml
  cores = 8
)
cond_400 <- eval_util_L(
  settings,
  fit,
  theta = predict(fit, newdata, "theta", type = "samples")[, , , 4], # selects 400 ml
  cores = 8
)

# Figure 3
plot_fig3(
  result1000 = cond_1000,
  result800  = cond_800,
  result600  = cond_600,
  result400  = cond_400
)

