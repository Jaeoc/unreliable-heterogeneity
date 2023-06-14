# Project: unexplained heterogeneity
# Script purpose: functions for simulations
# code: Anton Olsson-Collentine

#****************************************
# functions and packages
#****************************************

source("./code/02_functions.r") # for create_condition_cols function


#****************************************
# [1] Result section
#****************************************
dat <- read.csv("data/means_combined.csv")

dat$tau_hat <- sqrt(dat$tau2_hat)
dat$true_tau <- sqrt(as.numeric(dat$true_tau2))
dat$k <- factor(dat$k, levels = c("5", "20", "40", "200"))
dat$N <- factor(dat$sample_size, levels = c("50", "100", "150", "200"))

dat <- dat[dat$k == "20" &
           dat$N == "150" &
           dat$reliability_sd == 0.15 &
           dat$method == "HV" &
           dat$effect_type == "r",]

#*[1.1] bias when tau =0, R = 0.8 ----
#****************************************

dat_1.1 <- dat[dat$true_tau == 0 & dat$reliability_mean == 0.8 ,]
#mu = 0.2
mu_0.2 <- dat_1.1$tau_hat[dat_1.1$mu == 0.2] - dat_1.1$true_tau[dat_1.1$mu == 0.2]
mu_0.2 <- round(mu_0.2, 3)

#mu = 0.4
mu_0.4 <- dat_1.1$tau_hat[dat_1.1$mu == 0.4] - dat_1.1$true_tau[dat_1.1$mu == 0.4]
mu_0.4 <- round(mu_0.4, 3)

bias_tau_zero <- data.frame(mu_0.2, mu_0.4)

#*[1.2] bias when tau increasing, mu = 0.2, R = 0.8 ----
#****************************************
dat_1.2 <- dat[dat$mu == 0.2 & dat$reliability_mean == 0.8,]
#tau = 0.1
tau_hat <- dat_1.2$tau_hat[dat_1.2$true_tau == 0.1]
true_tau <- dat_1.2$true_tau[dat_1.2$true_tau == 0.1]
bias_0.1 <-  tau_hat - true_tau
prop_0.1 <- 1 - (tau_hat / true_tau)

res_0.1 <- paste0(round(bias_0.1, 3), " (", round(100*prop_0.1, 0), "%)")

#tau = 0.15
tau_hat <- dat_1.2$tau_hat[dat_1.2$true_tau == 0.15]
true_tau <- dat_1.2$true_tau[dat_1.2$true_tau == 0.15]
bias_0.15 <-  tau_hat - true_tau
prop_0.15 <- 1 - (tau_hat / true_tau)

res_0.15 <- paste0(round(bias_0.15, 3), " (", round(100*prop_0.15, 0), "%)")

#tau = 0.2
tau_hat <- dat_1.2$tau_hat[dat_1.2$true_tau == 0.2]
true_tau <- dat_1.2$true_tau[dat_1.2$true_tau == 0.2]
bias_0.2 <-  tau_hat - true_tau
prop_0.2 <- 1 - (tau_hat / true_tau)

res_0.2 <- paste0(round(bias_0.2, 3), " (", round(100*prop_0.2, 0), "%)")

bias_increasing_tau <- list(tau_0.1 = res_0.1,
                            tau_0.15 = res_0.15,
                            tau_0.2 = res_0.2)


#*[1.3] bias when mu = 0.2 and tau => 0.1 ----
#****************************************
dat02 <- dat[dat$mu == 0.2 & dat$true_tau >= 0.1,]
dat02$prop <- dat02$tau_hat / dat02$true_tau

dat02$underestimate <- 1 - dat02$prop

#"heterogeneity can be expected to be underestimated by XX%"
underestimate_min <- data.frame(dat02[which.min(dat02$underestimate),
 c("true_tau", "reliability_mean", "underestimate")])

underestimate_max <- data.frame(dat02[which.max(dat02$underestimate),
 c("true_tau", "reliability_mean", "underestimate")])

underestimate_mu02 <- list(underestimate_min = underestimate_min,
                           underestimate_max = underestimate_max)


## [1.4] sample size (N = 50 vs N = 200, mu = 0.2, tau >0)
#****************************************

dat <- read.csv("data/means_combined.csv")
dat$tau_hat <- sqrt(dat$tau2_hat)
dat$nominal_tau <- sqrt(dat$true_tau2)
dat$k <- factor(dat$k, levels = c("5", "20", "40", "200"))
dat$N <- factor(dat$sample_size, levels = c("50", "100", "150", "200"))

dat_3 <- dat[dat$k == 20 &
             dat$reliability_sd == 0.15 &
             dat$nominal_tau %in% c(0.1, 0.15, 0.2) &
             dat$sample_size %in% c(50, 200) &
             dat$method == "HV" &
             dat$effect_type == "r",]

### [1.4.1] differences when  reliability is 0.8
#****************************************

dat_3b <- dat_3[dat_3$reliability_mean == 0.8 &
              dat_3$mu == 0.2,]

res2 <- split(dat_3b, dat_3b$nominal_tau)
res2 <- lapply(res2, function(x) {
    x[1, "tau_hat"] - x[2, "tau_hat"]
})
res2 <- abs(unlist(res2))

max_bias_n50_0.8 <- round(max(res2), 2)

### [1.4.2] difference reliability is 0.6
#****************************************

dat_3a <- dat_3[dat_3$reliability_mean == 0.6 &
                dat_3$mu == 0.2,]

res1 <- split(dat_3a, dat_3a$nominal_tau)
res1 <- lapply(res1, function(x) {
    x[1, "tau_hat"] - x[2, "tau_hat"]
})
res1 <- abs(unlist(res1))

max_bias_n50_0.6 <- round(max(res1), 2)

### [1.4.2] bias when mu = 0.2, tau = 0.15, N = 50, R =  0.6
#****************************************
dat_142 <- dat_3a[dat_3a$nominal_tau == 0.15 &
                  dat_3a$N == 50,]

net_bias_n50 <-  dat_142$tau_hat - dat_142$nominal_tau
net_bias_n50  <- round(net_bias_n50, 3)

###*[1.4.4] bias when mu = 0.2 and tau => 0.1 ----
#****************************************
dat50 <- dat_3[dat_3$mu == 0.2 &
               dat_3$nominal_tau >= 0.1 &
               dat_3$N == 50,]

dat50$prop <- dat50$tau_hat / dat50$nominal_tau

dat50$underestimate <- 1 - dat50$prop

#"heterogeneity can be expected to be underestimated by XX%"
underestimate_min50 <- data.frame(dat50[which.min(dat50$underestimate),
 c("nominal_tau", "reliability_mean", "underestimate")])

underestimate_max50 <- data.frame(dat50[which.max(dat50$underestimate),
 c("nominal_tau", "reliability_mean", "underestimate")])

underestimate_50 <- list(underestimate_min = underestimate_min50,
                           underestimate_max = underestimate_max50)


#****************************************
# [2] Discussion
#****************************************
k  <- 12

sample_size = 150
mu <- 0
true_tau2 <- 0.17^2

quick_sim <- function(k, sample_size, mu, true_tau2){
    rho <- rnorm_truncated(n = k, mean = mu, sd = sqrt(true_tau2),
                            lower_bound = -1, upper_bound = 1)
    sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)

    # given study rho, draw sample rho
    r_se <- rnorm_truncated(k, mean = rho, sd = sqrt(sampling_var_rho),
                            lower_bound = -1, upper_bound = 1)
    r_var_estimate <- (1-r_se^2)^2 / (sample_size -1)
    fit <- metafor::rma(yi = r_se, vi = r_var_estimate)

    cis <- confint(fit)

    CI_range <- cis$random[2, 3] - cis$random[2, 2]

    CI_range
}

quick_sim(k = k, sample_size = sample_size, mu = mu, true_tau2 = true_tau2)

set.seed(1532)
res <- replicate(1e4, quick_sim(k = k, sample_size = sample_size, mu = mu, true_tau2 = true_tau2))
median(res)



#****************************************
# [3] Output
#****************************************

saveRDS(list(
    bias_tau_zero = bias_tau_zero,
    bias_increasing_tau = bias_increasing_tau,
    underestimate_mu02 = underestimate_mu02,
    bias_n50_0.8 = max_bias_n50_0.8,
    bias_n50_0.6 = max_bias_n50_0.6,
    net_bias_n50 = net_bias_n50,
    underestimate_50 = underestimate_50),
    file = "./manuscript/in-text-values.RDS"
)
