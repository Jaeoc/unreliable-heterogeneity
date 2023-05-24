# Project: unexplained heterogeneity
# Script purpose: functions for simulations
# code: Anton Olsson-Collentine

#****************************************
# functions and packages
#****************************************

source("./code/functions.r") # for create_condition_cols function


#****************************************
# [1] Result section
#****************************************
dat_r <- readRDS("data/means_r_tau_0-0.2.RDS")

dat_r$tau_hat <- sqrt(dat_r$tau2_hat)
dat_r$true_tau <- sqrt(as.numeric(dat_r$true_tau2))

dat <- dat_r[dat_r$k == "20" & dat_r$N == "150",]

#*[1.1] bias when tau =0, R = 0.8 ----
#****************************************

dat_1.1 <- dat[dat$true_tau == 0 & dat$mean_rel == 0.6 ,]
#mu = 0.2
mu_0.2 <- dat_1.1$tau_hat[dat_1.1$mu == 0.2] - dat_1.1$true_tau[dat_1.1$mu == 0.2]
mu_0.2 <- round(mu_0.2, 3)

#mu = 0.4
mu_0.4 <- dat_1.1$tau_hat[dat_1.1$mu == 0.4] - dat_1.1$true_tau[dat_1.1$mu == 0.4]
mu_0.4 <- round(mu_0.4, 3)

bias_tau_zero <- data.frame(mu_0.2, mu_0.4)

#*[1.2] bias when tau increasing, mu = 0.2, R = 0.8 ----
#****************************************
dat_1.2 <- dat[dat$mu == 0.2 & dat$mean_rel == 0.8,]
#tau = 0.1
tau_hat <- dat_1.2$tau_hat[dat_1.2$true_tau == 0.1]
true_tau <- dat_1.2$true_tau[dat_1.2$true_tau == 0.1]
bias_0.1 <-  tau_hat - true_tau
percent_0.1 <- 100*(tau_hat / true_tau)

res_0.1 <- paste0(round(bias_0.1, 3), " (", round(percent_0.1, 0), "%)")

#tau = 0.15
tau_hat <- dat_1.2$tau_hat[dat_1.2$true_tau == 0.15]
true_tau <- dat_1.2$true_tau[dat_1.2$true_tau == 0.15]
bias_0.15 <-  tau_hat - true_tau
percent_0.15 <- 100*(tau_hat / true_tau)

res_0.15 <- paste0(round(bias_0.15, 3), " (", round(percent_0.15, 0), "%)")

#tau = 0.2
tau_hat <- dat_1.2$tau_hat[dat_1.2$true_tau == 0.2]
true_tau <- dat_1.2$true_tau[dat_1.2$true_tau == 0.2]
bias_0.2 <-  tau_hat - true_tau
percent_0.2 <- 100*(tau_hat / true_tau)

res_0.2 <- paste0(round(bias_0.2, 3), " (", round(percent_0.2, 0), "%)")

bias_increasing_tau <- list(tau_0.1 = res_0.1,
                            tau_0.15 = res_0.15,
                            tau_0.2 = res_0.2)

saveRDS(list(
    bias_tau_zero = bias_tau_zero,
    bias_increasing_tau = bias_increasing_tau),
    file = "./manuscript/in-text-values.RDS"
)
