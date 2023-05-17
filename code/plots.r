# Project: Unexplained heterogeneity
# Code purpose: Creating plots from simulated data
# Author: Anton Olsson-Collentine (anton@olssoncollentine.com)


#****************************************
# functions
#****************************************


source("./code/functions.r") # for function compute_var_truncated
#****************************************
# plotting
#****************************************

library(ggplot2)

#****************************************
# General plot prep
#****************************************



dat_r <- readRDS("data/means_r_over_vs_underestimate.RDS")
dat_z <- readRDS("data/means_z_over_vs_underestimate.RDS")
dat_hs <- readRDS("data/means_r_HS_over_vs_underestimate.RDS")

dat <- c(dat_r, dat_z, dat_hs)
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

#create condition ids
dat <- create_condition_cols(dat)

# Add real truncated values
es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.1, 0.15, 0.2)

tau_trunc <- lapply(tau, compute_var_truncated, mu = es)
tau_trunc <- unlist(tau_trunc)
trunc <- data.frame(mu = es,
                    nominal_tau = rep(tau, each = length(es)),
                    true_tau = tau_trunc)


#Prepping plot
dat$true_tau2 <- as.numeric(dat$true_tau2)
dat$mu <- as.numeric(dat$mu)

dat$tau_hat <- sqrt(dat$tau2_hat)
dat$nominal_tau <- sqrt(dat$true_tau2)

dat$mean_rel <- factor(dat$mean_rel)
dat$k <- factor(dat$k, levels = c("5", "20", "40", "200"))
dat$N <- factor(dat$N, levels = c("50", "100", "150", "200"))


#****************************************
# Main plot 1 (N = 150, k = 20, r, r_z, and HS)
#****************************************
dat_20_150 <- dat[dat$k == 20 & dat$N == 150,]

# For fisher's z there is no truncation,
# create dataframe with all true heterogeneity levels
trunc$effect_type  <- "Pearson's r"
hs_true <- trunc
hs_true$effect_type <- "Hunter & Schmidt" #effect type is r, but named this way for plotting

z_true <- trunc
z_true$effect_type <- "Fisher's z"
dat_z <- dat_20_150[dat_20_150$effect_type == "r_z"]
z_tau <- rep(unique(dat_z$nominal_tau), each = 7) #4 tau-values across 4*7 conditions
z_true$true_tau <- z_tau
#NB! the z_true object has mu expressed in Pearson's r, not fisher z (this is useful for plotting later)

true_tau_line <- rbind(trunc, z_true, hs_true)
true_tau_line$effect_type <- factor(true_tau_line$effect_type,
                            levels = c("Fisher's z",
                                        "Pearson's r",
                                        "Hunter & Schmidt"))



# truncate the fisher z nominal tau-values to match the pearson r ones
#this is purely for plotting purposes. The true value can be seen from the black line in plot.
dat_20_150$nominal_tau <- ifelse(dat_20_150$effect_type == "r_z",
                              trunc(dat_20_150$nominal_tau *100) / 100,
                              dat_20_150$nominal_tau)


# What I've done in the above to lines of code is to
# 1) put the 'truncated' true_tau for r_z to fisher z nominal value
# 2) put the 'true tau' in the main data object to the nominal r value
# This helps with the facet_grid to make the plot look nice

# More plot polish
dat_20_150$effect_type <- ifelse(dat_20_150$effect_type == "r_z",
                                "Fisher's z",
                                "Pearson's r")
dat_20_150$effect_type <- ifelse(dat_20_150$effect_type == "Pearson's r" &
                                 dat_20_150$method == "HS",
                                "Hunter & Schmidt",
                                dat_20_150$effect_type)

dat_20_150$effect_type <- factor(dat_20_150$effect_type,
                            levels = c("Fisher's z",
                                        "Pearson's r",
                                        "Hunter & Schmidt"))

dat_20_150$mu <- ifelse(dat_20_150$effect_type == "Fisher's z",
                        tanh(dat_20_150$mu),
                        dat_20_150$mu)

ggplot(dat_20_150, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = mean_rel), show.legend = TRUE) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE)) +
geom_line(aes(x = mu, y = true_tau), data = true_tau_line, linetype = 1) +
expand_limits(y = 0) +
facet_grid(effect_type~nominal_tau, scales = "free", switch = "y") +
theme_bw()

ggsave("figures/z-r-hs-plot.png", width = 8.62, height = 9.93)

#****************************************
# Main plot 2 (r, N = 150, k variable, tau low)
#****************************************
dat_r_low <- readRDS("data/means_r_borderline_estimate_tau0.02-0.08.RDS")
dat_r <- readRDS("data/means_r_over_vs_underestimate.RDS")
dat_r_zero <- dat_r[grep("true_tau2 = 0;", names(dat_r))]

dat <- c(dat_r_zero, dat_r_low)
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

#create condition ids
dat <- create_condition_cols(dat)

#Prepping plot
dat$true_tau2 <- as.numeric(dat$true_tau2)
dat$mu <- as.numeric(dat$mu)

dat$tau_hat <- sqrt(dat$tau2_hat)
dat$true_tau <- sqrt(dat$true_tau2)

dat$mean_rel <- factor(dat$mean_rel)
dat$k <- factor(dat$k, levels = c("5", "20", "40", "200"))
dat$N <- factor(dat$N, levels = c("50", "100", "150", "200"))

#****************************************
# Truncated values

es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.02, 0.04, 0.06, 0.08)

tau_trunc <- lapply(tau, compute_var_truncated, mu = es)
tau_trunc <- unlist(tau_trunc)
trunc <- data.frame(mu = es,
                    nominal_tau = rep(tau, each = length(es)),
                    true_tau = tau_trunc)


trunc$true_tau <- trunc$nominal_tau
trunc$tau_hat <- trunc$tau

#****************************************

dat_n_150 <- dat[dat$N == 150 & dat$true_tau > 0,]
trunc <- trunc[trunc$true_tau > 0,]

dat_n_150 <- dat[dat$N == 150,]

#Plot with fixed N and varying k, and small tau-values
ggplot(dat_n_150, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = mean_reliability), show.legend = TRUE) +
geom_line(data = trunc, linetype = 1) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE)) +
expand_limits(y = 0) +
facet_grid(k~true_tau, scales = "free") +
theme_bw()

ggsave("figures/r_tau_0.02-0.08.png", width = 8.62, height = 9.93)


#****************************************
# r-plots for supplement (variable N, variable k)
#****************************************
# Plots with many different N and K
# Either for Pearson's r or for Fisher's z
# Tau and mu are equal in both cases

#1) Create all Pearson's r plots
dat_n_150 <- dat[dat$N == 150 & dat$effect_type == "r",]

#Plot with fixed and and varying k
ggplot(dat_n_150, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = mean_reliability), show.legend = TRUE) +
geom_line(data = trunc, linetype = 1) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE)) +
expand_limits(y = 0) +
facet_grid(k~true_tau, scales = "free") +
theme_bw()

ggsave("figures/N_150_r_0.6-0.9_tau_sd_0.15.png", width = 8.62, height = 9.93)

dat_k_20 <- dat[dat$k == 20 & dat$effect_type == "r",]

#Plot with fixed k and and varying N
ggplot(dat_k_20, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = mean_reliability), show.legend = TRUE) +
geom_line(data = trunc, linetype = 1) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE)) +
expand_limits(y = 0) +
facet_grid(N~true_tau, scales = "free") +
theme_bw()

ggsave("figures/k_20_r_0.6-0.9_tau_sd_0.15.png", width = 8.62, height = 9.93)


#****************************************
# z-plots for supplement (variable N, variable k)
#****************************************

dat_n_150 <- dat[dat$N == 150 & dat$effect_type == "r_z",]

#Plot with fixed and and varying k
ggplot(dat_n_150, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = mean_reliability), show.legend = TRUE) +
geom_line(data = trunc, linetype = 1) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE)) +
expand_limits(y = 0) +
facet_grid(k~true_tau, scales = "free") +
theme_bw()

ggsave("figures/z_N_150_0.6-0.9_tau_sd_0.15.png", width = 8.62, height = 9.93)


dat_k_20 <- dat[dat$k == 20 & dat$effect_type == "r_z",]


#Plot with fixed k and and varying N
ggplot(dat_k_20, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = mean_reliability), show.legend = TRUE) +
geom_line(data = trunc, linetype = 1) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE)) +
expand_limits(y = 0) +
facet_grid(N~true_tau, scales = "free") +
theme_bw()

ggsave("figures/z_k_20_0.6-0.9_tau_sd_0.15.png", width = 8.62, height = 9.93)
