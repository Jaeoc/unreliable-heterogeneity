# Project: Unexplained heterogeneity
# Code purpose: Creating plots from simulated data
# Author: Anton Olsson-Collentine (anton@olssoncollentine.com)


#****************************************
# Simulating data for plots
#****************************************
# see supercomputer_simulator.r


#****************************************
# plotting
#****************************************

library(ggplot2)

#****************************************
# Plots 1 - 4
#****************************************

dat <- readRDS("data/reliability_overestimate.RDS")
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

dat$reliability_sd <- gsub(".*& ", "", dat$mu)
dat$reliability_sd <- gsub(".*= ", "", dat$reliability_sd)
dat$mu <- gsub("&.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))


ggplot(data = dat) +
geom_line(aes(x = mu, y = sqrt(tau2_hat), color = reliability_sd))

#ggsave("figures/plot1.png", width = 8.62, height = 9.93)


# Plot 2 original

# line_labels <- dat[dat$r == "0.6",]

# ggplot(dat,aes(x = r, y = percent_tau2)) +
# geom_line(aes(linetype = reliability), show.legend = FALSE) +
# scale_linetype_manual(values = 8:1) +
# geom_label(data = line_labels, aes(x = r, y = percent_tau2, group = reliability), nudge_x = 0.025,
# label = line_labels$reliability) +
# #annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
# geom_hline(yintercept = 100, linetype = "dashed") +
# theme_bw()

# ggsave("figures/plot2.png", width = 8.62, height = 9.93)

## Plot 2 v2
# That is, with perfect reliability and tau2-values of
# c(0.002, 0.0035, 0.0055, 0.0085, 0.012, 0.0185, 0.031, 0.069)
# (i.e, 20 - 90% I2)

dat <- readRDS("data/k_20_over_vs_underestimate.RDS")
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

#create condition ids
conditions <- unlist(strsplit(dat$mu, split = ";"))

dat$true_tau2 <- grep("true_tau2 =", conditions, value = TRUE)
dat$true_tau2 <- as.numeric(gsub("true_tau2 = ", "", dat$true_tau2))
dat$mean_reliability <- grep("mean_rel", conditions, value = TRUE)
dat$mean_reliability <- gsub("mean_rel = ", "", dat$mean_reliability)
dat$mu <- gsub(";.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))

# Only extract reliability = 1
dat <- dat[dat$mean_reliability == 1,]

# Add real truncated values
trunc <- readRDS("data/truncated_tau2.RDS")
trunc$true_tau2 <- trunc$nominal_tau2

ggplot(dat,aes(x = mu, y = tau2_hat)) +
geom_line() +
geom_line(aes(x = mu, y = tau2), linetype = "dashed", data = trunc) +
facet_wrap(~true_tau2, scales = "free")

ggsave("figures/reliability_1_tau2.png", width = 8.62, height = 9.93)

# Plot2 with tau instead of tau2
dat$tau_hat <- sqrt(dat$tau2_hat)
dat$true_tau <- round(sqrt(dat$true_tau2),3)
trunc$true_tau <- round(sqrt(trunc$true_tau2), 3)
trunc$tau <- sqrt(trunc$tau2)

ggplot(dat,aes(x = mu, y = tau_hat)) +
geom_line() +
geom_line(aes(x = mu, y = tau), linetype = "dashed", data = trunc) +
expand_limits(y = 0) +
facet_wrap(~true_tau, scales = "free")+
theme_bw()


ggsave("figures/reliability_1_tau.png", width = 8.62, height = 9.93)

# Plot 3


dat2 <- readRDS("data/increasing_true_tau2.RDS")

dat2 <- lapply(dat2, function(x) as.data.frame(as.list(x)))

dat_tau2 <- data.table::rbindlist(dat2, idcol = "condition")

dat_tau2$true_tau2 <- gsub("&.*", "", dat_tau2$condition)
dat_tau2$true_tau2 <- as.numeric(gsub("tau2 = ", "", dat_tau2$true_tau2))

dat_tau2$reliability <- gsub(".* & reliability = ", "", dat_tau2$condition)

tau2_labels <- dat_tau2[dat_tau2$true_tau2 == 0.27,]

ggplot(dat_tau2, aes(x = true_tau2, y = tau2_hat)) +
geom_line(aes(linetype = reliability), show.legend = FALSE) +
#scale_linetype_manual(values = 10:1) +
geom_label(data = tau2_labels, aes(x = true_tau2, y = tau2_hat, group = reliability), nudge_x = 0.0025,
label = tau2_labels$reliability) +
annotate("text", x = 0.24, y = 0.3, label = "Reliability", hjust = 0) +
geom_abline() +
scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
theme_bw()

ggsave("figures/plot3.png", width = 8.62, height = 9.93)


# Plot 4

dat <- readRDS("data/k_5_k20_tau2_over_vs_underestimate.RDS")
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

#create condition ids
conditions <- unlist(strsplit(dat$mu, split = ";"))

dat$reliability_sd <- grep("reliability_sd", conditions, value = TRUE)
dat$reliability_sd <- gsub("reliability_sd = ", "", dat$reliability_sd)
dat$true_tau2 <- grep("true_tau2 =", conditions, value = TRUE)
dat$true_tau2 <- as.numeric(gsub("true_tau2 = ", "", dat$true_tau2))
dat$mean_reliability <- grep("mean_rel", conditions, value = TRUE)
dat$mean_reliability <- gsub("mean_rel = ", "", dat$mean_reliability)
dat$mu <- gsub(";.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))

# remove reliability = 1
dat <- dat[dat$mean_reliability < 1,]

line_labels <- dat[dat$mu == "0.6",]

# Add real truncated values
trunc <- readRDS("data/k_20_truncated_tau2.RDS")
trunc$true_tau2 <- trunc$nominal_tau2
trunc$tau2_hat <- trunc$tau2

# Supplement: should do one facet plot for each mean reliability
ggplot(dat,aes(x = mu, y = tau2_hat)) +
geom_line(aes(linetype = mean_reliability, color = reliability_sd), show.legend = FALSE) +
geom_line(data = trunc, color = "black", linetype = 4) +
#scale_linetype_manual(values = 10:1) +
geom_text(data = line_labels, aes(x = mu, y = tau2_hat, group = reliability_sd), nudge_x = 0.025,
label = line_labels$reliability_sd) +
expand_limits(y = 0) +
facet_wrap(~true_tau2, scales = "free") +
theme_bw()

ggsave("figures/k_20_R_0.6-0.9_tau2.png", width = 8.62, height = 9.93)

## Plot 4 but with only SD = 0.15
dat <- dat[dat$reliability_sd == "0.15",]
line_labels <- dat[dat$mu == "0",]

ggplot(dat,aes(x = mu, y = tau2_hat)) +
geom_line(aes(linetype = mean_reliability), show.legend = FALSE) +
geom_line(data = trunc, color = "black", linetype = 1) +
scale_linetype_manual(values = 2:5) +
geom_text(data = line_labels, aes(x = mu, y = tau2_hat, group = mean_reliability), nudge_x = 0.025,
label = line_labels$mean_reliability) +
expand_limits(y = 0) +
facet_wrap(~true_tau2, scales = "free") +
theme_bw()

ggsave("figures/k_20_R_0.6-0.9_tau2_sd_0.15.png", width = 8.62, height = 9.93)


# With tau instead of tau2
dat$tau_hat <- sqrt(dat$tau2_hat)
dat$true_tau <- round(sqrt(dat$true_tau), 3)
trunc$tau_hat <- sqrt(trunc$tau2_hat)
trunc$true_tau <- round(sqrt(trunc$true_tau2), 3)

line_labels$tau_hat <- sqrt(line_labels$tau2_hat)
line_labels$true_tau <- round(sqrt(line_labels$true_tau2), 3)



ggplot(dat,aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = mean_reliability), show.legend = TRUE) +
geom_line(data = trunc, linetype = 1) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE)) +
# geom_text(data = line_labels, aes(x = mu, y = tau_hat), nudge_x = -0.05,
# label = line_labels$mean_reliability, size = 3) +
expand_limits(y = 0) +
facet_wrap(~true_tau, scales = "free") +
theme_bw()


ggsave("figures/k_5_k20tau_R_0.6-0.9_tau_sd_0.15.png", width = 8.62, height = 9.93)




# Robbie plot
# That is: Plot2 with perfect reliability and the following conditions
#a1) tau2 = 0.078, Pearson's r
#a2) as above but using rho when estimating the sampling variance
#a3) as a1 but transforming to Fisher's z before meta-analysis
#b1) generate data as fisher's z
#b2) as above, but estimate heterogeneity with DL

dat <- readRDS("data/robbie.RDS")
names(dat) <- c("r", "r with rho", "r to z", "z REML", "z DL")

dat <- data.table::rbindlist(dat, idcol = "condition")
dat$mu <- rep(seq(from = 0, to = 0.6, by = 0.1), 5)

ggplot(dat,aes(x = mu, y = prop)) +
geom_line(aes(color = condition), show.legend = TRUE) +
#scale_linetype_manual(values = 8:1) +
#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = 1, linetype = "dashed") +
ylab("prop of true tau2") +
theme_bw()

dat2 <- dat[dat$condition != "r to z",]

ggplot(dat2,aes(x = mu, y = prop)) +
geom_line(aes(color = condition), show.legend = TRUE) +
#scale_linetype_manual(values = 8:1) +
#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = 1, linetype = "dashed") +
ylab("prop of true tau2") +
theme_bw()

#absolute values
dat3 <- dat[dat$condition == "r",]
ggplot(dat3,aes(x = mu, y = diff)) +
geom_line(aes(color = condition), show.legend = TRUE) +
#scale_linetype_manual(values = 8:1) +
#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = sqrt(0.078), linetype = "dashed") +
ylab("prop of true tau2") +
theme_bw()


#****************************************
# supplemental plots + fisher's z
#****************************************

# Plots with many different N and K
# Either for Pearson's r or for Fisher's z
# Tau and mu are equal in both cases

#1) Create all Pearson's r plots


dat_r <- readRDS("data/means_r_over_vs_underestimate.RDS")
dat_z <- readRDS("data/means_z_over_vs_underestimate.RDS")
dat <- c(dat_r, dat_z)
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

#create condition ids
grep_strip <- function(string, vec){
    #Uses a common string like "k = " to extract "k = 50" or "k = 200"
    #Then removes the "k = " part and returns as a numeric value
    a <- grep(string, vec, value = TRUE)
    b <- gsub(string, "", a)
}

conditions <- unlist(strsplit(dat$mu, split = ";"))

dat$k <- grep_strip("k = ", conditions)
dat$N <- grep_strip("N = ", conditions)
dat$reliability_sd <- grep_strip("reliability_sd = ", conditions)
dat$true_tau2 <- as.numeric(grep_strip("true_tau2 =", conditions)) #numeric for plotting later
dat$mean_reliability <- grep_strip("mean_rel = ", conditions)
dat$mu <- as.numeric(grep_strip("mu = ", conditions))
dat$effect_type <- grep_strip("effect_type = ", conditions)

# Add real truncated values
trunc <- readRDS("data/truncated_tau.RDS")
trunc$true_tau <- trunc$nominal_tau
trunc$tau_hat <- trunc$tau

#Prepping plot
dat$tau_hat <- sqrt(dat$tau2_hat)
dat$true_tau <- sqrt(dat$true_tau2)
dat$mean_reliability <- factor(dat$mean_reliability)
dat$k <- factor(dat$k, levels = c("5", "20", "40", "200"))
dat$N <- factor(dat$N, levels = c("50", "100", "150", "200"))

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



#2) Fisher's z


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
