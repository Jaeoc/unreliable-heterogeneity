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


# Plot 2

line_labels <- dat[dat$r == "0.6",]

ggplot(dat,aes(x = r, y = percent_tau2)) +
geom_line(aes(linetype = reliability), show.legend = FALSE) +
scale_linetype_manual(values = 8:1) +
geom_label(data = line_labels, aes(x = r, y = percent_tau2, group = reliability), nudge_x = 0.025,
label = line_labels$reliability) +
#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = 100, linetype = "dashed") +
theme_bw()

ggsave("figures/plot2.png", width = 8.62, height = 9.93)

## Plot 2 v2
# That is, with perfect reliability and tau2-values of
# c(0.002, 0.0035, 0.0055, 0.0085, 0.012, 0.0185, 0.031, 0.069)
# (i.e, 20 - 90% I2)

dat <- readRDS("data/reliability_1.RDS")
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

#create condition ids
conditions <- unlist(strsplit(dat$mu, split = ";"))

dat$true_tau2 <- grep("true_tau2 =", conditions, value = TRUE)
dat$true_tau2 <- as.numeric(gsub("true_tau2 = ", "", dat$true_tau2))
dat$mu <- gsub(";.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))

# Add real truncated values
trunc <- readRDS("data/truncated_tau2.RDS")
trunc$true_tau2 <- trunc$nominal_tau2

ggplot(dat,aes(x = mu, y = tau2_hat)) +
geom_line() +
geom_line(aes(x = mu, y = tau2), linetype = "dashed", data = trunc) +
facet_wrap(~true_tau2)



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

dat <- readRDS("data/over_vs_underestimate.RDS")
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

dat$prop_tau2 <-dat$tau2_hat /  dat$true_tau2
dat$percent_tau2 <- 100 * dat$prop_tau2


line_labels <- dat[dat$mu == "0.6",]

ggplot(dat,aes(x = mu, y = tau2_hat)) +
geom_line(aes(linetype = mean_reliability, color = reliability_sd), show.legend = FALSE) +
#scale_linetype_manual(values = 10:1) +
geom_text(data = line_labels, aes(x = mu, y = tau2_hat, group = reliability_sd), nudge_x = 0.025,
label = line_labels$reliability_sd) +
facet_wrap(~true_tau2, scales = "free")

#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = 100, linetype = "dashed") +
scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
ylim(c(0, 100)) +
theme_bw()

ggsave("figures/plot4.png", width = 8.62, height = 9.93)


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
2
