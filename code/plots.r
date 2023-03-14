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


# Plot perfect reliability

dat2 <- readRDS("data/tau2_increasing_rel.RDS")

dat2 <- lapply(dat2, function(x) as.data.frame(as.list(x)))

dat_tau2 <- data.table::rbindlist(dat2, idcol = "condition")

dat_tau2$true_tau2 <- gsub("&.*", "", dat_tau2$condition)
dat_tau2$true_tau2 <- as.numeric(gsub("tau2 = ", "", dat_tau2$true_tau2))

dat_tau2$reliability <- gsub(".* & reliability = ", "", dat_tau2$condition)
dat_tau2 <- dat_tau2[dat_tau2$reliability == 1,]

dat3 <- readRDS("data/tau2_increasing_fisherz.RDS")

dat3 <- lapply(dat3, function(x) as.data.frame(as.list(x)))

dat2_tau2 <- data.table::rbindlist(dat3, idcol = "condition")

dat2_tau2$true_tau2 <- gsub("&.*", "", dat2_tau2$condition)
dat2_tau2$true_tau2 <- as.numeric(gsub("tau2 = ", "", dat2_tau2$true_tau2))

dat2_tau2$reliability <- gsub(".* & reliability = ", "", dat2_tau2$condition)

dat4 <- readRDS("data/tau2_increasing_r_sampling_var_rho.RDS")

dat4 <- lapply(dat4, function(x) as.data.frame(as.list(x)))

dat3_tau2 <- data.table::rbindlist(dat4, idcol = "condition")

dat3_tau2$true_tau2 <- gsub("&.*", "", dat3_tau2$condition)
dat3_tau2$true_tau2 <- as.numeric(gsub("tau2 = ", "", dat3_tau2$true_tau2))

dat3_tau2$reliability <- gsub(".* & reliability = ", "", dat3_tau2$condition)

dat5 <- readRDS("data/tau2_increasing_fisherz_sampling_var_rho.RDS")

dat5 <- lapply(dat5, function(x) as.data.frame(as.list(x)))

dat4_tau2 <- data.table::rbindlist(dat5, idcol = "condition")

dat4_tau2$true_tau2 <- gsub("&.*", "", dat4_tau2$condition)
dat4_tau2$true_tau2 <- as.numeric(gsub("tau2 = ", "", dat4_tau2$true_tau2))

dat4_tau2$reliability <- gsub(".* & reliability = ", "", dat4_tau2$condition)

# Combining
dat_tau2 <- rbind(dat_tau2, dat2_tau2, dat3_tau2, dat4_tau2)
dat_tau2$condition <- rep(c("r & estimated variance",
                            "z & estimated variance w/ rho",
                            "r & true sampling variance",
                            "z & true sampling variance" ), each = 10)

library(ggplot2)
ggplot(dat_tau2, aes(x = true_tau2, y = tau2_hat)) +
geom_line(aes(color = condition)) +
#scale_linetype_manual(values = 10:1) +
geom_abline() +
scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
theme_bw()

#ggsave("figures/tau2_perfect_reliability.png", width = 8.62, height = 9.93)



#****************************************
# Plots 1 - 4
#****************************************
library(ggplot2)

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

#I need to add zero reliability variance to this plot!

# Plot 1 troubleshooting
library(ggplot2)

dat <- readRDS("data/reliability_overestimate_original.RDS")
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

dat$reliability_min <- gsub(".*& ", "", dat$mu)
dat$reliability_min <- gsub(".*= ", "", dat$reliability_min)
dat$mu <- gsub("&.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))

temp <- split(dat, dat$reliability_min)
temp1 <- temp[[1]][1:11,]
temp2 <- temp[[2]][12:22,]
temp3 <- temp[[3]][23:33,]
temp4 <- temp[[4]][34:44,]
temp5 <- temp[[5]][45:55,]
dat <- rbind(temp1, temp2, temp3, temp4, temp5)

ggplot(data = dat, aes(x = mu, y = sqrt(tau2_hat)))  +
geom_line(aes(x = mu, y = sqrt(tau2_hat), color = reliability_min,
group = reliability_min))

#ggsave("figures/plot1.png", width = 8.62, height = 9.93)


# Plot 2
dat <- readRDS("data/reliability_underestimate.RDS")
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

dat$reliability <- gsub(".*& ", "", dat$mu)
dat$reliability <- as.factor(gsub(".*= ", "", dat$reliability))
dat$mu <- gsub("&.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))

line_labels <- dat[dat$mu == 0.6,]

ggplot(dat,aes(x = mu, y = tau2_hat)) +
geom_line(aes(linetype = reliability), show.legend = FALSE) +
#scale_linetype_manual(values = 10:1) +
geom_text(data = line_labels, aes(x = mu, y = tau2_hat, group = reliability), nudge_x = 0.025,
label = line_labels$reliability) +
#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = 0.078, linetype = "dashed") +
scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
theme_bw()

ggsave("figures/plot2.png", width = 8.62, height = 9.93)

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

dat$reliability_sd <- gsub(".*& ", "", dat$mu)
dat$reliability_sd <- gsub(".*= ", "", dat$reliability_sd)
dat$mu <- gsub("&.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))

line_labels <- dat[dat$mu == 0.6,]

ggplot(dat,aes(x = mu, y = tau2_hat)) +
geom_line(aes(linetype = reliability_sd), show.legend = FALSE) +
#scale_linetype_manual(values = 10:1) +
geom_text(data = line_labels, aes(x = mu, y = tau2_hat, group = reliability_sd), nudge_x = 0.025,
label = line_labels$reliability_sd) +
#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = 0.00078, linetype = "dashed") +
scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
theme_bw()

ggsave("figures/plot4.png", width = 8.62, height = 9.93)
