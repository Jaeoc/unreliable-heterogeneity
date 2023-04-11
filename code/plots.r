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
dat$reliability <- gsub(".*= ", "", dat$reliability)
dat$mu <- gsub("&.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))

dat$prop_tau2 <-dat$tau2_hat /  0.069  #0.069 is the tau2 used in simulation
dat$percent_tau2 <- 100 * dat$prop_tau2


#convert back to Pearson's r
dat$r <- z_to_r(dat$mu)

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


ggplot(dat,aes(x = mu, y = tau2_hat)) +
geom_line() +
facet_wrap(~true_tau2, scales = "free")

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

dat$reliability_sd <- gsub(".*& ", "", dat$mu)
dat$reliability_sd <- gsub(".*= ", "", dat$reliability_sd)
dat$mu <- gsub("&.*", "", dat$mu)
dat$mu <- as.numeric(gsub("mu = ", "", dat$mu))

dat$prop_tau2 <-dat$tau2_hat /  0.069  #0.00078 is the tau2 used in simulation
dat$percent_tau2 <- 100 * dat$prop_tau2

#convert back to Pearson's r
dat$r <- z_to_r(dat$mu)

line_labels <- dat[dat$r == "0.6",]

ggplot(dat,aes(x = r, y = percent_tau2)) +
geom_line(aes(linetype = reliability_sd), show.legend = FALSE) +
#scale_linetype_manual(values = 10:1) +
geom_text(data = line_labels, aes(x = r, y = percent_tau2, group = reliability_sd), nudge_x = 0.025,
label = line_labels$reliability_sd) +
#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = 100, linetype = "dashed") +
scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
ylim(c(0, 100)) +
theme_bw()

ggsave("figures/plot4.png", width = 8.62, height = 9.93)




# Plot 4 v2
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
