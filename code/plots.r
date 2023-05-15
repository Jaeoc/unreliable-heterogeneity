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



#****************************************
# supplemental plots + fisher's z
#****************************************
dat_20_150 <- dat[dat$k == 20 & dat$N == 150,]

# For fisher's z there is no truncation,
# hack the trunc dataframe for the plot
trunc$effect_type  <- "Pearson's r"
trunc <- rbind(trunc, trunc)
trunc[29:56, "effect_type"] <- "Fisher's z"

dat_z <- dat_20_150[dat_20_150$effect_type == "r_z"]
trunc$tau_hat[29:56] <- rep(unique(dat_z$true_tau), each = 7)

dat_20_150$true_tau <- ifelse(dat_20_150$effect_type == "r_z",
                              trunc(dat_20_150$true_tau *100) / 100,
                              dat_20_150$true_tau)


# What I've done in the above to lines of code is to
# 1) put the 'truncated' true_tau for r_z to fisher z nominal value
# 2) put the 'true tau' in the main data object to the nominal r value
# This helps with the facet_grid to make the plot look nice
# But changes the interpretation of the black line for fisher z
# Does it? It is still the true heterogeneity

# More plot prep
dat_20_150$effect_type <- ifelse(dat_20_150$effect_type == "r_z",
                                "Fisher's z",
                                "Pearson's r")


ggplot(dat_20_150, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = mean_reliability), show.legend = TRUE) +
geom_line(data = trunc, linetype = 1) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE)) +
expand_limits(y = 0) +
facet_grid(effect_type~true_tau, scales = "free", switch = "y") +
theme_bw()
