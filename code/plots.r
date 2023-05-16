# Project: Unexplained heterogeneity
# Code purpose: Creating plots from simulated data
# Author: Anton Olsson-Collentine (anton@olssoncollentine.com)


#****************************************
# functions
#****************************************
grep_strip <- function(string, vec){
    #Uses a common string like "k = " to extract "k = 50" or "k = 200"
    #Then removes the "k = " part and returns as character
    a <- grep(string, vec, value = TRUE)
    gsub(string, "", a)
}

create_condition_cols <- function(data){

    #takes a data frame as input with a column named "mu"
    #in which all conditions are separated by a semicolon
    # Then creates new columns named after each condition

    #split all conditions into a character vector
    conditions <- unlist(strsplit(data$mu, split = ";"))
    #select all unique condition factors (e.g., mu, k, N)
    category <- unique(gsub("=.*", "= ", conditions))
    # Remove equal sign to column names (kept above in case some condition name include 'k')
    col_name <- gsub(" = ", "", category)

    for(g in 1:length(category)){

        #This loop is equivalent to the line below, but applied to all conditions
        # data$k <- grep_strip("k = ", conditions)

        data[[col_name[g]]] <- grep_strip(category[g], conditions)

    }

    data #out
}
#****************************************
# plotting
#****************************************

library(ggplot2)

#****************************************
# General plot prep
#****************************************

# Plots with many different N and K
# Either for Pearson's r or for Fisher's z
# Tau and mu are equal in both cases

#1) Create all Pearson's r plots


dat_r <- readRDS("data/means_r_over_vs_underestimate.RDS")
dat_z <- readRDS("data/means_z_over_vs_underestimate.RDS")
dat_hs <- readRDS("data/means_r_HS_over_vs_underestimate.RDS")

dat <- c(dat_r, dat_z, dat_hs)
dat <- lapply(dat, function(x) as.data.frame(as.list(x)))

dat <- data.table::rbindlist(dat, idcol = "mu")

#create condition ids
dat <- create_condition_cols(dat)


# Add real truncated values
trunc <- readRDS("data/truncated_tau.RDS")
trunc$true_tau <- trunc$nominal_tau
trunc$tau_hat <- trunc$tau

#Prepping plot
dat$true_tau2 <- as.numeric(dat$true_tau2)
dat$mu <- as.numeric(dat$mu)

dat$tau_hat <- sqrt(dat$tau2_hat)
dat$true_tau <- sqrt(dat$true_tau2)

dat$mean_rel <- factor(dat$mean_rel)
dat$k <- factor(dat$k, levels = c("5", "20", "40", "200"))
dat$N <- factor(dat$N, levels = c("50", "100", "150", "200"))


#****************************************
# Main plot 1 (N = 150, k = 20, r, r_z, and HS)
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

source("./code/functions.r")
es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.02, 0.04, 0.06, 0.08)

a <- lapply(tau, function(x){
        lapply(es, function(y) {
            res <- rnorm_truncated(n = 1e6, mean = y, sd = x,
            lower_bound = -1, upper_bound = 1)

            data.frame(tau = sd(res), mu = y, nominal_tau = x)
        })
})

b <- lapply(a, function(x) do.call(rbind, x))
trunc <- do.call(rbind, b)
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
