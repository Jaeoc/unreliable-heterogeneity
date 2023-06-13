# Project: Unexplained heterogeneity
# Code purpose: Creating plots from simulated data
# Author: Anton Olsson-Collentine (anton@olssoncollentine.com)


#****************************************
# functions and packages
#****************************************

source("./code/02_functions.r") # for plot prep functions

library(ggplot2) #for plotting

#****************************************
# General plot prep
#****************************************


dat <- read.csv("data/means_combined.csv")


# Add real truncated values
es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.1, 0.15, 0.2)

tau_trunc <- lapply(tau, compute_var_truncated, mu = es)
tau_trunc <- unlist(tau_trunc)
trunc <- data.frame(mu = es,
                    nominal_tau = rep(tau, each = length(es)),
                    true_tau = tau_trunc)


#Prepping plot
dat$tau_hat <- sqrt(dat$tau2_hat)
dat$nominal_tau <- sqrt(dat$true_tau2)

dat$k <- factor(dat$k, levels = c("5", "20", "40", "200"))
dat$N <- factor(dat$sample_size, levels = c("50", "100", "150", "200"))

#****************************************
# Main manuscript plots
#****************************************

# Figure 1 (N = 150, k = 20, r, r_z, and HS)
#****************************************
dat_20_150 <- dat[dat$k == 20 & dat$N == 150 &
                 dat$reliability_mean < 1 &
                 dat$reliability_sd == 0.15,]

dat_20_150$reliability_mean <- factor(dat_20_150$reliability_mean)


# For fisher's z there is no truncation,
# create dataframe with all true heterogeneity levels
trunc$effect_type  <- "Pearson's r"
hs_true <- trunc
hs_true$effect_type <- "Hunter & Schmidt" #effect type is r, but named this way for plotting

z_true <- trunc
z_true$effect_type <- "Fisher's z"
dat_z <- dat_20_150[dat_20_150$effect_type == "r_z",]
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
# This helps with the facet_grid to make the plot look nice and helps us subset the tau-values of interest

dat_20_150 <- dat_20_150[dat_20_150$nominal_tau %in% c(0, 0.1, 0.15, 0.2), ]

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
geom_line(aes(linetype = reliability_mean), show.legend = TRUE) +
scale_linetype_manual(values = 2:5) +
guides(linetype = guide_legend(reverse = TRUE, title = "Mean Reliability")) +
geom_line(aes(x = mu, y = true_tau), data = true_tau_line, linetype = 1) +
expand_limits(y = 0) +
ylab(expression("Between-studies standard deviation "~tau)) +
xlab(expression("Average effect size "~mu)) +
facet_grid(effect_type~nominal_tau, scales = "free", switch = "y") +
theme_bw()

#ggsave("figures/manuscript/z-r-hs-plot.png", width = 8.62, height = 9.93)


# Figure 2 (r, N = 150, k variable, tau low)
#****************************************

dat_low <- dat[dat$N == 150 &
               dat$method == "HV" &
               dat$effect_type == "r" &
               dat$nominal_tau %in% c(0.02, 0.04, 0.06, 0.08) &
               dat$reliability_mean < 1 &
               dat$reliability_sd == 0.15,]

dat_low$reliability_mean <- factor(dat_low$reliability_mean)

# Truncated values (different from Figure 1 as tau is different)
es <- seq(0, 0.6, 0.1)
tau <- c(0.02, 0.04, 0.06, 0.08)

tau_trunc <- lapply(tau, compute_var_truncated, mu = es)
tau_trunc <- unlist(tau_trunc)
trunc_low <- data.frame(mu = es,
                    nominal_tau = rep(tau, each = length(es)),
                    true_tau = tau_trunc)

trunc_low$true_tau <- trunc_low$nominal_tau
trunc_low$tau_hat <- trunc_low$tau



#Plot with fixed N and varying k, and small tau-values
ggplot(dat_low, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = reliability_mean), show.legend = TRUE) +
scale_linetype_manual(values = 2:5) +
geom_line(aes(x = mu, y = true_tau), data = trunc_low, linetype = 1) +
guides(linetype = guide_legend(reverse = TRUE, title = "Mean Reliability")) +
expand_limits(y = 0) +
ylab(expression("Between-studies standard deviation "~tau)) +
xlab(expression("Average effect size "~mu)) +
facet_grid(k~nominal_tau, scales = "free") +
theme_bw()


#ggsave("figures/manuscript/r_tau_0.02-0.08.png", width = 8.62, height = 9.93)

#****************************************
# Supplement A plots
#****************************************


# Supplement A Figure 1 ((variance in reliabilities))
#****************************************

dat_a <- dat[dat$k == 20 &
             dat$sample_size == 150 &
             dat$reliability_mean < 1 &
             dat$method == "HV" &
             dat$effect_type == "r" &
             dat$nominal_tau %in% c(0, 0.1, 0.15, 0.2),]

dat_a$reliability_mean <- factor(dat_a$reliability_mean)
dat_a$reliability_sd <- factor(dat_a$reliability_sd,
                               levels = c(0.15, 0.1, 0.05, 0))

# Truncated values (same as in Figure 1, row Pearson's r)
es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.1, 0.15, 0.2)

tau_trunc <- lapply(tau, compute_var_truncated, mu = es)
tau_trunc <- unlist(tau_trunc)
trunc_a <- data.frame(mu = es,
                    nominal_tau = rep(tau, each = length(es)),
                    true_tau = tau_trunc)

trunc_a$true_tau <- trunc_a$nominal_tau
trunc_a$tau_hat <- trunc_a$tau


# NB! the object true_tau_line below comes from Figure 1 main manuscript section

ggplot(dat_a, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = reliability_mean), show.legend = TRUE) +
scale_linetype_manual(values = 2:5) +
geom_line(aes(x = mu, y = true_tau), data = trunc_a, linetype = 1) +
guides(linetype = guide_legend(reverse = TRUE, title = "Mean Reliability")) +
expand_limits(y = 0) +
ylab(expression("Between-studies standard deviation "~tau)) +
xlab(expression("Average effect size "~mu)) +
facet_grid(reliability_sd~nominal_tau, scales = "free") +
theme_bw()

#ggsave("figures/supplement/supplement_A_fig1.png", width = 8.62, height = 9.93)


# Supplement A figure 2 (variable k)
#****************************************
dat_b <- dat[dat$sample_size == 150 &
             dat$reliability_sd == 0.15 &
             dat$nominal_tau %in% c(0, 0.1, 0.15, 0.2) &
             dat$method == "HV" &
             dat$effect_type == "r",]

dat_b$reliability_mean <- factor(dat_b$reliability_mean)
dat_b$k <- factor(dat_b$k)

# Truncated values (same as in supplement A)
es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.1, 0.15, 0.2)

tau_trunc <- lapply(tau, compute_var_truncated, mu = es)
tau_trunc <- unlist(tau_trunc)
trunc_b <- data.frame(mu = es,
                    nominal_tau = rep(tau, each = length(es)),
                    true_tau = tau_trunc)

trunc_b$true_tau <- trunc_b$nominal_tau
trunc_b$tau_hat <- trunc_b$tau


# NB! the object true_tau_line below comes from Figure 1 main manuscript section

ggplot(dat_b, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = reliability_mean), show.legend = TRUE) +
scale_linetype_manual(values = 2:5) +
geom_line(aes(x = mu, y = true_tau), data = trunc_b, linetype = 1) +
guides(linetype = guide_legend(reverse = TRUE, title = "Mean Reliability")) +
expand_limits(y = 0) +
ylab(expression("Between-studies standard deviation "~tau)) +
xlab(expression("Average effect size "~mu)) +
facet_grid(k~nominal_tau, scales = "free") +
theme_bw()

#ggsave("figures/supplement/supplement_A_fig2.png", width = 8.62, height = 9.93)


# Supplement A Figure 3 (variable sample size)
#****************************************
dat_c <- dat[dat$k == 20 &
             dat$reliability_sd == 0.15 &
             dat$nominal_tau %in% c(0, 0.1, 0.15, 0.2) &
             dat$method == "HV" &
             dat$effect_type == "r",]


dat_c$reliability_mean <- factor(dat_c$reliability_mean)
dat_c$N <- factor(dat_c$N)

# Truncated values (same as in supplement A)
es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.1, 0.15, 0.2)

tau_trunc <- lapply(tau, compute_var_truncated, mu = es)
tau_trunc <- unlist(tau_trunc)
trunc_c <- data.frame(mu = es,
                    nominal_tau = rep(tau, each = length(es)),
                    true_tau = tau_trunc)


# NB! the object true_tau_line below comes from Figure 1 main manuscript section

ggplot(dat_c, aes(x = mu, y = tau_hat)) +
geom_line(aes(linetype = reliability_mean), show.legend = TRUE) +
scale_linetype_manual(values = 2:5) +
geom_line(aes(x = mu, y = true_tau), data = trunc_c, linetype = 1) +
guides(linetype = guide_legend(reverse = TRUE, title = "Mean Reliability")) +
expand_limits(y = 0) +
ylab(expression("Between-studies standard deviation "~tau)) +
xlab(expression("Average effect size "~mu)) +
facet_grid(N~nominal_tau, scales = "free") +
theme_bw()

#ggsave("figures/supplement/supplement_A_fig3.png", width = 8.62, height = 9.93)

#****************************************
# Supplement B plots
#****************************************

# Figure B1. (z-r-hs, perfect reliability)
#****************************************


dat_perfect <- dat[dat$k == 20 &
                   dat$N == 150 &
                   dat$reliability_mean == 1,]

dat_perfect$reliability_mean <- factor(dat_perfect$reliability_mean)

# truncate the fisher z nominal tau-values to match the pearson r ones
#this is purely for plotting purposes. The true value can be seen from the black line in plot.
dat_perfect$nominal_tau <- ifelse(dat_perfect$effect_type == "r_z",
                              trunc(dat_perfect$nominal_tau *100) / 100,
                              dat_perfect$nominal_tau)


# What I've done in the above to lines of code is to
# 1) put the 'truncated' true_tau for r_z to fisher z nominal value
# 2) put the 'true tau' in the main data object to the nominal r value
# This helps with the facet_grid to make the plot look nice and helps us subset the tau-values of interest

dat_perfect <- dat_perfect[dat_perfect$nominal_tau %in% c(0, 0.1, 0.15, 0.2), ]

# More plot polish
dat_perfect$effect_type <- ifelse(dat_perfect$effect_type == "r_z",
                                "Fisher's z",
                                "Pearson's r")
dat_perfect$effect_type <- ifelse(dat_perfect$effect_type == "Pearson's r" &
                                 dat_perfect$method == "HS",
                                "Hunter & Schmidt",
                                dat_perfect$effect_type)

dat_perfect$effect_type <- factor(dat_perfect$effect_type,
                            levels = c("Fisher's z",
                                        "Pearson's r",
                                        "Hunter & Schmidt"))

dat_perfect$mu <- ifelse(dat_perfect$effect_type == "Fisher's z",
                        tanh(dat_perfect$mu),
                        dat_perfect$mu)

# NB! the object true_tau_line below comes from Figure 1 main manuscript section

ggplot(dat_perfect, aes(x = mu, y = tau_hat)) +
geom_line(linetype = 2) +
geom_line(aes(x = mu, y = true_tau), data = true_tau_line, linetype = 1) +
expand_limits(y = 0) +
ylab(expression("Between-studies standard deviation "~tau)) +
xlab(expression("Average effect size "~mu)) +
facet_grid(effect_type~nominal_tau, scales = "free", switch = "y") +
theme_bw()

#ggsave("figures/supplement/supplement_B_fig1.png", width = 8.62, height = 9.93)


# Figure B2
#****************************************
library(data.table)
library(parallel)

simulate_perfect_r <- function(

    k, #number of studies
    sample_size, #sample size, fixed across studies
    true_tau2, #variance of superpopulation
    mu){ #these are the default values in the function: https://www.metafor-project.org/doku.php/tips:convergence_problems_rma

    # Output is a dataframe with selected results from a metafor::rma object


    # Draw rho from mean rho and compute sampling variances

    rho <- rnorm_truncated(n = k, mean = mu, sd = sqrt(true_tau2),
                            lower_bound = -1, upper_bound = 1)
    sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)

    # given study rho, draw sample rho
    r_se <- rnorm_truncated(k, mean = rho, sd = sqrt(sampling_var_rho),
                                lower_bound = -1, upper_bound = 1)

    #compute observed sampling variance estimate given measurement error
    r_var_estimate <- (1-r_se^2)^2 / (sample_size -1)


    # fit meta-analysis on observed values
    fit <- rma(yi = r_se, vi = sample(r_var_estimate))



        data.frame(intercept_b = fit$b,
                   intercept_p = fit$pval,
                   tau2_hat = fit$tau2,
                   tau2_p = fit$QEp)

}

# Prep simulation
sample_size <- 150
k <- 20
mu <- seq(from = 0, to = 0.6, by = 0.1)
true_tau <- c(0, 0.1, 0.15, 0.2)
true_tau2 <- true_tau^2

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu)


reps <- 1e3
out_list <- vector("list", length = nrow(cond))

ncores <-parallel::detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers.

clusterEvalQ(cl, library(metafor))
clusterExport(cl, c("simulate_perfect_r", "rnorm_truncated", "try_run"))


# Run simulation

for(r in 1:nrow(cond)){

mes <- paste0("\n now on condition ", r)
    cat(mes)
    task <- function(iteration, cond_r){ #anonymous function needed when using for replications
                simulate_perfect_r(k = cond_r$k,
                            sample_size = cond_r$sample_size,
                            true_tau2 = cond_r$true_tau2,
                            mu = cond_r$mu)
        }

        out_list[[r]] <- parLapply( #function from parabar
                    cl = cl, # cluster
                    X = 1:reps, #looping over
                    fun = task,
                    cond_r = cond[r, ]
    )

    #Compute the mean across replications for the condition and add condition identifiers
    out_list[[r]] <- rbindlist(out_list[[r]]) #function from data.table
    out_list[[r]] <- out_list[[r]][, lapply(.SD, mean)] #data.table colMeans but returns a dataframe (well, data.table)
    out_list[[r]] <- cbind(out_list[[r]], cond[r,])

}
stopCluster(cl)


# prep plot


e <- rbindlist(out_list)
e$tau_hat <- sqrt(e$tau2_hat)
e$nominal_tau <- sqrt(e$true_tau2)

# Add real truncated values
es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.1, 0.15, 0.2)

tau_trunc <- lapply(tau, compute_var_truncated, mu = es)
tau_trunc <- unlist(tau_trunc)
trunc <- data.frame(mu = es,
                    nominal_tau = rep(tau, each = length(es)),
                    true_tau = tau_trunc)


# Plot D2.

ggplot(e) +
geom_line(aes(x = mu, y = tau_hat), linetype = 2) +
geom_line(aes(x = mu, y = true_tau), data = trunc, linetype = 1) +
facet_wrap(~nominal_tau) +
ylab(expression("Between-studies standard deviation "~tau)) +
xlab(expression("Average effect size "~mu))

#ggsave("figures/supplement/supplement_B_fig2.png", width = 8.62, height = 9.93)
