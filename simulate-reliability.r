# Project: Unexplained heterogeneity
# Code purpose: simulate the effect of reliability on heterogeneity
# Author: Anton Olsson-Collentine (anton@olssoncollentine.com)

#********************************************************************
# fixed effect
#********************************************************************

# We use the classical test theory framework for this

# Start with a single study
# Step 1: simulate a true score

# Latent score, translated into scale score?
# Or start with a binary score?
# Can also simulate multiple binary scores that sum to a scale score
# but then why not simulate the scale score directly?

#for reliability we must have more than one item..


#variance is prob * (1 - prob) = 0.25 (maximal) for each item

# Reliability is the proportion true score variance of the observed variance. Measurement error by definition has a mean of zero.

# Let's first imagine we are measuring some unrestrained variable
true_var <- 1
true_mean <- 0
n  <- 10
error_var <- 0.5


true <- rnorm(n = n, mean = true_mean, sd = sqrt(true_var))
e <- rnorm(n = n, mean = 0, sd = sqrt(error_var))
x <- true + e

var_x <- true_var + error_var
reliability <- true_var / var_x

# We would probably prefer to go from reliability to measurement error
#Assuming true var = 1 and
reliability <- 0.8

var_x <- true_var / reliability
#and
error_var <- var_x - true_var

calc_error_var <- function(true_var, reliability){
    var_x <- true_var / reliability
    error_var <- var_x - true_var
    data.frame(true_var, reliability, var_x, error_var)
}

reliability <- seq(from = 0, to = 1, by = 0.1)

reference_values <- calc_error_var(true_var, reliability)

#So for 0.8 reliability we have 0.25 measurement var and our observed data:
true <- rnorm(n = n, mean = true_mean, sd = sqrt(true_var))
e <- rnorm(n = n, mean = 0, sd = sqrt(error_var))
x <- true + e
dat <- data.frame(true, e, x)

#This gives us the scores of one group on a variable
#To have an effect size we must compare with another variable
#easiest perhaps is correlation


#assuming var = 1 for both variables and cor/covar = 0
n  <- 10
sigma <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
true <- MASS::mvrnorm(n = n, mu = rep(0, 2), Sigma = sigma)
#Given reliability = 0.8
error_var <- 0.25
e <- matrix(rnorm(2*n, mean = 0, sd = sqrt(error_var)), ncol = 2, byrow = TRUE)
x <- matrix(c(true[,1] + e[,1], true[,2] + e[,2]), ncol = 2)

dat <- cbind(true, e, x)
colnames(dat) <- paste0(rep(c("true", "e", "x"), each =  2), "_", rep(1:2, 3))
#we have now the data for one study, correlational
#true score, error and observed score for reliability = 0.8 and true score variance = 1

#Next, we would like to compute this data for multiple studies with the same
#reliability and check so that we recuperate the true effect

# After that we use a fixed effect, but vary the reliabilities to check if we can recuperate the true effect if we add the true reliabilities in a meta-regression

# Note that this does not so far take into account that in real life the estimate of reliability is affected by sampling error. This would require simulating item-level data. So far I have only simulated some kind of continuous score, to avoid problems of range restrictions etc.

simulate_study <- function(true_var = 1, covar, reliability, n){

#draw true scores
sigma <- matrix(c(true_var, covar, covar, true_var), nrow = 2, ncol = 2)
true <- MASS::mvrnorm(n = n, mu = rep(0, 2), Sigma = sigma)

#Add error
e_dat <- calc_error_var(true_var, reliability)
error_var <- e_dat$error_var

e <- matrix(rnorm(2*n, mean = 0, sd = sqrt(error_var)), ncol = 2, byrow = TRUE)
x <- matrix(c(true[,1] + e[,1], true[,2] + e[,2]), ncol = 2)

#combine and output
dat <- cbind(true, e, x)
colnames(dat) <- paste0(rep(c("true", "e", "x"), each =  2), "_", rep(1:2, 3))
dat <- as.data.frame(dat)

dat # out
}

sample_size <- 10
k  <- 5
dat <- replicate(n = k,
    simulate_study(covar = 0, reliability = 0.8, n = sample_size),
    simplify = FALSE)

# To meta-analyze the studies I need a) a correlation and b) its sampling variance
#the true sampling variance of a single mean is simply
mean_var <- true_var / sample_size
#Variance of a correlation is
var_r <- (1-r^2)^2 / (n1 + n2 -1) #Cooper and Hedges (2009), equation 12.27
#for the true sampling variance we would put r to the known value (zero)

dat_cors <- lapply(dat, function(x) cor(x[,c(1:2)]))
dat_cors <- lapply(dat_cors, function(x) x[1, 2])
dat_cors <- do.call(rbind, dat_cors)

var_cors <- (1-dat_cors^2)^2 / (sample_size*2 -1)
res <- data.frame(r = dat_cors, v_i = var_cors)

var_cors_true <-rep((1-0^2)^2 / (sample_size*2 -1), k)
res_true <- data.frame(r = dat_cors, v_i = var_cors_true)
#Meta-analyzing the true scores. That is,
# a meta-analysis with no measurement error, but with sampling error
library(metafor)

fit_true1 <- rma(yi = r, vi = v_i,  data = res)
fit_true2 <- rma(yi = r, vi = v_i,  data = res_true)

#If I use ~1e4 or so participants for each study I should have basically no
#sampling error either and should recuperate the ES perfectly.
# EDIT: mostly true, but it is better to increase the k? Unexpected.
sample_size <- 150
k  <- 100
dat <- replicate(n = k,
    simulate_study(covar = 0, reliability = 0.8, n = sample_size),
    simplify = FALSE)

dat_cors <- lapply(dat, function(x) cor(x[,c(1:2)]))
dat_cors <- lapply(dat_cors, function(x) x[1, 2])
dat_cors <- do.call(rbind, dat_cors)
var_cors <- (1-dat_cors^2)^2 / (sample_size*2 -1)

res <- data.frame(r = dat_cors, v_i = var_cors)

fit <- rma(yi = r, vi = v_i,  data = res)
# Heterogeneity is here indeed essentially zero, as well as the mean estimate.

#If I use ~1e4 participants but the correlations and estimated sampling variance based on the observed scores I should have only measurement error. Given this I should still recuperate the true effect well since the measurement error is random.
compute_cor <- function(data, data_type = c("true", "x")){

    cols <-grepl(pattern = data_type, x = names(data))
    dat_cor <- cor(data[, cols])

    dat_cor[1,2] #out a single cor instead of matrix
}

summarize_cors <- function(data_list, data_type = c("true", "x")){

dat_cors <- lapply(data_list, compute_cor, data_type = data_type)
dat_cors <- do.call(rbind, dat_cors)

var_cors <- (1-dat_cors^2)^2 / (sample_size -1)

data.frame(r = dat_cors, v_i = var_cors) #out

}

sample_size <- 50
k  <- 5
dat <- replicate(n = k,
    simulate_study(covar = 0, reliability = 0.8, n = sample_size),
    simplify = FALSE)

res <- summarize_cors(dat, data_type = "x")

fit <- rma(yi = r, vi = v_i,  data = res)

#Even with n = 50 and k = 5 we get pretty close to zero.
#Basically, with a fixed effect effect size and reliability we have very little problems in general. What happens if we decrease reliability?

sample_size <- 150
k  <- 20
dat <- replicate(n = k,
    simulate_study(covar = 0.6, reliability = 0.8, n = sample_size),
    simplify = FALSE)

res <- summarize_cors(dat, data_type = "x")

fit <- metafor::rma(yi = r, vi = v_i,  data = res)

# Even with low reliability, as long as the ES is zero, we find a more or less zero effect estimate. Hoever, if the correlation is > 0, the estimate is much further off, even if the reliability is 0.8.

#Now, let's focus on the heterogeneity, by adding variable reliability to the fixed effect.

sample_size <- 150
k <- 20
reliabilities <- runif(n = k, min = 0.3, max = 0.9)
reliabilities <- c(rep(0.3, 5), rep(0.5, 5), rep(0.7, 5), rep(0.8, 5))
dat <- lapply(reliabilities, function(x){
    simulate_study(covar = 0.8, n = sample_size, reliability = x)
    }
)

res <- summarize_cors(dat, data_type = "x")

fit <- metafor::rma(yi = r, vi = v_i,  data = res)
fit
#Given a set sample_size, k, and reliabilities, we can do a plot showing the tau for different effect sizes
#Although we don't actually need to simulate this, we can just compute it based on heterogeneity formula

sample_size <- 150
k <- 20
ES <- seq(from = 0, to = 1, by = 0.1)
reliabilities <- c(rep(0.3, 5), rep(0.5, 5), rep(0.7, 5), rep(0.9, 5))

reps <- 1e3
tau <- rep(NA, reps)

out <- vector("list", length = length(ES))
names(out) <- paste0("rho_", ES)

start <- Sys.time()
for(e in 1:length(ES)){

    for(r in seq_len(reps)){
    dat <- lapply(reliabilities, function(x){
     simulate_study(covar = ES, n = sample_size, reliability = x)
     }
    )
    res <- summarize_cors(dat, data_type = "x")
    fit <- metafor::rma(yi = r, vi = v_i,  data = res)
    tau[r] <- sqrt(fit$tau2)
    }
    out[[e]] <- mean(tau)
    cat(ES[e], "\n")
}
 end <- Sys.time()
 end - start

rel_0.3_0.9 <- out
#rel_0.6_0.8 <- out


rel_0.3_0.9_df <- do.call(rbind, rel_0.3_0.9)
rel_0.3_0.9_df <- cbind(rel_0.3_0.9_df, ES)
colnames(rel_0.3_0.9_df)[1] <- "tau_hat"
rel_0.3_0.9_df <- as.data.frame(rel_0.3_0.9_df)
rel_0.3_0.9_df$id <- "R between 0.3 - 0.9"

rel_0.6_0.8_df <- do.call(rbind, rel_0.6_0.8)
rel_0.6_0.8_df <- cbind(rel_0.6_0.8_df, ES)
colnames(rel_0.6_0.8_df)[1] <- "tau_hat"
rel_0.6_0.8_df <- as.data.frame(rel_0.6_0.8_df)
rel_0.6_0.8_df$id <- "R between 0.6 - 0.8"

plot_dat <- rbind(rel_0.6_0.8_df, rel_0.3_0.9_df)

library(ggplot2)

ggplot(data = plot_dat) +
geom_line(aes(x = ES, y = tau_hat, linetype = id))

#Strange results,not monotone..
# Hmm, I think it's because I only do one random draw of reliabilities, so there is too much variability in it. EDIT: no,I suspect something wrong with my calculations.

#Alternative version: 1) decide a rho, 2) draw an r for each study with measurement error
rho <- 0.5
r_me <- rho*sqrt(reliabilities)^2
data.frame(rho, r_me) #r_x is for infinite sample size atm
sample_size <- 150
sampling_var_rho <- (1-rho^2)^2 / (sample_size -1) #rho or r, or r_me?
r_me_se <- rnorm(k, mean = r_me, sd = sqrt(sampling_var_rho))
data.frame(rho, r_me, r_me_se)
r_var_estimate <- (1-r_me_se^2)^2 / (sample_size -1)
res <- data.frame(r = r_me_se, v_i = r_var_estimate)

fit <- metafor::rma(yi = r, vi = v_i,  data = res)

#Joran also suggested to just draw a rho from a normal distribution, then the
#variance is measurement error variance
#then can just use that as the mean to draw a new effect size with some sampling
#error? But how to decide the measurement error variance in the first step?

#****************************
# Effect size and tau hat
#****************************
#As a loop
sample_size <- 150
k <- 20
rho <- seq(from = 0, to = 1, by = 0.1)
#reliabilities <- c(rep(0.3, 5), rep(0.5, 5), rep(0.7, 5), rep(0.9, 5))
#reliabilities <- c(rep(0.6, 5), rep(0.65, 5), rep(0.75, 5), rep(0.8, 5))
#These weren't necessary after all
reliability_min <- 0.4
reliability_max <- 0.6
reps <- 1e3
tau <- rep(NA, reps)

out <- vector("list", length = length(rho))
names(out) <-  rho
library(metafor)

start <- Sys.time()
for(e in seq_along(rho)){

    sampling_var_rho <- (1-rho[e]^2)^2 / (sample_size -1) #rho or r, or r_me?

    for(r in seq_len(reps)){
    # given population rho, draw sample rho
    r_se <- rnorm(k, mean = rho[e], sd = sqrt(sampling_var_rho))

    #Compute observed r given measurement error
    reliabilities <- runif(k, min = reliability_min, max = reliability_max) #draw reliabilities
    r_se_me <- r_se*sqrt(reliabilities)^2 #number of reliabilities values gives k
    r_var_estimate <- (1-r_se_me^2)^2 / (sample_size -1)

    fit <- rma(yi = r_se_me, vi = r_var_estimate)
    tau[r] <- sqrt(fit$tau2)
}
    out[[e]] <- mean(tau)

    cat(rho[e], "\n")
}
end <- Sys.time()
end - start


#rel_0_1 <- out
#rel_0.1_0.9 <- out
#rel_0.2_0.8 <- out
#rel_0.3_0.7 <- out
rel_0.4_0.6 <- out
results <- list(rel_0_1, rel_0.1_0.9, rel_0.2_0.8, rel_0.3_0.7, rel_0.4_0.6)

#As far as I can tell, I have done something wrong when simulating individual
#level data, everything seems fine at the effect size level..

clean_res <- function(res, rho){
    a <- do.call(rbind, res)
    a <- cbind(a, rho)
    colnames(a)[1] <- "tau_hat"
    a <- as.data.frame(a)
    a #out
}
id_values <- paste0("reliability between ", 0:4 / 10, " - ", seq(from = 1, to = 0.6, by = -0.1))
id_values <- rep(id_values, each = 11)

results <- lapply(results, clean_res, rho = rho)

plot_dat <- do.call(rbind, results)
plot_dat$id <- id_values

saveRDS(plot_dat, "rho_tau.RDS")
library(ggplot2)

ggplot(data = plot_dat) +
geom_line(aes(x = rho, y = tau_hat, linetype = id))
ggsave("rho_tau_plot.png")

# SAME but with smaller k (EDIT: doesn't really change, deleted section)

#******************************
#fixed effect, no meta-regression----
#******************************
sample_size <- 150
k <- 20
rho <- 0.5
sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)
reliability_min <- 0.3
reliability_max <- 0.9

reps <- 1e2
out <- vector("list", length = reps)

start <- Sys.time()
for(r in seq_len(reps)){
    # given population rho, draw sample rho
     r_se <- rnorm(k, mean = rho, sd = sqrt(sampling_var_rho))

    #Compute observed r given measurement error
    reliabilities <- runif(k, min = reliability_min, max = reliability_max) #draw reliabilities to avoid some strange draw
    r_se_me <- r_se*sqrt(reliabilities)^2 #number of reliabilities values gives k
    r_var_estimate <- (1-r_se_me^2)^2 / (sample_size -1)

    fit <- rma(yi = r_se_me, vi = r_var_estimate)
    #Shows that if measurement reliability is zero then the estimate is zero, but if perfect, then estimate would be 0.5
    meta_b <- fit$b
    meta_p <- fit$pval
    meta_tau2 <- fit$tau2

    out[[r]] <- data.frame(intercept_b = meta_b[1],
                           intercept_p = meta_p[1],
                           meta_tau2)
}

end <- Sys.time()
end - start

a <- do.call(rbind, out)
(a_mean <- colMeans(a))
round(a_mean, 3)

#Fixed effect, estimate very little heterogeneity, tau = 0.058 and tau2 = 0.003
# This despite variance in reliabilities
#******************************
#Meta-regression, fixed effect----
#******************************
sample_size <- 150
k <- 20
rho <- 0.5
sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)
reliability_min <- 0.3
reliability_max <- 0.9

reps <- 1e2
out <- vector("list", length = reps)

start <- Sys.time()
for(r in seq_len(reps)){
    # given population rho, draw sample rho
     r_se <- rnorm(k, mean = rho, sd = sqrt(sampling_var_rho))

    #Compute observed r given measurement error
    reliabilities <- runif(k, min = reliability_min, max = reliability_max) #draw reliabilities to avoid some strange draw
    r_se_me <- r_se*sqrt(reliabilities)^2 #number of reliabilities values gives k
    r_var_estimate <- (1-r_se_me^2)^2 / (sample_size -1)

    fit <- rma(yi = r_se_me, vi = r_var_estimate, mods = ~ reliabilities)
    #Shows that if measurement reliability is zero then the estimate is zero, but if perfect, then estimate would be 0.5
    meta_b <- fit$b
    meta_p <- fit$pval
    meta_tau2 <- fit$tau2

    out[[r]] <- data.frame(intercept_b = meta_b[1],
                           intercept_p = meta_p[1],
                           reliability_b = meta_b[2],
                           reliability_p = meta_p[2],
                           meta_tau2)
}

end <- Sys.time()
end - start

a <- do.call(rbind, out)
colMeans(a)
mean(a$intercept_p < 0.05)
mean(a$reliability_p < 0.05)
#Evidence that works for fixed effect. Now we try for a random effect
#*****************************************
#RANDOM EFFECT, no meta-regression----
#*****************************************
sample_size <- 150
k <- 20
true_tau2 <- 0.1
mu <- 0.1
reliability_min <- 0.6
reliability_max <- 1

reps <- 1e2
out <- vector("list", length = reps)

start <- Sys.time()

for(r in seq_len(reps)){

    # Draw rho from mean rho
    rho <- rnorm(n = k, mean = mu, sd = sqrt(true_tau2))
    sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)

    # given study rho, draw sample rho
    r_se <- rnorm(k, mean = rho, sd = sqrt(sampling_var_rho))

    #Compute observed r given measurement error
    reliabilities <- runif(k, min = reliability_min, max = reliability_max) #draw reliabilities to avoid some strange draw
    r_se_me <- r_se*sqrt(reliabilities)^2

    r_var_estimate <- (1-r_se_me^2)^2 / (sample_size -1) #observed estimate

    fit <- rma(yi = r_se_me, vi = r_var_estimate)

    meta_b <- fit$b
    meta_p <- fit$pval
    meta_tau2 <- fit$tau2

    out[[r]] <- data.frame(intercept_b = meta_b[1],
                           intercept_p = meta_p[1],
                           meta_tau2)
}
end <- Sys.time()
end - start

b <- do.call(rbind, out)
(b_means <- colMeans(b))
round(b_means, 3)
#tau2 = 0.1
#but our estimate of tau2 =~ 0.054
#or in tau: tau =~0.31, but our estimate is 0.22

#That is, we underestimate the true heterogeneity

# the problem with heteroscedasticity is in the standard errors
mean(b$intercept_p < 0.05)
mean(b$reliability_p < 0.05)

# As we can see, the reliability is consistently identified as non-significant
# despite being a true explanation for heterogeneity
#*****************************************
# random effects, correct for measurement error
#****************************
#Here I will just pretend that we know the measurement error perfectly (i.e., meta-analyze the values without any measurement error).

#EDIT: haven't worked on this section, because I seem to be getting the same amount of heterogeneity regardless whether I have measurement error, which completely indicates something is wrong with my code..

# EDIT2: Now ran the code in different terminals, I think there was some problem with the same variables being saved... Now I am getting the expected results again: without measurement error we recuperate the average effect and the tau2 perfectly, with (uncorrected for) measurement error we underestimate both the effect size (expected from Schmidt & Hunter) and the heterogeneity (unexpected).

#That is, I have confirmed that heterogeneity seems to be underestimated under a (true) random effects model with varying reliabilities. Need to check if someone has made this point yet.

sample_size <- 150
k <- 20
true_tau2 <- 0.1
mu <- 0.5
reliability_min <- 0.7
reliability_max <- 0.9

reps <- 1e2
out <- vector("list", length = reps)

start <- Sys.time()

for(r in seq_len(reps)){

    # Draw rho from mean rho
    rho <- rnorm(n = k, mean = mu, sd = sqrt(true_tau2))
    sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)

    # given study rho, draw sample rho
    r_se <- rnorm(k, mean = rho, sd = sqrt(sampling_var_rho))

    r_var_estimate <- (1-r_se^2)^2 / (sample_size -1) #observed estimate

    fit <- rma(yi = r_se, vi = r_var_estimate)

    meta_b <- fit$b
    meta_p <- fit$pval
    meta_tau2 <- fit$tau2

    out[[r]] <- data.frame(intercept_b = meta_b[1],
                           intercept_p = meta_p[1],
                           meta_tau2)
}
end <- Sys.time()
end - start

e <- do.call(rbind, out)
(e_means <- colMeans(b))
round(e_means, 3)

#We get a very similar level of heterogeneity estimated (0.047)
#EDIT: in fact exactly identical. That cannot be right.. If so, measurement error does nothing..

#****************************************
# Turn into function
#****************************************
sample_size <- 150
k <- 20
true_tau2 <- 0.1
mu <- 0.1
reliability_min <- 0.6
reliability_max <- 1


simulate_rma <- function(k, sample_size, true_tau2, mu, reliability_min, reliability_max){
 # Draw rho from mean rho
    rho <- rnorm(n = k, mean = mu, sd = sqrt(true_tau2))
    sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)

    # given study rho, draw sample rho
    r_se <- rnorm(k, mean = rho, sd = sqrt(sampling_var_rho))

    #Compute observed r given measurement error
    reliabilities <- runif(k, min = reliability_min, max = reliability_max) #draw reliabilities to avoid some strange draw
    r_se_me <- r_se*sqrt(reliabilities)^2

    r_var_estimate <- (1-r_se_me^2)^2 / (sample_size -1) #observed estimate

    fit <- rma(yi = r_se_me, vi = r_var_estimate)

    data.frame(intercept_b = fit$b,
               intercept_p = fit$pval,
               tau2_hat = fit$tau2,
               tau2_p = fit$QEp)
}

#****************************************
# Simulations of underestimated heterogeneity
#****************************************
library(metafor)
# Suspicion 1, it is the average reliability that matters, not the range
sample_size <- 150
k <- 20
true_tau2 <- 0.1
mu <- 0.1
reliability_min <- 0.6
reliability_max <- 1

reps <- 1e3

start <- Sys.time()

# Simple replicated simulation, one condition
out <- replicate(n = reps,
                simulate_rma(k = k,
                            sample_size = sample_size,
                            true_tau2 = true_tau2,
                            mu = mu,
                            reliability_min = reliability_min,
                            reliability_max = reliability_max),
                simplify = FALSE)

end <- Sys.time()
end - start

b <- do.call(rbind, out)
(b_means <- colMeans(b))
round(b_means, 3)

# check if varying reliability matters
set.seed(126)

#a)
reliability_min <- 0.1
reliability_max <- 0.7

out <- replicate(n = reps,
                simulate_rma(k = k,
                            sample_size = sample_size,
                            true_tau2 = true_tau2,
                            mu = mu,
                            reliability_min = reliability_min,
                            reliability_max = reliability_max),
                simplify = FALSE)

b <- do.call(rbind, out)
(b_means <- colMeans(b))
a  <- round(b_means, 3)

#b)
reliability_min <- 0.4
reliability_max <- 0.4

out <- replicate(n = reps,
                simulate_rma(k = k,
                            sample_size = sample_size,
                            true_tau2 = true_tau2,
                            mu = mu,
                            reliability_min = reliability_min,
                            reliability_max = reliability_max),
                simplify = FALSE)

b <- do.call(rbind, out)
(b_means <- colMeans(b))
round(b_means, 3)

# Conclusion: Varying reliability increases heterogeneity, but only minimally
# It is largely still the same underestimate
#Possibly I should make a plot of this to drive home the point better

# Conclusion 2: As reliability increases, the underestimate of heterogeneity decreases


# Across conditions
# We simply make the conditions a vector of input and loop over them
# Checking the effect of increasing mu
sample_size <- 150
k <- 20
true_tau2 <- 0.1
mu <- c(0, 0.1, 0.3, 0.5)
reliability_min <- 0.5
reliability_max <- 1
reps <- 1e3

out_list <- vector("list", length = length(mu))


start <- Sys.time()

for(ES in seq_along(mu)){ #gives us a list of lists

    out_list[[ES]] <- replicate(n = reps,
                simulate_rma(k = k,
                            sample_size = sample_size,
                            true_tau2 = true_tau2,
                            mu = mu[ES],
                            reliability_min = reliability_min,
                            reliability_max = reliability_max),
                simplify = FALSE)

}
end <- Sys.time()
end - start

e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- paste0("mu = ", mu)
lapply(e_means, round, 3)

#First round: difference seems random, due to simulation error
# Second round supports this idea, increasing number of reps to 1e3 for third round
# third round: confirmed


# Finally, check the effect of different tau. Larger tau -> larger underestimates?

sample_size <- 150
k <- 20
true_tau2 <- c(0.1, 0.3, 0.5, 0.7)
mu <- 0.3
reliability_min <- 0.6
reliability_max <- 1
reps <- 1e2

out_list <- vector("list", length = length(true_tau2))


start <- Sys.time()

for(ES in seq_along(true_tau2)){ #gives us a list of lists

    out_list[[ES]] <- replicate(n = reps,
                simulate_rma(k = k,
                            sample_size = sample_size,
                            true_tau2 = true_tau2[ES],
                            mu = mu,
                            reliability_min = reliability_min,
                            reliability_max = reliability_max),
                simplify = FALSE)

}
end <- Sys.time()
end - start

e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- paste0("tau2 = ", true_tau2)
lapply(e_means, round, 3)
e_df <- do.call(rbind, .Last.value)
e_df[,3] - true_tau2

#conclusion: very clear pattern of increase, as tau2 increases the underestimate
# of tau2 increases as well

# In other words, all my predictions bore out.

#****************************************
# Simulating data for plots
#****************************************
# see supercomputer_simulator.r


#****************************************
# plotting
#****************************************

library(ggplot2)

dat1 <- readRDS("./data/mu_varying_rel.RDS")
dat1 <- lapply(dat1, function(x) as.data.frame(as.list(x)))

dat_mu <- dat1[grepl("mu", names(dat1))]
dat_mu <- data.table::rbindlist(dat_mu, idcol = "mu")

dat_mu$reliability <- gsub(".*& ", "", dat_mu$mu)
dat_mu$mu <- gsub("&.*", "", dat_mu$mu)
dat_mu$mu <- as.numeric(gsub("mu = ", "", dat_mu$mu))

dat_rel <- dat1[grepl("rel", names(dat1))]
dat_rel <- data.table::rbindlist(dat_rel, idcol = "Reliability")
dat_rel$reliability_range <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

#Mu plot
ggplot(dat_mu, aes(x = mu, y = tau2_hat)) +
geom_point() +
geom_text(label = round(dat_mu$tau2_hat, 3), nudge_y = 0.005) +
geom_line(aes(linetype = reliability)) +
coord_cartesian(ylim = c(0, 0.1)) +
geom_hline(yintercept = 0.1, linetype = "dashed") +
theme_bw() +
annotate("text", label = "tau2 = 0.1",
y = 0.06, x = 0, hjust = 0)

ggsave("figures/mu_tau2_rel.png", width = 9.3, height = 8.8)

#reliability plot
ggplot(dat_rel, aes(x = reliability_range, y = tau2_hat)) +
geom_point() +
geom_text(label = round(dat_rel$tau2_hat, 3), nudge_y = 0.005) +
geom_line(aes(group = 1)) +
coord_cartesian(ylim = c(0, 0.1)) +
scale_x_continuous(breaks = dat_rel$reliability_range) +
geom_hline(yintercept = 0.1, linetype = "dashed") +
theme_bw() +
annotate("text", label = "mu = 0.3 \n tau2 = 0.1 \n Mean reliability = 0.5",
y = 0.075, x = 0, hjust = 0)

# ggsave("figures/reliability_range_tau2.png", width = 1.6*6, height = 6)



# plot varying true tau2 and increasing reliability mean

dat3 <- readRDS("data/mu_rel_variability.RDS")

dat3 <- lapply(dat3, function(x) as.data.frame(as.list(x)))

dat_mu_rel <- data.table::rbindlist(dat3, idcol = "condition")

dat_mu_rel$mu <- gsub(" &.*", "", dat_mu_rel$condition)
dat_mu_rel$mu <- as.numeric(gsub("mu = ", "", dat_mu_rel$mu))

dat_mu_rel$reliability <- gsub(".* & min reliability = ", "", dat_mu_rel$condition)
dat_mu_rel$reliability <- paste0(dat_mu_rel$reliability, "-1")

line_labels <- dat_mu_rel[dat_mu_rel$mu == 0.9,]

ggplot(dat_mu_rel, aes(x = mu, y = tau2_hat)) +
geom_line(aes(linetype = reliability), show.legend = FALSE) +
#scale_linetype_manual(values = 10:1) +
geom_text(data = line_labels, aes(x = mu, y = tau2_hat, group = reliability), nudge_x = 0.025,
label = line_labels$reliability) +
#annotate("text", x = 0.885, y = 0.086, label = "Reliability", hjust = 0) +
geom_hline(yintercept = 0.1, linetype = "dashed") +
scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
theme_bw()

ggsave("figures/mu_increasing_reliability.png", width = 8.62, height = 9.93)


# plot varying mu and varying reliability

dat2 <- readRDS("data/tau2_increasing_rel.RDS")

dat2 <- lapply(dat2, function(x) as.data.frame(as.list(x)))

dat_tau2 <- data.table::rbindlist(dat2, idcol = "condition")

dat_tau2$true_tau2 <- gsub("&.*", "", dat_tau2$condition)
dat_tau2$true_tau2 <- as.numeric(gsub("tau2 = ", "", dat_tau2$true_tau2))

dat_tau2$reliability <- gsub(".* & reliability = ", "", dat_tau2$condition)

tau2_labels <- dat_tau2[dat_tau2$true_tau2 == 0.9,]

ggplot(dat_tau2, aes(x = true_tau2, y = tau2_hat)) +
geom_line(aes(linetype = reliability), show.legend = FALSE) +
#scale_linetype_manual(values = 10:1) +
geom_label(data = tau2_labels, aes(x = true_tau2, y = tau2_hat, group = reliability), nudge_x = 0.025,
label = tau2_labels$reliability) +
annotate("text", x = 0.885, y = 0.86, label = "Reliability", hjust = 0) +
geom_abline() +
scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
theme_bw()

ggsave("figures/tau2_increasing_reliability.png", width = 8.62, height = 9.93)
