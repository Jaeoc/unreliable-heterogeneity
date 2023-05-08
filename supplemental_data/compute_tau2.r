#****************************************
# Project: unexplained heterogeneity
# script purpose: Compute I2 values for set tau values and sample sizes
# for Pearson's r. Then use these I2 values and sample sizes to compute
# corresponding tau values for Fisher's z.
# code: Anton Olsson-Collentine
#****************************************

#****************************************
# Background
#****************************************

# Given N and tau, compute I2
#I2 = tau2 / s2 + tau2
#s2 = sum(w_i)*(k - 1) / (sum(wi_i)^2) - sum(w_i^2)

# Because of using Pearson's r, the sampling variance depends on the
# effect size and N. Consequently, I2 depends on rho, N, tau
# But because we have fixed sample size across studies, not on K

# We only want one tau-value for each N and tau, not for each ES as well
# Option 1: Compute I2 for ES = 0 and use that
# Option 2: Compute I2 for all ES and take average/median.
# I prefer option 1 as the simpler option, and also because I2 is not linear

#****************************************
# Functions
#****************************************
compute_s2 <- function(sigma2, k){

    wi <- rep((1 / sigma2), k)

    s2_num <- sum(wi)*(k - 1)
    s2_denom <- (sum(wi)^2) - sum(wi^2)
    s2_num / s2_denom #s2 out
}

compute_I2 <- function(vec = c(rho, tau, sample_size)){
    # #rho = average true effect size
    # #tau = SD in true effect sizes
    # #sample_size = fixed within-study sample size

    rho <- vec[["rho"]]
    tau <- vec[["tau"]]
    sample_size <- vec[["sample_size"]]

    tau2 <- tau^2
    k <- 20 #constant, random value, doesn't matter when N is fixed across studies

    sigma2  <- (1-rho^2)^2 / (sample_size -1)
    s2 <- compute_s2(sigma2, k)

    tau2 / (s2 + tau2) #I2

}

compute_tau2_z <- function(vec = c(I2, sample_size)){
    #I2 = I2 from pearson R
    #sample_size = fixed within-study sample size
    I2 <- vec[["I2"]]
    sample_size <- vec[["sample_size"]]
    k <- 20 #constant, random value, doesn't matter when N is fixed across studies

    sigma2  <- 1 / (sample_size - 3)
    s2 <- compute_s2(sigma2, k)

    I2*s2 / (1 - I2) #tau2
}

#****************************************
# Computations
#****************************************

#Conditions

rho  <-  0
tau_r <- c(0.1, 0.15, 0.2) #Tau values we use for Pearson's r
sample_size <- c(50, 100, 150, 200)

cond <- expand.grid(rho = rho,
                    tau = tau_r,
                    sample_size = sample_size)


#Compute first I2 for Pearson's r, then tau2 for Fisher's z based on that I2
cond$I2 <- apply(cond, 1, compute_I2)
cond$z_tau2 <- apply(cond, 1, compute_tau2_z)

cond$z_tau <- sqrt(cond$z_tau2)
