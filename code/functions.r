# Project: unexplained heterogeneity
# Script purpose: functions for simulations
# code: Anton Olsson-Collentine

#****************************************
# Functions
#****************************************

rnorm_truncated <- function(
    n,
    mean,
    sd,
    lower_bound,
    upper_bound){

    # Output is a the same as from rnorm but from a truncated distribution

    #Draw probababilities from a uniform distribution with the boundaries based on the relevant normal distribution
    trunc_norm <- runif(n = n,
                        min = pnorm(q = lower_bound, mean = mean, sd = sd),
                        max = pnorm(q = upper_bound, mean = mean, sd = sd)
    )
    # Transform probabilities to x-values based on the normal distribution
    qnorm(trunc_norm, mean = mean, sd = sd)
}

simulate_rma <- function(

    effect_type = c("r", "r_z"),
    method = c("Hedges", "HS"),
    reliability_distribution = c("uniform", "normal"),
    k, #number of studies
    sample_size, #sample size, fixed across studies
    true_tau2, #variance of superpopulation
    mu, #mean of superpopulation
    reliability_min, #for uniform distribution
    reliability_max,
    reliability_mean = NULL, #for normal distribution
    reliability_sd = NULL,
    steplength = 1, #these are for controlling the fisher algorithm in rma
    maxiter = 100){ #these are the default values in the function: https://www.metafor-project.org/doku.php/tips:convergence_problems_rma

    # Output is a dataframe with selected results from a metafor::rma object


    # Draw rho from mean rho and compute sampling variances
    if(effect_type == "r") {

        rho <- rnorm_truncated(n = k, mean = mu, sd = sqrt(true_tau2),
                              lower_bound = -1, upper_bound = 1)
        sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)

        # given study rho, draw sample rho
        r_se <- rnorm_truncated(k, mean = rho, sd = sqrt(sampling_var_rho),
                                lower_bound = -1, upper_bound = 1)

    }else if(effect_type == "r_z") { #fisher's z

        rho <- rnorm(n = k, mean = mu, sd = sqrt(true_tau2))
        sampling_var_rho  <- 1 / (sample_size - 3)

        # given study rho, draw sample rho
        r_se <- rnorm(k, mean = rho, sd = sqrt(sampling_var_rho))

    }else{
        stop("specify effect type")
    }

    #draw reliability of each study k
    if(reliability_distribution == "uniform"){

        reliabilities <- runif(k, min = reliability_min, max = reliability_max)

    } else if(reliability_distribution == "normal"){

        reliabilities <- rnorm_truncated(n = k,
         mean = reliability_mean,
         sd = reliability_sd,
         lower_bound = 0,
         upper_bound = 1)
    } #we sample from a truncated normal distribution, inducing some bias in the range of reliabilities, at least on the upper limit

    #Compute observed r given measurement error
    ##If Fisher's z, transform to r before adding measurement error
    if(effect_type == "r_z"){
        r_se <- tanh(r_se)
    }

    ##Add the measurement error
    ## Assume equal reliability for X and Y, then
    ## r_xy = rho_xy * sqrt(R_xx')*sqrt(R_yy') = rho_xy * R_xx'
    r_se_me <- r_se*reliabilities

    if(effect_type == "r_z"){
        #Transform back to Fisher's z
        r_se_me <- atanh(r_se_me)
    }

    #compute observed sampling variance estimate given measurement error
    if(effect_type == "r") {
        r_var_estimate <- (1-r_se_me^2)^2 / (sample_size -1)
    }else if(effect_type == "r_z") { #fisher's z
        r_var_estimate  <- 1 / (sample_size - 3) #= sampling_var_rho
    }

    # fit meta-analysis on observed values
    if(method == "Hedges"){
        fit <- rma(yi = r_se_me, vi = r_var_estimate,
                  control=list(stepadj=steplength, maxiter=maxiter))

    else if(method == "HS"){ #https://www.metafor-project.org/doku.php/tips:hunter_schmidt_method
        r_var_estimate <- escalc(measure="COR", ri=r_se_me,
                      ni=rep(sample_size, k), vtype="AV")[["vi"]]

        fit <- rma(yi = r_se_me, vi = r_var_estimate, weights = rep(sample_size, k),
              control=list(stepadj=steplength, maxiter=maxiter))

    }


    data.frame(intercept_b = fit$b,
               intercept_p = fit$pval,
               tau2_hat = fit$tau2,
               tau2_p = fit$QEp)
}


# The below function is an adapted version of above, specifically for
# estimating tau2 values that correspond to I2 levels
# Hence, it assumes perfect reliability has adapted output

simulate_rma_I2 <- function(effect_type = c("r", "r_z"), k, sample_size, true_tau2, mu){
    #adjusted input and output compared to my primary function
    #Also, only for perfect reliability, i.e., uses sampling_var_rho to estimate meta-analysis

 # Draw rho from mean rho and compute sampling variances

    if(effect_type == "r") {
        rho <- rnorm_truncated(n = k, mean = mu, sd = sqrt(true_tau2),
                              lower_bound = -1, upper_bound = 1)
        sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)

        # given study rho, draw sample rho
        r_se <- rnorm_truncated(k, mean = rho, sd = sqrt(sampling_var_rho),
                                lower_bound = -1, upper_bound = 1)

    }else if(effect_type == "r_z") { #fisher's z
        rho <- rnorm(n = k, mean = mu, sd = sqrt(true_tau2))
        sampling_var_rho  <- 1 / (sample_size - 3)

    # given study rho, draw sample rho
    r_se <- rnorm(k, mean = rho, sd = sqrt(sampling_var_rho))

    }else{
        stop("specify effect type")
    }


    fit <- rma(yi = r_se, vi = sampling_var_rho)

    data.frame(I2 = fit$I2,
               true_tau2 = true_tau2,
               tau2_hat = fit$tau2,
               tau2_p = fit$QEp)
}
