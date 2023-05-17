# Project: unexplained heterogeneity
# Script purpose: functions for simulations
# code: Anton Olsson-Collentine

#****************************************
# Simulation functions
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

    }else if(method == "HS"){ #https://www.metafor-project.org/doku.php/tips:hunter_schmidt_method
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



#****************************************
# Plot prep functions
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

# Compute variance of a truncated normal distribution
compute_var_truncated <- function(mu, nominal_tau, a = -1, b = 1){
# Formula from wikipedia (double-check): https://en.wikipedia.org/wiki/Truncated_normal_distribution

    #mu = mean of distribution before truncation
    #nominal_tau = standard deviation of distribution before truncation
    #a  = lower bound of truncation
    #b = upper bound of truncation

    if(nominal_tau == 0){
        return(rep(0, length(mu)))
    }

    alpha <- (a - mu) / nominal_tau
    beta <-  (b - mu) / nominal_tau

    denom <- pnorm(beta) - pnorm(alpha)
    num1 <- beta*dnorm(beta) - alpha*dnorm(alpha)
    num2 <- dnorm(beta) - dnorm(alpha)
    true_tau2 <- nominal_tau^2 * (1 - num1/denom - (num2 / denom)^2)

    sqrt(true_tau2) #out true_tau

}
