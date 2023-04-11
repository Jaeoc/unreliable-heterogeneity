# Project: unexplained heterogeneity
# Script purpose: run simulations on supercomputer
# code: Anton Olsson-Collentine



#****************************************
# libraries and functions
#****************************************
library(parallel) #comes pre-installed with R
#library(metafor) #loaded when setting up simulations

rnorm_truncated <- function(n, mean, sd, lower_bound, upper_bound){
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
    fit <- rma(yi = r_se_me, vi = r_var_estimate,
              control=list(stepadj=steplength, maxiter=maxiter))

    data.frame(intercept_b = fit$b,
               intercept_p = fit$pval,
               tau2_hat = fit$tau2,
               tau2_p = fit$QEp)
}


#****************************************
# generate data for all conditions
#****************************************

# conditions
sample_size <- 150
k <- 20
true_tau2 <- c(0.002, 0.0035, 0.0055, 0.0085, 0.012, 0.0185, 0.031, 0.069)
#Using Pearson's r as the effect size
mu <- seq(from = 0, to = 0.6, by = 0.1) #

#Based on Flake et al., average alpha was .79, SD = .13, range .17 - .87
#Based on Sanchez-Meca, mean across 5 meta-analysis was 0.767 - 0.891 and SD ranged between 0.034 - 0.133
reliability_mean <- c(0.8, 0.9, 1)
reliability_sd <- seq(from = 0, to = 0.15, by = 0.05)


cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_mean = reliability_mean,
                    reliability_sd = reliability_sd)

#drop cases where perfect reliability and reliability > 0
remove_cond <- cond[["reliability_mean"]] == 1 & cond[["reliability_sd"]] > 0
cond <- cond[!remove_cond,]

#these below are just appendices to makes sure the function runs as expected
cond$reliability_max <- cond$reliability_mean
cond$reliability_min <- cond$reliability_mean
cond$effect_type <- "r"
cond$reliability_distribution <- "normal"


reps <- 1e4

out_list <- vector("list", length = nrow(cond))

start <- Sys.time()

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers.

clusterEvalQ(cl, library(metafor))
clusterExport(cl, c("simulate_rma", "rnorm_truncated"))

for(r in 1:nrow(cond)){ #gives us a list of lists

    mes <- paste0("\n now on condition ", r)
    cat(mes)

cond_r <- cond[r,]
#need to add because otherwise the nodes don't find r. Since this happens outside the nodes.

clusterExport(cl, "cond_r") #exported everything in environment to each node, otherwise they don't have access to all functions in the master global environment

    out_list[[r]] <- parallel::parLapply(
                    cl = cl, # cluster
                    1:reps, #looping over
                    function(iteration){ #anonymous function needed when using for replications

                        simulate_rma(effect_type = cond_r$effect_type,
                                    reliability_distribution = cond_r$reliability_distribution,
                                    k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu,
                                    reliability_min = cond_r$reliability_min,
                                    reliability_max = cond_r$reliability_max,
                                    reliability_mean = cond_r$reliability_mean,
                                    reliability_sd = cond_r$reliability_sd)
        }
    )
}
stopCluster(cl) # Shut down the nodes

e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("mu = ", cond$mu, ";reliability_sd = ", cond$reliability_sd,
                           ";mean_rel = ", cond$reliability_mean, ";true_tau2 = ", cond$true_tau2))

saveRDS(e_means, "../data_new/over_vs_underestimate.RDS")
#lapply(e_means, round, 3)

end <- Sys.time()
end - start
