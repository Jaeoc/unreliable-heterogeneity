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
    k,
    sample_size,
    true_tau2,
    mu,
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

    }else if(effect_type == "r_z") { #fisher's z

        rho <- rnorm(n = k, mean = mu, sd = sqrt(true_tau2))
        sampling_var_rho  <- 1 / (sample_size - 3)

    }else{
        stop("specify effect type")
    }

    # given study rho, draw sample rho
    r_se <- rnorm(k, mean = rho, sd = sqrt(sampling_var_rho))

    #draw new reliabilities for each run to avoid some strange draw
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
    #If Fisher's z transform to r before adding measurement error
    if(effect_type == "r_z"){
        r_se <- tanh(r_se)
    }

    r_se_me <- r_se*sqrt(reliabilities)^2

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
# Plot 1) effect size and variance in reliability estimates increase----
#****************************************
# For a fixed effects model this leads to an overestimate of reliability

sample_size <- 150
k <- 20
true_tau2 <- 0 #Note, fixed effect model
mu <- seq(from = 0, to = 0.9, by = 0.1)
mu <- atanh(mu)

#Based on Flake et al., average alpha was .79, SD = .13, range .17 - .87
#Based on Sanchez-Meca, mean across 5 meta-analysis was 0.767 - 0.891 and SD ranged between 0.034 - 0.133
reliability_mean <- 0.8
reliability_sd <- seq(from = 0, to = 0.15, by = 0.05)
reps <- 1e4

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_mean = reliability_mean,
                    reliability_sd = reliability_sd)

#these below are just appendices to makes sure the function runs as expected
cond$reliability_max <- cond$reliability_mean
cond$reliability_min <- cond$reliability_mean


out_list <- vector("list", length = nrow(cond))



start <- Sys.time()

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers.

clusterEvalQ(cl, library(metafor))

for(r in 1:nrow(cond)){ #gives us a list of lists

    mes <- paste0("\n now on condition ", r)
    cat(mes)

cond_r <- cond[r,]
#need to add because otherwise the nodes don't find r. Since this happens outside the nodes.

clusterExport(cl, ls(.GlobalEnv)) #exported everything in environment to each node, otherwise they don't have access to all functions in the master global environment

    out_list[[r]] <- parallel::parLapply(
                    cl = cl, # cluster
                    1:reps, #looping over
                    function(iteration){ #anonymous function needed when using for replications

                        simulate_rma(effect_type = "r_z",
                                    reliability_distribution = "normal",
                                    k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu,
                                    reliability_min = cond_r$reliability_min,
                                    reliability_max = cond_r$reliability_max,
                                    reliability_mean = cond_r$reliability_mean,
                                    reliability_sd = cond_r$reliability_sd,
                                    steplength = 0.5,
                                    maxiter = 1e3)
        }
    )
}
stopCluster(cl) # Shut down the nodes


e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("mu = ", cond$mu, " & reliability_sd = ", cond$reliability_sd))

#saveRDS(e_means, "reliability_overestimate.RDS")
lapply(e_means, round, 3)

end <- Sys.time()
end - start

#****************************************
# Plot 2) as reliability decreases the underestimate becomes worse----
#****************************************
#Reliability is fixed in each study here

sample_size <- 150
k <- 20
true_tau2 <- 0.078 #corresponds to 95% I2 according to my simulations for Pearson's r and the above sample size and k
#true_tau2 <- 0.069 #corresponds to 90% I2 according to my simulations for Fisher's z and the above sample size and k
mu <- seq(from = 0, to = 0.6, by = 0.1)
#Based on the 'upper median' (83.35% quantile) from Schäfer & Schwarz for non-preregistered studies.
reliability_min <- seq(from = 1, to = 1, by = 0.1)
reps <- 1e2

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_min = reliability_min)

cond$reliability_max <- cond$reliability_min

out_list <- vector("list", length = nrow(cond))



start <- Sys.time()

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers.

clusterEvalQ(cl, library(metafor))

for(r in 1:nrow(cond)){ #gives us a list of lists

    mes <- paste0("\n now on condition ", r)
    cat(mes)

cond_r <- cond[r,]
#need to add because otherwise the nodes don't find r. Since this happens outside the nodes.

clusterExport(cl, ls(.GlobalEnv)) #exported everything in environment to each node, otherwise they don't have access to all functions in the master global environment

    out_list[[r]] <- parallel::parLapply(
                    cl = cl, # cluster
                    1:reps, #looping over
                    function(iteration){ #anonymous function needed when using for replications

                        simulate_rma(effect_type = "r",
                                    reliability_distribution = "uniform",
                                    k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu,
                                    reliability_min = cond_r$reliability_min,
                                    reliability_max = cond_r$reliability_max)
        }
    )
}
stopCluster(cl) # Shut down the nodes


e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("mu = ", cond$mu, " & reliability = ", cond$reliability_min))

#saveRDS(e_means, "reliability_underestimate.RDS")
#lapply(e_means, round, 3)

end <- Sys.time()
end - start

#****************************************
# Plot 3) As true tau2 increases, the underestimate becomes worse
#****************************************
# This is only in the absolute sense, relatively stays the same


sample_size <- 150
k <- 20
true_tau2 <- seq(from = 0, to = 0.27, by = 0.01)
mu <- 0.24 #median effect size correlational studies (Michèles meta-meta-analysis)
#Using Pearson's r as the effect size
true_tau2 <- atanh(true_tau2) #transformed to fisher's z
mu <- atanh(mu) #transformed to fisher's z
reliability_min <- seq(from = 0.1, to = 1, by = 0.1)
reps <- 1e4

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_min = reliability_min)

cond$reliability_max <- cond$reliability_min



out_list <- vector("list", length = nrow(cond))



start <- Sys.time()

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers.

clusterEvalQ(cl, library(metafor))

for(r in 1:nrow(cond)){ #gives us a list of lists

    mes <- paste0("\n now on condition ", r)
    cat(mes)

cond_r <- cond[r,]
#need to add because otherwise the nodes don't find r. Since this happens outside the nodes.

clusterExport(cl, ls(.GlobalEnv)) #exported everything in environment to each node, otherwise they don't have access to all functions in the master global environment

    out_list[[r]] <- parallel::parLapply(
                    cl = cl, # cluster
                    1:reps, #looping over
                    function(iteration){ #anonymous function needed when using for replications

                        simulate_rma(effect_type = "r_z",
                                    k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu,
                                    reliability_min = cond_r$reliability_min,
                                    reliability_max = cond_r$reliability_max)
        }
    )
}
stopCluster(cl) # Shut down the nodes


e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("tau2 = ", cond$true_tau2, " & reliability = ", cond$reliability_min))

#saveRDS(e_means, "tau2_increasing_rel.RDS")
lapply(e_means, round, 3)

end <- Sys.time()
end - start


#****************************************
# Plot 4) overestimate is compensated by the underestimate
#****************************************


sample_size <- 150
k <- c(5, 20)
true_tau2 <- 0.069 #90% I2 in Fisher's z given above N and K
#Using Pearson's r as the effect size
mu <- seq(from = 0, to = 0.9, by = 0.1) #
mu <- atanh(mu) #Fisher's z

#Based on Flake et al., average alpha was .79, SD = .13, range .17 - .87
#Based on Sanchez-Meca, mean across 5 meta-analysis was 0.767 - 0.891 and SD ranged between 0.034 - 0.133
reliability_mean <- c(0.8, 0.9)
reliability_sd <- seq(from = 0, to = 0.15, by = 0.05)
reps <- 1e1

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_mean = reliability_mean,
                    reliability_sd = reliability_sd)

#these below are just appendices to makes sure the function runs as expected
cond$reliability_max <- cond$reliability_mean
cond$reliability_min <- cond$reliability_mean


out_list <- vector("list", length = nrow(cond))



start <- Sys.time()

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers.

clusterEvalQ(cl, library(metafor))

for(r in 1:nrow(cond)){ #gives us a list of lists

    mes <- paste0("\n now on condition ", r)
    cat(mes)

cond_r <- cond[r,]
#need to add because otherwise the nodes don't find r. Since this happens outside the nodes.

clusterExport(cl, ls(.GlobalEnv)) #exported everything in environment to each node, otherwise they don't have access to all functions in the master global environment

    out_list[[r]] <- parallel::parLapply(
                    cl = cl, # cluster
                    1:reps, #looping over
                    function(iteration){ #anonymous function needed when using for replications

                        simulate_rma(effect_type = "r_z",
                                    reliability_distribution = "normal",
                                    k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu,
                                    reliability_min = cond_r$reliability_min,
                                    reliability_max = cond_r$reliability_max,
                                    reliability_mean = cond_r$reliability_mean,
                                    reliability_sd = cond_r$reliability_sd,
                                    steplength = 0.5,
                                    maxiter = 1e3)
        }
    )
}
stopCluster(cl) # Shut down the nodes


e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("mu = ", cond$mu, ";reliability_sd = ", cond$reliability_sd, ";k = ", cond$k, ";mean_rel = ", cond$reliability_mean))

saveRDS(e_means, "./data/test_over_vs_underestimate.RDS")
#lapply(e_means, round, 3)

end <- Sys.time()
end - start
