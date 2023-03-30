# Project: unexplained heterogeneity
# script purpose: Simulate tau2 values corresponding to an I2 of 90 - 95% for sample size 150 and number of studies equal to 20
# code: Anton Olsson-Collentine

#****************************************
# functions and packages
#****************************************

simulate_rma_I2 <- function(effect_type = c("r", "r_z"), k, sample_size, true_tau2, mu){
    #adjusted input and output compared to my primary function
    #Also, only for perfect reliability, i.e., uses sampling_var_rho to estimate meta-analysis

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

    fit <- rma(yi = r_se, vi = sampling_var_rho)

    data.frame(I2 = fit$I2,
               true_tau2 = true_tau2,
               tau2_hat = fit$tau2,
               tau2_p = fit$QEp)
}

library(parallel)
#library(metafor) #is loaded through the simulation setup below

#****************************************
# run simulation
#****************************************

#NB! sample size and K
sample_size <- 150
k <- 20
#true_tau2 <- seq(from = 0.05, to = 0.2, by = 0.001)
#Initial experimentation suggests 90% I2 is around 0.8 and 95% around 0.15
#EDIT: for Fisher's z the above tau2 values are ok, for Pearson's r they need to be lower, about 0.03 to 0.1. I have here used mu = 0.3/
# To try to find low tau2 values for Pearson's r I update these values to  the below. I also set mu = 0.
true_tau2 <- seq(from = 0, to = 0.05, by = 0.001)
mu <- 0
reps <- 1e3

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu)


out_list <- vector("list", length = nrow(cond))



start <- Sys.time()

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers. This because mclapply doesn't work on windows

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

                        simulate_rma_I2(effect_type = "r_z",
                                    k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu)
        }
    )
}
stopCluster(cl) # Shut down the nodes

end <- Sys.time()
end - start


e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
lapply(e_means, round, 3)

#saveRDS("I2_high_fisherz.RDS")
#saveRDS("I2_high_r.RDS")
