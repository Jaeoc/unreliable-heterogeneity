# Project: unexplained heterogeneity
# Script purpose: run simulations on supercomputer
# code: Anton Olsson-Collentine



#****************************************
# libraries and functions
#****************************************
library(parallel) #comes pre-installed with R
#library(metafor) #loaded when setting up simulations

source("code/functions.r") #for rnorm_truncated and simulate_rma

#****************************************
# Conditions
#****************************************
sample_size <- c(50, 100, 150, 200)
k <- c(5, 20, 40, 200)
true_tau <- c(0, 0.1, 0.15, 0.2)
true_tau2 <- true_tau^2
#Using Pearson's r as the effect size. Correspond to I2 20 - 90%.
# Below corresponding values for k = 5
# k <- 5
#true_tau2 <- c(0, 0.0015, 0.0045, 0.008, 0.013, 0.0195, 0.0315, 0.0550, 0.1245)
mu <- seq(from = 0, to = 0.6, by = 0.1) #

#Based on Flake et al., average alpha was .79, SD = .13, range .17 - .87
#Based on Sanchez-Meca, mean across 5 meta-analysis was 0.767 - 0.891 and SD ranged between 0.034 - 0.133
reliability_mean <- c(0.6, 0.7, 0.8, 0.9)
reliability_sd <- seq(from = 0.15, to = 0.15, by = 0.05)


cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_mean = reliability_mean,
                    reliability_sd = reliability_sd)

#drop cases where perfect reliability and reliability > 0
remove_cond <- cond[["reliability_mean"]] == 1 & cond[["reliability_sd"]] > 0
cond <- cond[!remove_cond,]

#these below are higher level control input to the function
cond$reliability_max <- cond$reliability_mean
cond$reliability_min <- cond$reliability_mean
cond$effect_type <- "r"
cond$reliability_distribution <- "normal"
cond$step_length  <-  0.5 #decrease from 1 to 0.5 to improve convergence for low N
cond$maxiterations  <-  100 #default values are steplength = 1, and maxiter = 100

#****************************************
# Setup parallel simulation
#****************************************
reps <- 1e2
out_list <- vector("list", length = nrow(cond))

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers.

clusterEvalQ(cl, library(metafor))
clusterExport(cl, c("simulate_rma", "rnorm_truncated"))

#****************************************
# Run simulation
#****************************************
start <- Sys.time()
for(r in 1:nrow(cond)){ #gives us a list of lists

    mes <- paste0("\n now on condition ", r)
    cat(mes)

cond_r <- cond[r,]
#need to add because otherwise the nodes don't find r. Since this happens outside the nodes.

clusterExport(cl, "cond_r") #exported each row to the nodes otherwise they don't have access

task <- function(iteration){ #anonymous function needed when using for replications
        simulate_rma(effect_type = cond_r$effect_type,
                    reliability_distribution = cond_r$reliability_distribution,
                    k = cond_r$k,
                    sample_size = cond_r$sample_size,
                    true_tau2 = cond_r$true_tau2,
                    mu = cond_r$mu,
                    reliability_min = cond_r$reliability_min,
                    reliability_max = cond_r$reliability_max,
                    reliability_mean = cond_r$reliability_mean,
                    reliability_sd = cond_r$reliability_sd,
                    steplength = cond_r$step_length,
                    maxiter = cond_r$maxiterations)
    }


    out_list[[r]] <- parallel::parLapply(
                    cl = cl, # cluster
                    1:reps, #looping over
                    task
    )
}
stopCluster(cl) # Shut down the nodes and end simulation

end <- Sys.time()
end - start

#****************************************
# Clean up and save results
#****************************************

e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("mu = ", cond$mu, ";reliability_sd = ", cond$reliability_sd,
                           ";mean_rel = ", cond$reliability_mean, ";true_tau2 = ", cond$true_tau2))

lapply(e_means, round, 3)
#saveRDS(e_means, "../data_new/over_vs_underestimate.RDS")
