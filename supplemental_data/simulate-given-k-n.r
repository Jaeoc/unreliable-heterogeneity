# Project: unexplained heterogeneity
# script purpose: Simulate tau2 values corresponding to an I2 of 90 - 95% for sample size 150 and number of studies equal to 20
# code: Anton Olsson-Collentine

#****************************************
# functions and packages
#****************************************


library(parallel)
#library(metafor) #is loaded through the simulation setup below

source("code/functions.r") #for rnorm_truncated and simulate_rma_I2
#****************************************
# run simulation
#****************************************

#NB! sample size and K
sample_size <- 150
k <- 20
# To try to find low tau2 values for Pearson's r.
true_tau2 <- seq(from = 0, to = 0.14, by = 0.0005)
mu <- 0 #based on my simulations mu does not matter. Set to zero to avoid truncation.
reps <- 1e4

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu)


out_list <- vector("list", length = nrow(cond))



start <- Sys.time()

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers. This because mclapply doesn't work on windows

clusterEvalQ(cl, library(metafor))
clusterExport(cl, c("simulate_rma_I2", "rnorm_truncated"))

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

                        simulate_rma_I2(effect_type = "r",
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
f_means <- do.call(rbind, e_means)
