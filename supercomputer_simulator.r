# Project: unexplained heterogeneity
# Script purpose: run simulations on supercomputer
# code: Anton Olsson-Collentine



#****************************************
# libraries and functions
#****************************************
library(parallel) #comes pre-installed with R

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
# 1) Mu and varying reliability, to show lack of effect
#****************************************
sample_size <- 150
k <- 20
true_tau2 <- 0.1
mu <- seq(from = 0, to = 0.9, by = 0.1)
reliability_min <- 0
reliability_max <- 1
reps <- 1e4

cond_mu <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_min = reliability_min,
                    reliability_max = reliability_max)


reliability_min <- 0.5
reliability_max <- 0.5

cond_mu2 <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_min = reliability_min,
                    reliability_max = reliability_max)



mu  <- 0.3
reliability_min <- seq(from = 0.5, to = 0, by = -0.1)
reliability_max <- seq(from = 0.5, to = 1, by = 0.1)

cond_rel <- data.frame(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_min = reliability_min,
                    reliability_max = reliability_max)


cond <- rbind(cond_mu, cond_mu2, cond_rel) #my initial conditions



out_list <- vector("list", length = nrow(cond))



start <- Sys.time()

ncores <- detectCores()
cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers. This because mclapply doesn't work on windows

clusterEvalQ(cl, library(metafor))

for(r in 1:nrow(cond)){ #gives us a list of lists

cond_r <- cond[r,]
#need to add because otherwise the nodes don't find r. Since this happens outside the nodes.

clusterExport(cl, ls(.GlobalEnv)) #exported everything in environment to each node, otherwise they don't have access to all functions in the master global environment

    out_list[[r]] <- parallel::parLapply(
                    cl = cl, # cluster
                    1:reps, #looping over
                    function(iteration){ #anonymous function needed when using for replications

                        simulate_rma(k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu,
                                    reliability_min = cond_r$reliability_min,
                                    reliability_max = cond_r$reliability_max)
        }
    )
}
stopCluster(cl) # Shut down the nodes

end <- Sys.time()
end - start


e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("mu = ", cond$mu[1:10], "& Reliability = 0 - 1"),
 paste0("mu = ", cond$mu[11:20], "& Reliability = 0.5"),
 paste0("reliability = ", reliability_min, " - ", reliability_max))

saveRDS(e_means, "mu_varying_rel.RDS")
#lapply(e_means, round, 3)


#****************************************
# 2) tau2 and average reliability
#****************************************

# Tau2

sample_size <- 150
k <- 20
true_tau2 <- seq(from = 0, to = 0.9, by = 0.1)
mu <- 0.3
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

                        simulate_rma(k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu,
                                    reliability_min = cond_r$reliability_min,
                                    reliability_max = cond_r$reliability_max)
        }
    )
}
stopCluster(cl) # Shut down the nodes

end <- Sys.time()
end - start


e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("tau2 = ", cond$true_tau2, " & reliability = ", cond$reliability_min))

saveRDS(e_means, "tau2_increasing_rel.RDS")
#lapply(e_means, round, 3)


#****************************************
# 1) Mu and varying reliability, to show effect...
#****************************************
sample_size <- 150
k <- 20
true_tau2 <- 0.1
mu <- seq(from = 0, to = 0.9, by = 0.1)
reliability_min <- seq(from = 0.9, to = 0, by = -0.1)
reliability_max <- 1
reps <- 1e3

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_min = reliability_min,
                    reliability_max = reliability_max)



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

                        simulate_rma(k = cond_r$k,
                                    sample_size = cond_r$sample_size,
                                    true_tau2 = cond_r$true_tau2,
                                    mu = cond_r$mu,
                                    reliability_min = cond_r$reliability_min,
                                    reliability_max = cond_r$reliability_max)
        }
    )
}
stopCluster(cl) # Shut down the nodes

end <- Sys.time()
end - start


e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- c(paste0("mu = ", cond$mu, " & min reliability = ", cond$reliability_min))

saveRDS(e_means, "mu_rel_variability.RDS")
#lapply(e_means, round, 3)
