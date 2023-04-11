# Project: unexplained heterogeneity
# Code purpose: demonstrate how to add a progress bar in parallel computing using the 'parabar' package


sample_size <- 150
k <- 20
true_tau2 <- c(0.002, 0.0035, 0.0055, 0.0085, 0.012, 0.0185, 0.031, 0.069) #corresponds to 20 -90% I2 according to my simulations for Pearson's r and the above sample size and k
mu <- seq(from = 0, to = 0.6, by = 0.1)
reliability_min <- seq(from = 1, to = 1, by = 0.1)
reps <- 1e4

cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_min = reliability_min)

cond$reliability_max <- cond$reliability_min

out_list <- vector("list", length = nrow(cond))


library(parabar)

start <- Sys.time()

ncores <- detectCores()
# cl <- makePSOCKcluster(ncores) # Create cluster based on nworkers.
cl <- start_backend(ncores)

# clusterEvalQ(cl, library(metafor))
# clusterExport(cl, c("simulate_rma", "rnorm_truncated"))

evaluate(cl, library(metafor))
export(cl, c("simulate_rma", "rnorm_truncated"))
peek(cl)

for(r in 1:nrow(cond)){ #gives us a list of lists

    progress_bar_format <- paste0(
        "Condition ", r, "/", nrow(cond), ". [:bar] :percent [:elapsed]"
    )

    configure_bar(type = "modern", format = progress_bar_format)

    # mes <- paste0("\n now on condition ", r)
    # cat(mes)

    cond_r <- cond[r,]
    #need to add because otherwise the nodes don't find r. Since this happens outside the nodes.

    # clusterExport(cl, "cond_r") #exported row to each node
    export(cl, "cond_r")

    task <- function(iteration){ #anonymous function needed when using for replications
        simulate_rma(effect_type = "r",
            reliability_distribution = "uniform",
            k = cond_r$k,
            sample_size = cond_r$sample_size,
            true_tau2 = cond_r$true_tau2,
            mu = cond_r$mu,
            reliability_min = cond_r$reliability_min,
            reliability_max = cond_r$reliability_max)
    }

    # out_list[[r]] <- parallel::parLapply(
    #     cl = cl, # cluster
    #     1:reps, #looping over
    #     task
    # )

    out_list[[r]] <- parabar::par_sapply(
        backend = cl, # cluster
        x = 1:reps, #looping over
        fun = task
    )

}

# stopCluster(cl) # Shut down the nodes
stop_backend(cl)

e <- lapply(out_list, function(x) do.call(rbind, x))
e_means <- lapply(e, colMeans)
names(e_means) <- paste0("mu = ", cond$mu, ";true_tau2 = ", cond$true_tau2)

#saveRDS(e_means, "../data_new/reliability_1,.RDS")
#lapply(e_means, round, 3)

end <- Sys.time()
end - start
