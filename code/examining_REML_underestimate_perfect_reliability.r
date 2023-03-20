# Project: unexplained heterogeneity
# Script purpose: Explore the underestimate of tau2 by REML for large tau2 values


#****************************************
# Replicating figure 2 in Langan et al. (2019)
#****************************************
# langan: a comparison of heterogeneity variance estimators in simulated random-effects
#https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1316


simulate_rma_langan <- function(version = c("original", "efficient"), k = 20, true_tau2 = 0.0444, mu = 0.5, sample_size = 20){ #see section 3.1

    mu_i <- rnorm(n = k, mean = mu, sd = sqrt(true_tau2))
    n_c = n_t = rep(sample_size, k)

    if(version == "original"){
    n1 <- replicate(k, rnorm(n = sample_size, mean = 0, sd = 1), simplify = FALSE)
    n2 <- lapply(1:k, function(i) rnorm(n = sample_size, mean = mu_i[i], sd = 1))

    avg_c <- lapply(n1, mean) #function below from heterogeneity project #1
    var_c <- lapply(n1, var)
    avg_c <- do.call(rbind, avg_c)
    var_c <- do.call(rbind, var_c)

    avg_t <- lapply(n2, mean)
    var_t <- lapply(n2, var)
    avg_t <- do.call(rbind, avg_t)
    var_t <- do.call(rbind, var_t)

    fit <- rma(measure = "SMD",  m1i = avg_t, m2i = avg_c, sd1i = sqrt(var_t),
                   sd2i = sqrt(var_c), n1i = n_t, n2i = n_c, method = "REML") #fit meta-analysis transforming into standardized mean difference (Hedge's g)

    }else if(version == "efficient"){

        avg_c <- rnorm(n = k, mean = 0, sd = 1 / sqrt(n_c)) #Draw means from sampling distribution control group
        avg_t <- rnorm(n = k, mean = mu_i, sd = 1 / sqrt(n_t)) #Draw means from sampling distribution treatment group
        var_c <- rchisq(n = k, df = n_c - 1) / (n_c - 1) #draw variances from sampling distribution control group
        var_t <- rchisq(n = k, df = n_t - 1) / (n_t - 1) #draw variances from sampling distribution treatment group

        fit <- rma(measure = "SMD",  m1i = avg_t, m2i = avg_c, sd1i = sqrt(var_t),
                   sd2i = sqrt(var_c), n1i = n_t, n2i = n_c, method = "REML") #fit meta-analysis transforming into standardized mean difference (Hedge's g)

    }

        data.frame(intercept_b = fit$b,
               intercept_p = fit$pval,
               tau2_hat = fit$tau2,
               tau2_p = fit$QEp)
}

# Run simulation
k = 20
group_sample_size = 20
mu = 0.5
true_tau2 = 0.991 #see Figure 1
reps <- 1e4
#

cond <- expand.grid(sample_size = group_sample_size,
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

                        simulate_rma_langan(version = "original",
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

#Gives similar results to Langan et al.

#****************************************
# Exploring different variances
#****************************************

#using the main function simulate_rma we
# a) use rho instead of r_se_me in estimating variance
# b) change Pearson's r for Fisher's z
# c) Use Fisher's z and the true sampling_var_rho instead of its estimate
# d) Use Pearson's r and the true sampling_var_rho instead of its estimate

# c) gives correct estimates for any tau2 values, d) gives correct estimates as long as tau2 is not very large. Generally, as long as tau2 is not large (i.e., above about 0.35) any method gives correct results.

simulate_rma <- function(effect_type = c("r", "r_z"), k, sample_size, true_tau2, mu, reliability_min, reliability_max){
 # Draw rho from mean rho
    rho <- rnorm(n = k, mean = mu, sd = sqrt(true_tau2))

    if(effect_type == "r") {
        sampling_var_rho  <- (1-rho^2)^2 / (sample_size -1)

    }else if(effect_type == "r_z") { #fisher's z

        sampling_var_rho  <- 1 / (sample_size - 3)
    }else{
        stop("specify effect type")
    }

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
