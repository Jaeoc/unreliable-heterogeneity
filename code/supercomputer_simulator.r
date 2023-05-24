# Project: unexplained heterogeneity
# Script purpose: run simulations on supercomputer
# code: Anton Olsson-Collentine



#****************************************
# libraries and functions
#****************************************
library(parabar) #for a progress-bar in parallel computing
# library(parallel) #comes pre-installed with R
#library(metafor) #loaded when setting up simulations

source("code/functions.r") #for rnorm_truncated and simulate_rma

#****************************************
# Conditions
#****************************************
sample_size <- c(50, 100, 150, 200)
k <- c(5, 20, 40, 200)
mu <- seq(from = 0, to = 0.6, by = 0.1)

#Based on Flake et al., average alpha was .79, SD = .13, range .17 - .87
#Based on Sanchez-Meca, mean across 5 meta-analysis was 0.767 - 0.891 and SD ranged between 0.034 - 0.133
reliability_mean <- c(0.6, 0.7, 0.8, 0.9, 1)
reliability_sd <- seq(from = 0.15, to = 0.15, by = 0.05)

effect_type  <- "r" # or r_z
meta_method  <- "Hedges" #or HS

#Using Pearson's r as the effect size.
true_tau <- c(0, 0.1, 0.15, 0.2)
true_tau2 <- true_tau^2

if(effect_type == "r_z"){
    #Using Fisher's z as the effect size (see 'compute_tau_fisher_z.r')

    true_tau_50 <- c(0.1031263, 0.1566007, 0.2123514) #For N = 50
    true_tau_100 <- c(0.1020357, 0.1549445, 0.2101057) #For N = 100
    true_tau_150 <- c(0.1016845, 0.1544113, 0.2093826) #For N = 150
    true_tau_200 <- c(0.1015112, 0.1541480, 0.2090256) #For N = 200
    true_tau <- c(0, true_tau_50, true_tau_100, true_tau_150, true_tau_200)
    true_tau2 <- true_tau^2

    mu <- atanh(mu)
}


cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_mean = reliability_mean,
                    reliability_sd = reliability_sd)

#drop cases where perfect reliability and reliability_sd > 0
remove_cond <- cond[["reliability_mean"]] == 1 & cond[["reliability_sd"]] > 0
cond <- cond[!remove_cond,]

# fisher' z
if(effect_type == "r_z"){
    #drop incorrect tau-values for fisher's z
    keep_cond <- cond[["true_tau2"]] == 0 |
                (cond[["sample_size"]] == 50 & cond[["true_tau2"]] %in% true_tau_50^2) |
                (cond[["sample_size"]] == 100 & cond[["true_tau2"]] %in% true_tau_100^2) |
                (cond[["sample_size"]] == 150 & cond[["true_tau2"]] %in% true_tau_150^2) |
                (cond[["sample_size"]] == 200 & cond[["true_tau2"]] %in% true_tau_200^2)
    cond <- cond[keep_cond,]
}

#these below are higher level control input to the function
cond$effect_type <- effect_type
cond$method  <- meta_method
cond$reliability_distribution <- "normal"
cond$step_length  <-  0.5 #decrease from 1 to 0.5 to improve convergence for low N
cond$maxiterations  <-  100 #default values are steplength = 1, and maxiter = 100
cond$reliability_max <- cond$reliability_mean
cond$reliability_min <- cond$reliability_mean

#****************************************
# Setup parallel simulation
#****************************************
reps <- 1e1
out_list <- vector("list", length = nrow(cond))

ncores <-parallel::detectCores()
cl <- parabar::start_backend(ncores) # Create cluster based on nworkers.

parabar::evaluate(cl, library(metafor))
parabar::export(cl, c("simulate_rma", "rnorm_truncated"))

#****************************************
# Run simulation
#****************************************

start <- Sys.time()
for(r in 1:nrow(cond)){ #gives us a list of lists

    progress_bar_format <- paste0(
        "Condition ", r, "/", nrow(cond), ". [:bar] :percent [:elapsed]"
    )

    parabar::configure_bar(type = "modern", format = progress_bar_format)

    task <- function(iteration, cond_r){ #anonymous function needed when using for replications
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


    out_list[[r]] <- parabar::par_lapply(
                    backend = cl, # cluster
                    x = 1:reps, #looping over
                    fun = task,
                    cond_r = cond[r, ]
    )
}
stop_backend(cl) # Shut down the nodes and end simulation

end <- Sys.time()
end - start

#****************************************
# Clean up and save results
#****************************************

e <- lapply(out_list, data.table::rbindlist)
names(e) <- c(paste0("mu = ", cond$mu,
                     ";k = ", cond$k,
                     ";N = ", cond$sample_size,
                     ";reliability_sd = ", cond$reliability_sd,
                     ";mean_rel = ", cond$reliability_mean,
                     ";true_tau2 = ", cond$true_tau2,
                     ";effect_type = ", cond$effect_type,
                     ";method = ", cond$method))


#saveRDS(e, "../data_new/raw_r_tau_0-0.2.RDS")
