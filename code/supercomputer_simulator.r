# Project: unexplained heterogeneity
# Script purpose: run simulations on supercomputer
# code: Anton Olsson-Collentine



#****************************************
# libraries and functions
#****************************************
library(parabar) #for a progress-bar in parallel computing
library(data.table) #to bind lists into dataframe efficiently
# I load these here because using :: slows down the simulation substantially
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
reliability_sd <- seq(from = 0.05, to = 0.15, by = 0.05)

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
cond$step_length  <-  0.5 #decrease from 1 to 0.5 to improve convergence for low N
cond$maxiterations  <-  100 #default values are steplength = 1, and maxiter = 100

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

    progress_bar_format <- paste0( #for parabar
        "Condition ", r, "/", nrow(cond), ". [:bar] :percent [:elapsed]"
    )

    configure_bar(type = "modern", format = progress_bar_format) #from parabar

    task <- function(iteration, cond_r){ #anonymous function needed when using for replications
                simulate_rma(effect_type = cond_r$effect_type,
                            method = cond_r$method,
                            k = cond_r$k,
                            sample_size = cond_r$sample_size,
                            true_tau2 = cond_r$true_tau2,
                            mu = cond_r$mu,
                            reliability_mean = cond_r$reliability_mean,
                            reliability_sd = cond_r$reliability_sd,
                            steplength = cond_r$step_length,
                            maxiter = cond_r$maxiterations)
        }


    out_list[[r]] <- par_lapply( #function from parabar
                    backend = cl, # cluster
                    x = 1:reps, #looping over
                    fun = task,
                    cond_r = cond[r, ]
    )

    #Compute the mean across replications for the condition and add condition identifiers
    out_list[[r]] <- rbindlist(out_list[[r]]) #function from data.table
    out_list[[r]] <- out_list[[r]][, lapply(.SD, mean)] #data.table colMeans but returns a dataframe (well, data.table)
    out_list[[r]] <- cbind(out_list[[r]], cond[r,])


}

stop_backend(cl) # Shut down the nodes and end simulation

end <- Sys.time()
end - start


#****************************************
# Clean up and save results
#****************************************

e <- rbindlist(out_list)

#saveRDS(e, "../data_new/raw_r_tau_0-0.2.RDS")

    if(r %% 1000 == 0 | r == nrow(cond)){ #if even thousand or simulation finished

        e <- rbindlist(out_list[r-999:r])

        file_name <- paste0("../data/simulation_means_cond_", r-999,"-", r, ".csv")
        fwrite(x = e, file = file_name)

        rm(e)
    }
