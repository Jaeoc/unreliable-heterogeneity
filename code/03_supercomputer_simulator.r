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


source("code/02_functions.r") #for rnorm_truncated and simulate_rma

#****************************************
# Full conditions
#****************************************
sample_size <- c(50, 100, 150, 200)
k <- c(5, 20, 40, 200)
reliability_mean <- c(0.6, 0.7, 0.8, 0.9, 1)
reliability_sd <- seq(from = 0, to = 0.15, by = 0.05)

mu <- seq(from = 0, to = 0.6, by = 0.1)

#Using Pearson's r as the effect size.
effect_type  <- "r" #
r_method  <- c("HV", "HS") #Hedges & Vevea and Hunter & Schmidt

true_tau <- c(0, 0.1, 0.15, 0.2)
true_tau2 <- true_tau^2

r_cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_mean = reliability_mean,
                    reliability_sd = reliability_sd,
                    effect_type = effect_type,
                    method = r_method)

# Using Fisher' z as the effect size
effect_type  <- "r_z"
z_method  <- "HV"

mu <- atanh(mu)

#Tau for Fisher's z (see 'compute_tau_fisher_z.r')
true_tau_50 <- c(0.1031263, 0.1566007, 0.2123514) #For N = 50
true_tau_100 <- c(0.1020357, 0.1549445, 0.2101057) #For N = 100
true_tau_150 <- c(0.1016845, 0.1544113, 0.2093826) #For N = 150
true_tau_200 <- c(0.1015112, 0.1541480, 0.2090256) #For N = 200
true_tau <- c(0, true_tau_50, true_tau_100, true_tau_150, true_tau_200)
true_tau2 <- true_tau^2

z_cond <- expand.grid(sample_size = sample_size,
                    k = k,
                    true_tau2 = true_tau2,
                    mu = mu,
                    reliability_mean = reliability_mean,
                    reliability_sd = reliability_sd,
                    effect_type = effect_type,
                    method = z_method)

#drop incorrect tau-values for fisher's z
keep_cond <- z_cond[["true_tau2"]] == 0 |
            (z_cond[["sample_size"]] == 50 & z_cond[["true_tau2"]] %in% true_tau_50^2) |
            (z_cond[["sample_size"]] == 100 & z_cond[["true_tau2"]] %in% true_tau_100^2) |
            (z_cond[["sample_size"]] == 150 & z_cond[["true_tau2"]] %in% true_tau_150^2) |
            (z_cond[["sample_size"]] == 200 & z_cond[["true_tau2"]] %in% true_tau_200^2)
z_cond <- z_cond[keep_cond,]

#combine Pearson's r and Fisher'z conditions into one
cond <- rbind(r_cond, z_cond)

#drop cases where perfect reliability and reliability_sd > 0
remove_cond <- cond[["reliability_mean"]] == 1 & cond[["reliability_sd"]] > 0
cond <- cond[!remove_cond,]

#Higher level control input to the function
cond$step_length  <-  0.5 #decrease from 1 to 0.5 to improve convergence for low N
cond$maxiterations  <-  100 #default values are steplength = 1, and maxiter = 100

#****************************************
# Subset conditions
#****************************************
#1) Figure 1 main manuscript
cond1 <- cond[cond$sample_size == 150 &
              cond$k == 20 &
              cond$reliability_sd ==0.15,]

#2) Figure 2 main manuscript
true_tau <- c(0.02, 0.04, 0.06, 0.08) #NB! different tau2!
true_tau2 <- true_tau^2

cond2 <- cond[cond$sample_size == 150 &
              cond$reliability_sd ==0.15 &
              cond$method == "HV" &
              cond$effect_type == "r",]
cond2$true_tau2 <- rep(true_tau2, each = 4)

# Supplement A1 (variance in reliabilities)
cond3 <- cond[cond$k == 20 &
              cond$sample_size == 150 &
              cond$reliability_mean < 1 &
              cond$method == "HV" &
              cond$effect_type == "r",]

# Supplement A2 (variable K)
cond4 <- cond[cond$sample_size == 150 &
              cond$reliability_sd == 0.15 &
              cond$method == "HV" &
              cond$effect_type == "r",]

# Supplement A3 (variable sample size)
cond5 <- cond[cond$k == 20 &
              cond$reliability_sd == 0.15 &
              cond$method == "HV" &
              cond$effect_type == "r",]

# Supplement B (Perfect reliability)
cond6 <- cond[cond$k == 20 &
              cond$sample_size == 150 &
              cond$reliability_mean == 1,]


#****************************************
# Conditions run
#****************************************
cond <- rbind(cond1, cond2, cond3, cond4, cond5, cond6)

cond <- cond[!duplicated(cond),]
#save_points_subset <- cumsum(sapply(list(cond1, cond2, cond3, cond4, cond5, cond6), nrow))

#****************************************
# Setup parallel simulation
#****************************************
reps <- 1e4
out_list <- vector("list", length = nrow(cond))

last_save <- 0 #last condition row that was saved, see end of simulation below
intermediate_save <- FALSE
save_points <- 1000 #save simulation results after how many conditions have run
#save_points <- save_points_subset

ncores <-parallel::detectCores()
cl <- parabar::start_backend(ncores) # Create cluster based on nworkers.

parabar::evaluate(cl, library(metafor))
parabar::export(cl, c("simulate_rma", "rnorm_truncated", "try_run"))

#****************************************
# Run simulation
#****************************************

start <- Sys.time()
for(r in last_save+1:nrow(cond)){ #gives us a list of lists

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
    non_converged <- reps - nrow(out_list[[r]])
    out_list[[r]] <- out_list[[r]][, lapply(.SD, mean)] #data.table colMeans but returns a dataframe (well, data.table)
    out_list[[r]] <- cbind(out_list[[r]], cond[r,], non_converged)

    if(intermediate_save){
        if(r %in% save_points | r == nrow(cond)){
            #if reached a save point or simulation finished
            #turn results into dataframe and save as csv

            save_range <- (last_save + 1):r

            e <- rbindlist(out_list[save_range])

            file_name <- paste0("./data/means_cond_",
                                save_range[1],"-", save_range[r],
                                ".csv")

            fwrite(x = e, file = file_name)
            rm(e)
            last_save <- r
        }
    }

}

stop_backend(cl) # Shut down the nodes and end simulation

end <- Sys.time()
end - start


#****************************************
# Clean up and save results
#****************************************

#e <- rbindlist(out_list)

#data.table::fwrite(x = e, file = "/data/means_combined.csv")
