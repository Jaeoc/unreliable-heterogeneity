# Project: Unexplained heterogeneity
# Code purpose: cleaning raw data
# Author: Anton Olsson-Collentine (anton@olssoncollentine.com)


#****************************************
# functions and packages
#****************************************

source("./code/functions.r") # for create_condition_cols function

#****************************************
# Clean data raw data
#****************************************

#With r
dat_r <- readRDS("data/raw_r_tau_0-0.2.RDS")

#Turn list into dataframe using list-names as idcol
dat_r <- data.table::rbindlist(dat_r, idcol = "mu")

#create condition ids
dat_r <- create_condition_cols(dat_r)

saveRDS(dat_r, "data/cleaned_r_tau_0-0.2.RDS")


# With z
dat_z <- readRDS("data/raw_z_tau_0-0.2.RDS")

#Turn list into dataframe using list-names as idcol
dat_z <- data.table::rbindlist(dat_z, idcol = "mu")

#create condition ids
dat_z <- create_condition_cols(dat_z)

saveRDS(dat_z, "data/cleaned_z_tau_0-0.2.RDS")

e_means <- lapply(dat_z, colMeans)
saveRDS(e_means, "data/means_z_tau_0-0.2.RDS")

#e_medians <- lapply(e, function(x) apply(x, 2, median))

# For HS method

dat_hs <- readRDS("data/raw_r_HS_tau_0-0.2.RDS")

#Turn list into dataframe using list-names as idcol
dat_hs <- data.table::rbindlist(dat_hs, idcol = "mu")

#create condition ids
dat_hs <- create_condition_cols(dat_hs)

saveRDS(dat_hs, "data/cleaned_hs_tau_0-0.2.RDS")

e_means <- lapply(dat_hs, colMeans)
saveRDS(e_means, "data/means_hs_tau_0-0.2.RDS")

# For low tau

dat_low <- readRDS("data/raw_r_tau_0.02-0.08.RDS")

#Turn list into dataframe using list-names as idcol
dat_low <- data.table::rbindlist(dat_low, idcol = "mu")

#create condition ids
dat_low <- create_condition_cols(dat_low)

saveRDS(dat_low, "data/cleaned_r_tau_0.02-0.08.RDS")

e_means <- lapply(dat_low, colMeans)
saveRDS(e_means, "data/means_r_tau_0.02-0.08.RDS")


#****************************************
# Create means
#****************************************

#With r
dat_r <- readRDS("data/raw_r_tau_0-0.2.RDS")
r_means <- lapply(dat_r, colMeans)

#Turn list into dataframe using list-names as idcol
r_means <- lapply(r_means, function(x) as.data.frame(as.list(x)))
r_means <- data.table::rbindlist(r_means, idcol = "mu")

#create condition ids
r_means <- create_condition_cols(r_means)

saveRDS(r_means, "data/means_r_tau_0-0.2.RDS")

#With z
dat_z <- readRDS("data/raw_z_tau_0-0.2.RDS")
z_means <- lapply(dat_z, colMeans)

#Turn list into dataframe using list-names as idcol
z_means <- lapply(z_means, function(x) as.data.frame(as.list(x)))
z_means <- data.table::rbindlist(z_means, idcol = "mu")

#create condition ids
z_means <- create_condition_cols(z_means)

saveRDS(z_means, "data/means_z_tau_0-0.2.RDS")


#With hs
dat_hs <- readRDS("data/raw_r_hs_tau_0-0.2.RDS")
hs_means <- lapply(dat_hs, colMeans)

#Turn list into dataframe using list-names as idcol
hs_means <- lapply(hs_means, function(x) as.data.frame(as.list(x)))
hs_means <- data.table::rbindlist(hs_means, idcol = "mu")

#create condition ids
hs_means <- create_condition_cols(hs_means)

saveRDS(hs_means, "data/means_hs_tau_0-0.2.RDS")

#With low tau
dat_low <- readRDS("data/raw_r_tau_0.02-0.08.RDS")
low_means <- lapply(dat_low, colMeans)

#Turn list into dataframe using list-names as idcol
low_means <- lapply(low_means, function(x) as.data.frame(as.list(x)))
low_means <- data.table::rbindlist(low_means, idcol = "mu")

#create condition ids
low_means <- create_condition_cols(low_means)

saveRDS(low_means, "data/means_r_tau_0.02-0.08.RDS")
