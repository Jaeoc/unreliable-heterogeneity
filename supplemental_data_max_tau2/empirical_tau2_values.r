# Project: Unexplained heterogeneity
# Script purpose: estimate common tau2 values for correlations in psychology
# based on van Erp et al. (2017): 10.5334/jopd.33
# Code: Anton Olsson-Collentine

#**********************************************
# Setup
#**********************************************
# Data downloaded from https://osf.io/p6ecw

dat <- read.csv("supplemental_data_I2/Data 1990-2013 with tau values.csv")

unique(dat$Type.of.ES)
#also contains corrected Pearson's r

dat_r <- dat[dat$Type.of.ES == "Pearson's r",]

hist(dat_r$tau.2, breaks = 30)
boxplot(dat_r$tau.2)
# We see that 50% of tau2 values are between about ~0.02- 0.06,
#with another 25% between 0.06- 0.125 and 25% below 0.02
quantile(dat_r$tau.2)
#Actually, the 25th quantile is 0.01, median 0.03, and the 75th is 0.0576
#the maximum value is 0.2704, which gives tau = 0.52


# For corrected pearson's rs
dat_r_corrected <- dat[dat$Type.of.ES == "Corrected Pearson's r",]
hist(dat_r_corrected$tau.2, breaks = 30)
boxplot(dat_r_corrected$tau.2)
quantile(dat_r_corrected$tau.2)
# The corrected rs have a lower distribution
# 25th quantiles = 0.001, median = 0.01, 75th = 0.0225,
#The maximum is 0.2209
