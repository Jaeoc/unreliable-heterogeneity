# Project: Unexplained heterogeneity
# Script purpose: estimate common tau2 values for correlations in psychology
# based on van Erp et al. (2017): 10.5334/jopd.33
# Code: Anton Olsson-Collentine

#**********************************************
# Setup
#**********************************************
# Data downloaded from https://osf.io/p6ecw

dat <- read.csv("supplemental_data/Data 1990-2013 with tau values.csv")

unique(dat$Type.of.ES)
#also contains corrected Pearson's r

dat_r <- dat[dat$Type.of.ES == "Pearson's r",]

hist(dat_r$tau.2, breaks = 30)
boxplot(dat_r$tau.2)
# We see that 50% of tau2 values are between about ~0.02- 0.06,
#with another 25% between 0.06- 0.125 and 25% below 0.02
quantile(dat_r$tau)
#Actually, the 25th quantile is 0.01, median 0.03, and the 75th is 0.0576
#the maximum value is 0.2704, which gives tau = 0.52


# For corrected pearson's rs
dat_r_corrected <- dat[dat$Type.of.ES == "Corrected Pearson's r",]
hist(dat_r_corrected$tau.2, breaks = 30)
boxplot(dat_r_corrected$tau.2)
quantile(dat_r_corrected$tau)
# The corrected rs have a lower distribution
# 25th quantiles = 0.001, median = 0.01, 75th = 0.0225,
#The maximum is 0.2209

quantile(dat_r_corrected$X..of.effect.sizes)
dat_r_corrected2 <- dat_r_corrected[dat_r_corrected$X..of.effect.sizes < 1000,]
quantile(dat_r_corrected2$X..of.effect.sizes)


# Check K
dat_r_any <- dat[grepl("Pearson", dat$Type.of.ES),]
unique(dat_r_any$Type.of.ES)
quantile(dat_r_any$X..of.effect.sizes)
hist(dat_r_any$X..of.effect.sizes)

quantile(dat_r$X..of.effect.sizes)
hist(dat_r$X..of.effect.sizes)

dat_r_any2 <- dat_r_any[dat_r_any$X..of.effect.sizes < 1000,]
quantile(dat_r_any2$X..of.effect.sizes)

# Check correction or not
dat_r = dat[dat$Type.of.ES %in% c("Pearson's r", "Weighted Pearson's r", "Unweighted Pearson's r"),]
unique(dat_r$Reference) #13

dat_r_corrected = dat[dat$Type.of.ES %in% c("Corrected Pearson's r", "Weighted corrected Pearson's r"),]
unique(dat_r_corrected$Reference) #13

# Corrected generally
dat_cor <- dat[dat$Type.of.ES %in% c("Corrected Pearson's r",
       "Weighted corrected Pearson's r",
       "Corrected Cohen's d"),]
unique(dat_cor$Reference) #14

length(unique(dat$Reference)) - length(unique(dat))
#49 uncorrected and 14 corrected
