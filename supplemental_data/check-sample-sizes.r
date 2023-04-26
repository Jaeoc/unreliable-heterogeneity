# Check sample sizes in https://journals.sagepub.com/doi/full/10.1177/2515245919838781



library(haven)

dat <- haven::read_sav("supplemental_data_max_tau2/data.sav")
quantile(dat$Sample)
#1300 rows

dat2 <- dat[dat$Sample < 1e4,]
#1291 rows
quantile(dat2$Sample, na.rm = TRUE)

dat3 <- dat2[dat2$Sample < 5e3,]
#1290 rows
quantile(dat3$Sample)
