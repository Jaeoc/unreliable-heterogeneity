# Project: Unexplained heterogeneity
# Code purpose: Checking other's corrections for unreliability
# Author: Anton Olsson-Collentine (anton@olssoncollentine.com)




# Metafor example
# https://www.metafor-project.org/doku.php/tips:hunter_schmidt_method

library(metafor)
dat <- dat.mcdaniel1994
dat <- escalc(measure="COR", ri=ri, ni=ni, data=dat, vtype="AV")
res <- rma(yi, vi, weights=ni, data=dat, method="HS")
res

# NB! Wolfgang shows increase in heterogeneity after corrected
# for unreliability!


# Psych-meta example

# Data from: https://cran.r-project.org/web/packages/psychmeta/vignettes/ma_r.html
# But only for X and Y (not Z)
dat <- matrix(
    c("1", "1", "X", "Y", "416", "0.49", "0.79", "0.77", "Watson2005",
    "2", "1", "X", "Y", "241", "0.54", "0.82", "0.84", "Watson2005",
    "3", "1", "X", "Y", "479", "0.34", "0.73", "0.87", "Zellars2006",
    "4", "1", "X", "Y", "167", "0.37", "0.78", "0.79", "Bandura1980"),
 ncol = 9, byrow = TRUE)

dat <- as.data.frame(dat)
names(dat) <- c("sample_id", "moderator","x_name","y_name","n","rxyi", "rxxi","ryyi", "citekey")
dat[,5] <- as.numeric(dat[,5])
dat[,6] <- as.numeric(dat[,6])
dat[,7] <- as.numeric(dat[,7])
dat[,8] <- as.numeric(dat[,8])


install.packages("psychmeta")
library(psychmeta)

# Without correction for unreliability
fit1 <- ma_r(rxyi = rxyi,
               n = n,
               construct_x = x_name,
               construct_y = y_name,
               rxx = rxxi,
               ryy = ryyi,
               sample_id = sample_id,
               data = dat
               )

summary(fit1)
# tau  0.0842 (sd_res)

# https://cran.r-project.org/web/packages/psychmeta/vignettes/overview.html+
# section "Individual-correction meta-analyses"

fit2 <- ma_r(ma_method = "ic",
               rxyi = rxyi,
               n = n,
               construct_x = x_name,
               construct_y = y_name,
               rxx = rxxi,
               ryy = ryyi,
               clean_artifacts = FALSE,
               impute_artifacts = FALSE,
               data = dat)
summary(fit2)

#tau = 0.102 (sd_rho)

# Also here we see an increase in heterogeneity after correction

# Borenstein example (p. 346)

dat <- matrix(c("University_1","130","0.24","0.75",
"University_2","90","0.11","0.75",
"Private_1","30","0.05","0.60",
"Private_2","25","0.17","0.60",
"Volunteer_1","50","0.38","0.90",
"Volunteer_2","65","0.50","0.90"), nrow = 6, byrow = TRUE)

dat <- as.data.frame(dat)
names(dat) <- c("study", "n", "r", "ryy")
dat$rxx <-  1

dat$r <- as.numeric(dat$r)
dat$ryy <- as.numeric(dat$ryy)
dat$n <- as.numeric(dat$n)

# Metafor approach
dat <- escalc(measure = "COR", ri = r, ni = n, data = dat, vtype = "AV")
res <- rma(yi, vi, weights = n, method = "HS", data = dat)

# tau = 0.0828 (tau2 = 0.0069)

dat$rho <- dat$r / sqrt(dat$ryy)
dat$rho_var <- dat$vi / dat$ryy
res2 <- rma(rho, rho_var, weights=1/rho_var, data=dat, method="HS")

# tau = 0.0649 (tau2 = 0.0042)

# These values correspond to what Borenstein compute on p. 347-348

# Hunter & Schmidt
# Tend to assume tau = 0 in their examples, for example p. 130
#(2nd edition, also 3rd edition, but note sure what page)

# However, the artifact distribution chapter has a non-zero example (3rd ed.)
# See obsidian notes
