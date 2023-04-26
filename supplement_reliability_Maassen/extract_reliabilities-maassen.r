# Project: Unexplained heterogeneity
# Code purpose: estimate typical reliability and its SD from E. Maassen's data
# Author: Anton Olsson-Collentine (anton@olssoncollentine.com)


#******************************************
library(openxlsx)


dat <- openxlsx::read.xlsx("supplement_reliability_Maassen/codebook-main-step4_rel.xlsx")

#From Esther:
#Your columns of interest are rel_rep and rel_est.
#(i.e., reported and estimated reliabilities)

#Note that this reliability may mean nothing, given that it’s on all the items, and therefore one factor is assumed, even though there may be more factors (f_total) underlying the data.

mean(dat$rel_rep)
#[1] 0.8392653
sd(dat$rel_rep)
#[1] 0.0997923

mean(dat$rel_est, na.rm = TRUE)
#[1] 0.8657143
sd(dat$rel_est, na.rm = TRUE)
#[1] 0.07702407

# Higher estimates and with less SD when recomputed. Likely a because of a correlation with who made their data openly available.


dat$comp_rel <- dat$rel_rep - dat$rel_est
mean(dat$comp_rel, na.rm = TRUE)
