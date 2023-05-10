


dat <- matrix(c("University_1","130","0.24","0.75",
"University_2","90","0.11","0.75",
"Private_1","30","0.05","0.60",
"Private_2","25","0.17","0.60",
"Volunteer_1","50","0.38","0.90",
"Volunteer_2","65","0.50","0.90"), nrow = 6, byrow = TRUE)

dat <- as.data.frame(dat)
names(dat) <- c("study", "n", "r", "rel")

dat$r <- as.numeric(dat$r)
dat$rel <- as.numeric(dat$rel)
dat$n <- as.numeric(dat$n)
dat$r_var <- (1 - dat$r^2)^2 / (dat$n - 1)

dat$rho <- dat$r / sqrt(dat$rel)
dat$rho_var <- (1 - dat$rho^2)^2 / (dat$n - 1)

dat$rho_var2 <- dat$r_var / sqrt(dat$rel)
#estimate rho var based on observed r

library(metafor)
rma(yi = dat$r, vi = dat$r_var)
rma(yi = dat$rho, vi = dat$rho_var)

rma(yi = dat$rho, vi = dat$rho_var2)
# If I use the variance based on the estimated r
# I get smaller heterogeneity with true effect sizes!

# 2 ideas that might be the difference
# a) add measurement error before sampling error
#   - This does not explain it, I just tried with my
#     function and it still underestimates tau by about 0.03
# b) Compute observed correlations from a sample
# (NB. b) interacts with a))
#   - This also doesn't explain it, I just used Brannick's
# code below and there is still an underestimate when using
# meta-analyzing unreliabile data

# Brannick et al., simulation
TxTy <- rnorm(10, mean = .26, sd = .13)
rel_x <- runif(10, min = 0, max = 1)
rel_y <- runif(10, min = 0, max = 1)

TxX <- sqrt(rel_X)
TyY <- sqrt(rel_y)
xy <- TxTy*TxX*TyY

#Compute attenuated correlations
n <- 142

rhobox <- matrix(0,4,4)+diag(4) # placeholder for correlation matrix for sampling
mu <- rep(0,4)          # multivariate means of zero
names.sim <- cbind("Tx", "Ty", "X", "Y")
colnames(rhobox)<- names.sim
rownames(rhobox)<- names.sim

obs <- vector("double", 10)
res <- vector("double", 1e3)

for(rep in 1:1e3){

for(j in 1:10){

ryyi <- 0.8
rxxi <- 0.8
rho.r <- rnorm(1, mean = .26, sd = .13)
att.rho.r <- rho.r*sqrt(rxxi*ryyi)              # attenuate rho for reliability

# find the population correlation matrix from which to sample observations
rhobox[2,1] <- rho.r                            # correlation between true scores
rhobox[1,2] <- rho.r
rhobox[3,1] <- sqrt(rxxi)                       # true x and observed X
rhobox[1,3] <- sqrt(rxxi)
rhobox[3,2] <- rho.r*sqrt(rxxi)                 # true Y and observed X
rhobox[2,3] <- rho.r*sqrt(rxxi)
rhobox[4,2] <- sqrt(ryyi)                       # true Y and observed Y
rhobox[2,4] <- sqrt(ryyi)
rhobox[4,1] <- rho.r*sqrt(ryyi)                 # true X and observed Y
rhobox[1,4] <- rho.r*sqrt(ryyi)
rhobox[4,3] <- att.rho.r                        # observed X and Y at pop level
rhobox[3,4] <- att.rho.r


sample.rs <-mvrnorm(n,mu=mu,Sigma=rhobox) # sample ni data for obs r from attenuated rho
cor1 <-cor(sample.rs)[3, 4]                           # compute correlation from
obs[j] <- cor1
}

obs_var <- (1 - obs^2)^2 / (n - 1)

res[rep] <- metafor::rma(yi = obs, vi = obs_var)$tau2
mes <- paste0(rep, "\n")
cat(mes)
}

mean(sqrt(res))
