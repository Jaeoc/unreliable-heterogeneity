
source("./code/functions.r")
es <- seq(0, 0.6, 0.1)
tau <- c(0, 0.1, 0.15, 0.2)

a <- lapply(tau, function(x){

lapply(es, function(y) {
res <- rnorm_truncated(n = 1e6, mean = y, sd = x,
lower_bound = -1, upper_bound = 1)

data.frame(tau = sd(res), mu = y, nominal_tau = x)
})

})

b <- lapply(a, function(x) do.call(rbind, x))
b <- do.call(rbind, b)
saveRDS(b, "data/truncated_tau.RDS")
