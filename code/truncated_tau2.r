

es <- seq(0, 0.6, 0.1)
vars <- c(0.002, 0.0035, 0.0055, 0.0085, 0.012, 0.0185, 0.031, 0.069)

a <- lapply(vars, function(x){

lapply(es, function(y) {
res <- rnorm_truncated(n = 1e6, mean = y, sd = sqrt(x),
lower_bound = -1, upper_bound = 1)
data.frame(tau2 = var(res), mu = y, nominal_tau2 = x)
})

})

b <- lapply(a, function(x) do.call(rbind, x))
b <- do.call(rbind, b)
saveRDS(b, "data/truncated_tau2.RDS")
