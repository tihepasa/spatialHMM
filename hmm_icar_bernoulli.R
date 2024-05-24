# estimate the spatial HMM with S states
library(rstan)
library(dplyr)
library(tidyr)

# data includes:
# T, number of time points
# N, number of regions
# n_neighbours, number of neighbours for each region as a vector
# neighbours, matrix of neighbours: row represents a region, each cell holds the id number of the neighbour
# is_measles, matrix of deaths: row represents a region, column a date
# missing, matrix indicating missingness: row represents a region, column a date; 0 = present, 1 = missing
par <- readRDS("measlesdata.rds")

# matrix of all neighbour pairs with two columns holding the site ids
edges <- vector("list", par$N)
for(i in 1:par$N) {
  n <- par$neighbours[i, 1:par$n_neighbours[i]]
  n <- n[n > i]
  if(length(n) > 0) {
    edges[[i]] <- cbind(i, n)
  } else {
    edges[[i]] <- NULL
  }
}
edges <- do.call("rbind", edges)


# indices for observed cells (row, column)
ind_temp <- which(as.matrix(par$missing) == 0, arr.ind = TRUE)
ind <- ind_temp[, 1]
# numbers of observed cells in every column (sites in each time point)
ns <- c("0" = 0, table(ind_temp[, 2]))
# cumulative number of observations
cum_n <- cumsum(ns)

# number of observed cells
n_obs <- sum(as.matrix(par$missing != 1))

# matrix form of the response
par$is_measles <- as.matrix(par$is_measles)

# replace NAs with zeros, this is only for stan to approve the data, these zeros are not used
par$is_measles[par$missing == 1] <- 0

# number of states, could be 4 or 6 or something else
S <- 5

# input data for the model
pars <- list('T' = par$T,
             N = par$N,
             S = S,
             deaths = t(par$is_measles),
             N_edges = nrow(edges),
             node1 = edges[, 1],
             node2 = edges[, 2],
             n_obs = n_obs,
             ind = ind,
             cum_n = cum_n)

# use different initial values for all 4 chains (at least for most of the parameters)
# just for an additional proof of some robustness in terms of multimodality
inits <- replicate(4, list(A = diag(1 - 0.05 * pars$S, pars$S, pars$S) + 0.05,
                           rho = rep(1 / pars$S, pars$S),
                           mu_1 = runif(1, -4.5, -4),
                           mu_S = runif(1, -2, -1.5),
                           m = c(MCMCpack::rdirichlet(1, rep(20, pars$S - 1))),
                           lambda_raw = rnorm(pars$N - 1, sd = 0.1),
                           sigma_lambda = runif(1, 0.1, 1),
                           sigma_phi = runif(1, 0.2, 1.2),
                           phi = matrix(0, nrow = pars$N, ncol = pars$S),
                           gamma_raw = rnorm(11, sd = 0.1)),
                   simplify = FALSE)

library(cmdstanr)
model <- cmdstan_model("hmm_icar_bernoulli.stan", stanc_options=list("O1"))
set.seed(1)

# fit the model,  10000 posterior samples after the warmup
fit <- model$sample(data = pars, init = inits,
                    iter_sampling = 2500, iter_warmup = 5000,
                    chains = 4, parallel_chains = 4, refresh = 100, save_warmup = FALSE)

# save the model with a suitable name for future use as a CmdStan object or as a Stan object
# the stan fit is needed for the changepoint analysis
#fit$save_object(paste0("cmdstanfit_hmm_icar_bernoulli.rds"))
out <- rstan::read_stan_csv(fit$output_files())
saveRDS(out, file = paste0("fit_hmm_icar_bernoulli.rds"))


## some results

# some parameter estimates
print(fit, pars = c("A", "rho", "mu", "sigma_phi", "sigma_lambda"))

# effective sample sizes
e <- rstan::extract(fit, pars = c("x", "log_omega"), include = FALSE, permuted = FALSE)

eb <- apply(e, 3, ess_bulk)
summary(eb)

et <- apply(e, 3, ess_tail)
summary(et)


# the most likely states
x <- rstan::extract(fit, "x")$x
map <- as.integer(apply(x, 2, function(x) names(table(x))[which.max(table(x))])) # modes

# the probabilities in the most likely states (computation takes a couple of minutes)
month <- rep(1:12, 101)
prob <- matrix(NA, nrow = par$N, ncol = par$T)
for (i in 1:par$T) {
  prob[, i] <- colMeans(plogis(e$mu[, map[i]] + e$lambda + c(e$sigma_phi) * e$phi[, , map[i]] + e$gamma[month[i]]))
}

ts.plot(colMeans(prob))