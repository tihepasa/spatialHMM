# estimate the spatial HMM with S states
library(cmdstanr)
library(rstan)
library(dplyr)
library(tidyr)
library(loo)

#model <- cmdstan_model("hmm_icar_bernoulli_loo.stan", stanc_options=list("O1"))
set.seed(1)

# number of states, could be 4 or 6 or something else
S1 <- 4
S2 <- 5

# read the saved model object
fit1 <- readRDS(paste0("hmm_icar_bern_", S1, "s.rds"))
fit2 <- readRDS(paste0("hmm_icar_bern_", S2, "s.rds"))

# extract draws
draws1 <- fit1$draws(variables = "log_lik")
draws2 <- fit2$draws(variables = "log_lik")

# from the above draws it is convenient to estimate the ELPDs
# relative effective sample sizes have to be added manually
r_eff1 <- relative_eff(exp(draws1))
r_eff2 <- relative_eff(exp(draws2))

# compute the ELPDs
l1 <- loo::loo(draws1, r_eff = r_eff1, cores = 4) # adjust the number of cores for the system at hand
l2 <- loo::loo(draws2, r_eff = r_eff2, cores = 4)

# results of the individual models
l1
l2

# compare the models
loo::loo_compare(l1, l2)
