# estimate the change point with a two-sate left-to-right hmm
library(rstan)
library(cmdstanr)

# read the fit holding the state trajectories
fit_measles <- readRDS("hmm_icar_bernoulli.rds")

# extract and choose 1000 trajectory samples
states <- extract(fit_measles, "x")[[1]][1:1000,]


# LR-HMM model
model <- cmdstan_model("hmm_changepoint.stan", stanc_options = list("O1"))

# initial values
init_B <- rbind(prop.table(table(head(t(states), ncol(states) / 2))),
                prop.table(table(tail(t(states), ncol(states) / 2))))

# number of states in the trajectory
K = 5

# run the model
fit_hmm <- model$sample(data = list(obs = states, N = nrow(states), T = ncol(states), K = K), 
                        init = replicate(4, list(B = init_B, move_probability = 1/nrow(states)), simplify = FALSE),
                        parallel_chains = 4, chains = 4, refresh = 10)

changepoints <- changepoint_fit$draws("changepoints")


# change point as a date
barplot(table(c(changepoints)))
prop.table(table(c(changepoints)))

library(lubridate)
dates <- ordered(seq(ymd("1750-01-01"), ymd("1850-12-01"), by = "months"))
quantile(dates[c(changepoints)], probs = c(0.025, 0.5, 0.975),type = 1)
#      2.5%        50%      97.5% 
# 1812-09-01 1812-11-01 1813-02-01 


# emission matrix
summary(B)

apply(B, 3, mean)
apply(B, 3, quantile, .025)
apply(B, 3, quantile, .975)

