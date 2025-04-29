// spatial hidden Markov model for one disease and S states
functions {
  // Function for computing likelihood
  real loglik(matrix log_py, int M, int T, matrix log_A, vector log_rho) {
    vector[M] log_alpha;
    vector[M] log_alpha_new;
    vector[M] tmp;
    // start from the first state
    log_alpha = log_rho + log_py[ : , 1];
    for (t in 2 : T) {
      for (k in 1 : M) {
        int n_finite = 0;
        tmp = log_alpha + log_A[ : , k] + log_py[k, t];
        // reorder elements of tmp so that all finite values are at the start of the vector
        for (i in 1 : M) {
          if (tmp[i] > negative_infinity()) {
            n_finite += 1;
            tmp[n_finite] = tmp[i];
          }
        }
        if (n_finite == 0) {
          log_alpha_new[k] = negative_infinity();
        } else {
          log_alpha_new[k] = log_sum_exp(tmp[1 : n_finite]);
        }
      }
      log_alpha = log_alpha_new;
    }
    return log_sum_exp(log_alpha);
  }
  // Function for computing posterior probabilities of states in log-scale
  matrix log_posterior_prob(matrix log_py, int S, int T, matrix log_A, vector log_rho) {
    matrix[S, T] log_alpha;
    matrix[S, T] log_beta;
    vector[S] tmp;
    // forward
    log_alpha[, 1] = log_rho + log_py[ : , 1];
    for (t in 2:T) {
      for (s in 1:S) {
        log_alpha[s, t] = log_sum_exp(log_alpha[, t - 1] + log_A[:, s] + log_py[s, t]);
      }
    }
    // backward
    log_beta[, T] = rep_vector(0, S);
    for (tt in 1:(T - 1)) {
      int t = T - tt;
      for (s in 1:S) {
        log_beta[s, t] = log_sum_exp(log_A[s, :]' + log_beta[:, t + 1] + log_py[, t + 1]);
      }
    }
    return log_alpha + log_beta - log_sum_exp(log_alpha[, T]);
  }
  // Function for computing Q of QR decomposition
  vector create_Q(int N) {
    vector[2 * N] Q;
    for (i in 1 : N) {
      Q[i] = -sqrt((N - i) / (N - i + 1.0));
      Q[i + N] = inv_sqrt((N - i) * (N - i + 1));
    }
    return Q;
  }
  // sum-to-zero constrained x from x_raw
  vector sum_to_zero(vector x_raw, vector Q, int N) {
    vector[N] x;
    real x_aux = 0;
    for (i in 1 : (N - 1)) {
      x[i] = x_aux + x_raw[i] * Q[i];
      x_aux = x_aux + x_raw[i] * Q[i + N];
    }
    x[N] = x_aux;
    return x;
  }
}
data {
  int<lower = 0> T; // number of time points
  int<lower = 0> N; // number of regions
  int<lower = 1> S; // number of hidden states
  array[T, N] int<lower = 0, upper = 1> deaths; // response variable, here dichotomous deaths
  // helper variables to index the neighbors
  int<lower = 0> N_edges;
  array[N_edges] int<lower = 1, upper = N> node1;
  array[N_edges] int<lower = 1, upper = N> node2;
  int<lower = 0> n_obs; // number of observations
  array[n_obs] int<lower = 0> ind; // (site)indeces of non-missing observations
  array[T + 1] int<lower = 0> cum_n; // cumulative sum of nationwide monthly non-missing observations, i.e. cumsum(n_t)
}
transformed data {
  // matrix of the deaths array, initialize to contain 2s for generated quantities block
  // otherwise fill the observed cells with observed values
  matrix[N, T] deaths_matrix = rep_matrix(2, N, T);
  for(t in 1:T) {
    deaths_matrix[ind[(cum_n[t] + 1):cum_n[t + 1]], t] = to_vector(deaths[t, ind[(cum_n[t] + 1):cum_n[t + 1]]]);
  }
  
  // scale parameter for the sum constrained prior
  real scaleN = inv(sqrt(1 - inv(N)));
  real scale12 = inv(sqrt(1 - inv(12)));
  
  // create QR decomposition for sum-to-zero constraint
  vector[2 * N] QN = create_Q(N);
  vector[2 * 12] Q12 = create_Q(12);
  
  // indices for seasonal effect
  array[T] int month;
  for (t in 1:T) {
    month[t] = 1 + (t - 1) % 12;
  }
  
  // helper matrix for giving parameteres for the dirichlet prior of the transition matrix A
  matrix[S, S] prior_A = rep_matrix(0.5, S, S);
  for(s in 1:S) prior_A[s, s] = 2 * S;
  // favors a bit cases where we stay in one state at least for a while
  // apply(MCMCpack::rdirichlet(10000, c(2 * S, rep(0.5, S - 1))), 2, quantile, probs = c(0.005,0.5,0.995)) 
  
  vector[S - 1] prior_m = rep_vector(5, S - 1);
  // apply(MCMCpack::rdirichlet(1e5,5*rep(1, 4)), 2, quantile, probs = c(0.005, 0.5, 0.995))
}
parameters {
  // SxS transition matrix, no prior => uniform prior for each row
  array[S] simplex[S] A;
  
  // initial probabilities
  simplex[S] rho;
  
  // just one sigma_phi, try to make the model a bit less flexible (multimodal)
  real<lower = 0> sigma_phi;
  
  // ICAR components, one per state
  matrix[N, S] phi;
  
  // parts of the ordered common intercept, one per state (phi's are constrained to mean zero)
  real mu_1; // intercept for the first state
  real<lower = mu_1> mu_S; // intercept for the last state
  simplex[S - 1] m; // for the intercepts between
  
  // site related intercepts and their deviation
  vector[N - 1] lambda_raw;
  real<lower=0> sigma_lambda;
  
  // raw parameters for seasonal effects with sum-to-zero constraint,
  // common for all states, same reasoning as with sigma_phi
  vector[11] gamma_raw;
}
transformed parameters {
  // log(p(y_1,t, ..., y_N,t | state = s))
  matrix[S, T] log_omega;
  
  // hmm_marginal needs transition _matrix_
  matrix[S, S] Gamma;
  for(s in 1:S) {
    Gamma[s,] = A[s]';
  }
  
  // set sum-to-zero constraints for seasonal effects
  vector[12] gamma = sum_to_zero(gamma_raw, Q12, 12);
  
  // fix mean of lambda to 0
  vector[N] lambda = sigma_lambda * sum_to_zero(lambda_raw, QN, N);
  
  // actual state related constants
  vector[S] mu = mu_1 + (mu_S - mu_1) * cumulative_sum(append_row(0, m));
  
  {
    array[12] matrix[N, S] logit_p;
    array[12] matrix[N, S] logit_p2;
    // logit(p) for each region when in state s and month i
    for (i in 1:12) {
      for (s in 1:S) {
        logit_p[i, , s] = gamma[i] + mu[s] + lambda + sigma_phi * phi[, s]; // linear predictor
        logit_p2[i, , s] = log1p_exp(logit_p[i, , s]); // other term needed for bernoulli log-pmf
      }
    }
    for (t in 1:T) {
      for (s in 1:S) {
        // log-bernoulli = y*logit_p - log(1 + exp(logit_p))
        // death_matrix is the original deaths array transformed to matrix in the transformed data block
        log_omega[s, t] = sum(deaths_matrix[ind[(cum_n[t] + 1):cum_n[t + 1]], t] .* 
        logit_p[month[t], ind[(cum_n[t] + 1):cum_n[t + 1]], s] - 
        logit_p2[month[t], ind[(cum_n[t] + 1):cum_n[t + 1]], s]);
      }
    }
  }
}
model {
  // seasonality
  gamma_raw ~ normal(0, scale12);
  
  // prior for the average level for the first state (which due to the ordering is the one with the smallest death probability)
  // 1% change (logit(0.01) = -4.59
  // These are based on 1750-1850 for measles
  // from yearly means:
  // qlogis(quantile(deaths, probs = seq(0.1, 0.9, length = 5)))
  // 10%  30%  50%  70%  90% 
  // 0.01 0.02 0.03 0.07 0.15 -> logit -4.85 -3.99 -3.36 -2.63 -1.74 
  mu_1 ~ normal(-4.5, 0.25);
  // this corresponds to prior quantiles 
  //  0.5%   50% 99.5% 
  // 0.006 0.011 0.021 
  mu_S ~ normal(-1.75, 0.5)T[mu_1, ]; 
  // this corresponds to prior quantiles 
  //  0.5%   50% 99.5% 
  // 0.045 0.148 0.384 
  // favor cases where all simplex components >> 0
  m ~ dirichlet(prior_m);
  
  lambda_raw ~ normal(0, scaleN);
  sigma_lambda ~ std_normal();
  
  sigma_phi ~ std_normal();
  
  // phi and soft-constraint
  for (s in 1:S) {
    A[s] ~ dirichlet(prior_A[, s]); 
    target += -0.5 * dot_self(phi[node1, s] - phi[node2, s]);
    sum(phi[, s]) ~ normal(0, 0.001 * N);
  }
  
//  target += loglik(log_omega, S, T, log(Gamma), log(rho)); // self implemented version to avoid numerical issues
  target += hmm_marginal(log_omega, Gamma, rho);
}
generated quantities {
  // samples of latent state trajectories
  array[T] int x = hmm_latent_rng(log_omega, Gamma, rho);
  
  // get the log likelihood values, formula derived in supplements
  // comment out if not needed since this takes time and memory
  vector[n_obs] log_lik;
  {
    // no need to save b, or X; hence the local environment
    vector[N] ll;
    vector[S] log_b; // 1 / P(y_it = m | z_t = s), s = 1,..., S
    // marginal log-posterior probabilities of each hidden state value
    matrix[S, T] log_X = log_posterior_prob(log_omega, S, T, log(Gamma), log(rho));
    //matrix[S, T] log_X = log(hmm_hidden_state_prob(log_omega, Gamma, rho)); // numerical issues, does not work
    for (t in 1:T) {
      for (i in 1:N) {
        if (deaths_matrix[i, t] < 2) {
          for (s in 1:S) {
            log_b[s] = bernoulli_logit_lpmf(deaths[t, i] | gamma[month[t]] + mu[s] + lambda[i] + sigma_phi * phi[i, s]);
          }
          ll[i] = -log_sum_exp(log_X[, t] - log_b);
        }
      }
      log_lik[(cum_n[t] + 1):cum_n[t + 1]] = ll[ind[(cum_n[t] + 1):cum_n[t + 1]]];
    }
  }
}
