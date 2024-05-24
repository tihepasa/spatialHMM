// two-state left-to-right HMM for finding the change point
functions {
  // likelihood
  real loglik(matrix log_py, matrix log_A) {
    int T = cols(log_py);
    vector[2] log_alpha;
    vector[2] log_alpha_new;
    // start from the first state
    log_alpha[1] = log_py[1, 1];
    log_alpha[2] = negative_infinity();
    for (t in 2:T) {
      // first state, P(z_{t-1} = 1)P(z_t=1 | z_{t-1} = 1)P(y_t | z_t = 1)
      log_alpha_new[1] = log_alpha[1] + log_A[1, 1] + log_py[1, t];
      // move to 2 or stay at 1
      log_alpha_new[2] = log_sum_exp(log_alpha + log_A[, 2] + log_py[2, t]);
      log_alpha = log_alpha_new;
    }
    return log_sum_exp(log_alpha);
  }
}
data {
  int<lower=0> T; // number of time points
  int<lower=1> N; // number of samples
  int<lower=1> K; // number of observation categories
  array[N, T] int<lower = 1, upper = K> obs; // observations
}
// start from first state
transformed data {
  vector[2] rho = [1, 0]';
}
parameters {
  // probability of transitioning
  real<lower=0,upper=1> move_probability;
  // emission probabilities
  array[2] simplex[K] B;
}
transformed parameters {
  // log_p(y), could be faster to modify the loglik function so that it uses log_B directly instead of computing log_omega
  // although we need log_omega for hmm_latent_rng so that should be computed in generated quantities anyway
  array[N] matrix[2, T] log_omega;
  matrix[2, 2] A;
  {
    matrix[2, K] log_B;
    for(s in 1:2) {
      log_B[s, ] = log(B[s]');
    }
    A[1, 1] = 1 - move_probability;
    A[1, 2] = move_probability;
    A[2, 1] = 0;
    A[2, 2] = 1;
    for(i in 1:N) {
      log_omega[i] = log_B[, obs[i]];
    }
  }
}
model {
  // move_probability ~ beta(2, 2);
  {
    matrix[2, 2] log_A = log(A);
    for(i in 1:N) {
      target += loglik(log_omega[i], log_A);
    }
  }
}
generated quantities {
  // samples of change points (discard the state trajectories)
  array[N] int changepoints;
  for(i in 1:N) {
    array[T] int states = hmm_latent_rng(log_omega[i], A, rho);
    changepoints[i] = 1;
    int t = 0;
    while (t < T) {
      t += 1;
      if (states[t] == 2) {
        changepoints[i] = t;
        break;
      }
    }
  }
}
