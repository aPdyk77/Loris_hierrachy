// Bornean slow loris hierarchical model: primary transect random-intercept model
// Ecological state = seasonal 200-m segment use
// Observation process = positive protocol-level observation conditional on segment use

data {
  int<lower=1> N;
  int<lower=1> R;
  int<lower=1> K;
  int<lower=1, upper=N> site_id[R];
  int<lower=0, upper=1> y[R];
  int<lower=1, upper=K> transect_id[N];
  vector[N] dist_z;
  vector[R] protocol_ves;
  vector[R] visit_sept;
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_slope_sd;
  real<lower=0> prior_re_sd;
}

parameters {
  real beta0;
  real beta_dist;
  real gamma0;
  real gamma_ves;
  real gamma_sept;
  real<lower=0> sigma_transect;
  vector[K] b_raw;
}

transformed parameters {
  vector[K] b;
  vector[N] psi;
  vector[R] p;

  b = sigma_transect * b_raw;

  for (i in 1:N)
    psi[i] = inv_logit(beta0 +
                       beta_dist * dist_z[i] +
                       b[transect_id[i]]);

  for (r in 1:R)
    p[r] = inv_logit(gamma0 +
                     gamma_ves * protocol_ves[r] +
                     gamma_sept * visit_sept[r]);
}

model {
  beta0 ~ normal(0, prior_intercept_sd);
  gamma0 ~ normal(0, prior_intercept_sd);
  beta_dist ~ normal(0, prior_slope_sd);
  gamma_ves ~ normal(0, prior_slope_sd);
  gamma_sept ~ normal(0, prior_slope_sd);

  sigma_transect ~ normal(0, prior_re_sd);
  b_raw ~ std_normal();

  for (i in 1:N) {
    real log_obs;
    int any_detection;

    log_obs = 0;
    any_detection = 0;

    for (r in 1:R) {
      if (site_id[r] == i) {
        log_obs += bernoulli_lpmf(y[r] | p[r]);
        if (y[r] == 1)
          any_detection = 1;
      }
    }

    if (any_detection == 1)
      target += log(psi[i]) + log_obs;
    else
      target += log_sum_exp(log(psi[i]) + log_obs, log1m(psi[i]));
  }
}

generated quantities {
  vector[N] post_prob_use;
  int n_positive_rep;
  int n_detected_segments_rep;
  int detected_rep[N];

  n_positive_rep = 0;
  n_detected_segments_rep = 0;

  for (i in 1:N) {
    int any_detection_obs;
    real log_zero_history;
    int z_rep;

    any_detection_obs = 0;
    log_zero_history = 0;
    detected_rep[i] = 0;

    for (r in 1:R) {
      if (site_id[r] == i) {
        if (y[r] == 1)
          any_detection_obs = 1;
        log_zero_history += bernoulli_lpmf(0 | p[r]);
      }
    }

    if (any_detection_obs == 1) {
      post_prob_use[i] = 1;
    } else {
      real log_num;
      real log_den;
      log_num = log(psi[i]) + log_zero_history;
      log_den = log_sum_exp(log_num, log1m(psi[i]));
      post_prob_use[i] = exp(log_num - log_den);
    }

    z_rep = bernoulli_rng(psi[i]);

    for (r in 1:R) {
      if (site_id[r] == i) {
        int y_rep;
        if (z_rep == 1)
          y_rep = bernoulli_rng(p[r]);
        else
          y_rep = 0;

        n_positive_rep += y_rep;
        if (y_rep == 1)
          detected_rep[i] = 1;
      }
    }

    n_detected_segments_rep += detected_rep[i];
  }
}
