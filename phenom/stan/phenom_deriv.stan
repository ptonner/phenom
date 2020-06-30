data {
  int<lower=1> N;
  int<lower=1> P; // number of replicates
  int<lower=1> K; // number of latent functions
  int<lower=1> L; // number of priors
  int<lower=1, upper=L> prior[K]; // prior assignment for each function
  real alpha_prior[L,2];
  real lengthscale_prior[L,2];
  real sigma_prior[2];
  real ls_min;
  real ls_max;

  matrix[P,K] design;
  row_vector[N] y[P];
  real x[N];
}
parameters {
  real<lower=ls_min, upper=ls_max> lengthscale[L];
  real<lower=0> alpha[L];
  real<lower=0> sigma;
  vector[N] f_eta[K];
}
transformed parameters {
  matrix[K,N] f;

  for (l in 1:L)
  {
    matrix[N, N] L_cov;
    matrix[N, N] cov;
    cov = cov_exp_quad(x, alpha[l], lengthscale[l]);
    for (n in 1:N)
      cov[n, n] += 1e-12;
    L_cov = cholesky_decompose(cov);

    for (k in 1:K)
      {
        if (prior[k] == l)
          f[k] = (L_cov * f_eta[k])';
      }
  }
}
model {

  for (l in 1:L)
  {
    lengthscale[l] ~ inv_gamma(lengthscale_prior[l,1], lengthscale_prior[l,2]);
    alpha[l] ~ gamma(alpha_prior[l,1], alpha_prior[l,2]);
  }

  sigma ~ gamma(sigma_prior[1], sigma_prior[2]);

  for (i in 1:K)
    f_eta[i] ~ normal(0, 1);

  for (i in 1:P)
    y[i] ~ normal(design[i]*f, sigma);
}
generated quantities{
  matrix[K,N] df;

  for (l in 1:L){
      vector[N] fobs;
      vector[N] df_pred;
      real lsInv = 1./lengthscale[l]/lengthscale[l];
      matrix[N, N] L_Sigma;
      vector[N] K_div_f;
      matrix[N, N] dK;
      matrix[N, N] ddK;
      matrix[N, N] v_pred;
      vector[N] df_pred_mu;
      matrix[N, N] cov_df_pred;
      matrix[N, N] nug_pred;
      matrix[N, N] Sigma;
      real diff;

      nug_pred = diag_matrix(rep_vector(1e-8,N));

      // cov(f)
      Sigma = cov_exp_quad(x, alpha[l], lengthscale[l]);
      for (n in 1:N)
        Sigma[n, n] += 1e-8;

      // prepare cholesky for operations
      L_Sigma = cholesky_decompose(Sigma);

      // compute dK: cov(df, f), and ddK: cov(df, df)
      dK = cov_exp_quad(x, alpha[l], lengthscale[l]);
      ddK = cov_exp_quad(x, alpha[l], lengthscale[l]);
      for (i in 1:N){
        for (j in 1:N){
          diff = x[i] - x[j];

          dK[i,j] *= (-lsInv * diff);
          ddK[i,j] *= (1.-lsInv*diff*diff) * lsInv;
        }
      }

      //compute df/dt for functions in this prior
      for (k in 1:K)
        {
          if (prior[k] == l)
            {
              fobs = f[k]';

              // solve for Sigma^{-1} f
              K_div_f = mdivide_left_tri_low(L_Sigma, fobs);
              K_div_f = mdivide_right_tri_low(K_div_f',L_Sigma)';

              df_pred_mu = (dK * K_div_f);

              v_pred = mdivide_left_tri_low(L_Sigma, dK');
              cov_df_pred = ddK - v_pred' * v_pred;

              df[k] = multi_normal_rng(df_pred_mu, cov_df_pred + nug_pred)';
            }
        }
  }
}
