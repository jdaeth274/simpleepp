
functions{
  
  matrix simpleepp_no_art(vector kappa, real iota, real mu, vector sigma,
                          matrix delta, vector mu_i, real dt){
    //// So in this case we have @kappa as our R parameter, @iota as our initial @mu is our population death rate, @sigma is our progression rate through CD4 stages @mu_i is our death rate for each of the cd4 stages, the length of this also determines the number of stages and finally @dt is our time step, in Euler this is one tenth of the time step
    vector[rows(kappa) + 1] S;
    matrix[rows(mu_i) - 1, rows(kappa)+1] I;
    matrix[rows(mu_i), rows(kappa)+1] D;
    vector[rows(kappa) + 1] rho;
    vector[rows(kappa) + 1] diagns;
    // rho is our I/N so when multipled by kappa it becomes our force of infection
    vector[rows(kappa)+1] lambda;         // This is our force of infection
    int DS;
    vector[rows(kappa) + 1] undiag;
    vector[rows(kappa) + 1] aids_deaths;
    DS = rows(sigma);
    
    
    // initial values
    S[1] = 1000000 * (1 - iota);
    I[1, 1] = 1000000 * iota;
    D[1, 1] = 0;
    for(m in 2:DS){
      I[m, 1] = 0;
      D[m, 1] = 0;
    }
    D[DS+1,1] = 0;
    
    rho[1] = 1 - S[1] / (S[1] + sum(I[,1]));
    lambda[1] = 0;
    
    for(t in 1:rows(kappa)) {
      real artcov;
      real deaths;
      real infections;
      real diagnoses;
      real Dt;
      real It;
      real aids_death;
      
      
      It = sum(I[,t]);
      Dt = sum(D[,t]);
      
      undiag[t] = It;
      
      lambda[t+1] = kappa[t] * rho[t];
      
      deaths = mu * (S[t] + It + Dt) + sum(mu_i[1:4] .* I[,t]) + sum(mu_i .* D[,t]);
      aids_death = mu_i[5] * D[5,t];
      aids_deaths[t] = aids_death;
      S[t+1] = S[t] + dt*(-lambda[t+1] * S[t] - mu * S[t] + deaths);                            // I think + deaths here is to 
      // keep the pop size constant
      
      I[1, t+1] = I[1, t] + dt*(lambda[t+1] * S[t] - (mu + mu_i[1] + sigma[1] + delta[1, t]) * I[1, t]);
      for(m in 2:(DS))
        I[m, t+1] = I[m, t] + dt*(sigma[m-1] * I[m-1, t] - (mu + mu_i[m] + sigma[m]+ delta[m,t]) * I[m, t]);
      
      D[1, t+1] = D[1, t] + dt*(delta[1, t] * I[1, t] - (mu + mu_i[1] + sigma[1]) * D[1, t]);
      for (m in 2:DS)
        D[m, t+1] = D[m, t] + dt*((delta[m, t] * I[m, t]) + (sigma[m-1] * D[m-1, t]) - (mu + mu_i[m] + sigma[m]) * D[m, t]);
      
      D[DS + 1, t+1] = D[DS + 1, t] + dt*((sigma[4] * I[DS, t]) + (sigma[4] * D[4, t])  - (mu + mu_i[DS + 1]) * D[DS + 1, t]);
      
      rho[t+1] = 1.0 - S[t+1] / (S[t+1] + sum(I[ ,t+1]) + sum(D[,t+1]));
      
      diagns[t] = Dt;
      if (t == rows(kappa)){
        diagns[t + 1] = sum(D[,rows(kappa) + 1]);
      }
      
    }
    
    
    return(append_col(append_col(append_row(kappa[1], kappa),
                                 append_col(lambda, rho)),append_col(diagns, append_col(undiag, aids_deaths))));
  }
  
  matrix delta_getter(vector delta_curve, vector multiplicative_vals){
    
    matrix[rows(multiplicative_vals) + 1, rows(delta_curve)] delta_out;
    int num_rows;
    row_vector[rows(delta_curve)] delta_c;
    
    
    num_rows = rows(multiplicative_vals) + 1;
    delta_c = to_row_vector(delta_curve);
    
    delta_out[1,] = delta_c;
    for (n in 2:num_rows)
      delta_out[n,] = delta_c * multiplicative_vals[ n - 1 ];
      
    return(delta_out);
  }
  
}

data {
  
  int<lower = 1> n_obs;                                                  // The number of time points we observe
  
  int<lower = 0> y[n_obs];                             // This is the number of poeple infected from our random draw
  
  int time_steps_euler;                                                  // This is our number of points over which to evaulate the function
  
  int penalty_order;                                                     // This is our first or second order penalty
  
  int time_steps_year;                                                   // This is our number of years we evaluate
  
  int estimate_period;                                                   // Number of years we predict data for
  
  matrix[time_steps_euler - 1, time_steps_year  ] X_design;                // This is our spline design matirx that we are modelling kappa with.
  
  matrix[time_steps_euler, time_steps_year] X_design_diag;
  
  matrix[time_steps_year - penalty_order , time_steps_year ] D_penalty;  // This is our penalty matrix, can be first or second order depending on the R code
  
  real mu;                                                               // This is our population level death rate
  
  vector[4] sigma;                                                       // This is our vector of transition rates
  
  vector[5] mu_i;                                                        // This is our vector of death rates
  
  
  real dt;                                                               // This is our time step
  
  real dt_2;                                                             // this is our second time step for generating the output from the fitted beta parameters
  
  int rows_to_interpret[n_obs];                             // This is a vector of the rows to use for the prevalence stats in y_hat. Corresponds to whole years
  
  int poisson_or_negbin;                                    // This variable is 0 for poisson or 1 for using the negbin dist
  
  vector[3] multiplicative_vals;
  
}



parameters{
  
  real<lower = 0> iota;                                         // The proportion of the population initially infected
  
  vector[cols(X_design)] beta;                                  // This is the knot point values
  
  real<lower = 0> sigma_pen;                                    // This is the penalty to apply to the spline to make it smooth
  
  
  real<lower= 0> phi_pen;                                       //This is the variance paramter for the neg binom model 
  
  vector<upper = 0.1>[cols(X_design)] delta_beta;                            // These are the knot values for the diagnoses curve
  
  real<lower= 0> sigma_pen_delta;                               // This is our sigma_penalty for the delta RW
  
  
  
}


transformed parameters{
  
  vector[size(rows_to_interpret)] y_hat;
  
  y_hat = simpleepp_no_art( exp(X_design * beta) , iota, mu, sigma, delta_getter((exp(X_design_diag * delta_beta)), multiplicative_vals), mu_i, dt_2)[rows_to_interpret, 4];
}


model {
  
  sigma_pen ~ normal(0, 5);
  sigma_pen_delta ~ normal(0, 10);
  delta_beta ~ normal(0.1, 3);
  
  
  beta ~ normal(0, 2.5);
  phi_pen ~ normal(0, 20);
  
  iota ~ normal(0, 0.1);                                                  // This models our initial population size
  
  target += normal_lpdf( D_penalty * beta | 0, sigma_pen);                // So this models our penalized spline with a slightly altered 
  
  target += normal_lpdf( D_penalty * delta_beta | 0, sigma_pen_delta);
  
  
  if(poisson_or_negbin == 1){
    
    target += neg_binomial_2_lpmf (y | y_hat, phi_pen);                                          // This fits it to the binomially sampled data
  }else{
    target += poisson_lpmf (y | y_hat);
    
  }
}

generated quantities{
  
  matrix[time_steps_euler , 6] fitted_output;
  vector[rows(X_design)] fitted_kappa;
  matrix[4, rows(X_design_diag)] fitted_delta;
  
  fitted_kappa = exp(X_design * beta);
  fitted_delta = delta_getter(exp(X_design_diag * delta_beta), multiplicative_vals);
  fitted_output = simpleepp_no_art(fitted_kappa, iota, mu, sigma, fitted_delta, mu_i, dt_2);       // This models the data after our fitting procedure
  
}
