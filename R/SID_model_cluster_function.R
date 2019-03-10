###############################################################################
## Function for the SID model on the cluster ##################################
###############################################################################

fitting_SID_model <- function(samples_data_frame,iteration_number,data_about_sampling,params,simulated_true_df){
  
  rstan::rstan_options(auto_write = TRUE)                               ## Need these options like this for parrallel going
  options(mc.cores = parallel::detectCores())
  
  data_about_sampling$check_rhat <- TRUE
  
  
  xout<-seq(1970,2019.9,0.1)
  spline_matrix<-splines::splineDesign(1969:2021,xout,ord = 2)                                   ## This matrix is the spline design one 
  penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=data_about_sampling$penalty_order)        ## This matrix creates the differences between your kappa values 
  
  plot_stan_model_fit<-function(model_output,sim_sample,sim_output,plot_name,xout){
    
    posts_hiv <- rstan::extract(model_output)
    
    
    iota_dist<-posts_hiv$iota
    params<-stats::median(posts_hiv$iota)
    params_low<-stats::quantile(posts_hiv$iota,c(0.025))
    params_high<-stats::quantile(posts_hiv$iota,c(0.975))
    
    params_df<-rbind.data.frame(params_low,params,params_high)
    names(params_df)<-c("iota")
    
    sigma_pen_dist<-posts_hiv$sigma_pen
    sigma_values<-stats::median(posts_hiv$sigma_pen)
    sigma_low<-stats::quantile(posts_hiv$sigma_pen,c(0.025))
    sigma_high<-stats::quantile(posts_hiv$sigma_pen,probs=c(0.975))
    sigma_df<-rbind.data.frame(sigma_low,sigma_values,sigma_high)
    names(sigma_df)<-c("sigma_pen")
    
    
    
    # These should match well. 
    
    #################
    # Plot model fit:
    
    # Proportion infected from the synthetic data:
    
    #sample_prop = sample_y / sample_n
    
    # Model predictions across the sampling time period.
    # These were generated with the "fake" data and time series.
    #mod_median = apply(posts_hiv$fake_I[,,2], 2, median)
    #mod_low = apply(posts_hiv$fake_I[,,2], 2, quantile, probs=c(0.025))
    #mod_high = apply(posts_hiv$fake_I[,,2], 2, quantile, probs=c(0.975))
    mod_time = xout
    
    prev_median<-(apply(posts_hiv$fitted_output[,,3],2,stats::median))*100
    prev_low<-(apply(posts_hiv$fitted_output[,,3],2,stats::quantile,probs=c(0.025)))*100
    prev_high<-(apply(posts_hiv$fitted_output[,,3],2,stats::quantile,probs=c(0.975)))*100
    
    
    incidence_median<-apply(posts_hiv$fitted_output[,,2],2,stats::median)
    incidence_low<-apply(posts_hiv$fitted_output[,,2],2,stats::quantile,probs=c(0.025))
    incidence_high<-apply(posts_hiv$fitted_output[,,2],2,stats::quantile,probs=c(0.975))
    
    r_median<-apply(posts_hiv$fitted_output[,,1],2,stats::median)
    r_low<-apply(posts_hiv$fitted_output[,,1],2,stats::quantile,probs=c(0.025))
    r_high<-apply(posts_hiv$fitted_output[,,1],2,stats::quantile,probs=c(0.975))
    
    diag_median<-apply(posts_hiv$fitted_output[,,4],2,stats::median)
    diag_low<-apply(posts_hiv$fitted_output[,,4],2,stats::quantile,probs=c(0.025))
    diag_high<-apply(posts_hiv$fitted_output[,,4],2,stats::quantile,probs=c(0.975))
    
    
    # Combine into two data frames for plotting
    #df_sample = data.frame(sample_prop, sample_time)
    df_fit_prevalence = data.frame(prev_median, prev_low, prev_high, xout )
    names(df_fit_prevalence)<-c("median","low","high","time")
    df_fit_prevalence$credible_skew<-(df_fit_prevalence$high - df_fit_prevalence$median) - (df_fit_prevalence$median - df_fit_prevalence$low)
    
    df_fit_incidence<-data.frame(incidence_low,incidence_median,incidence_high,xout)
    names(df_fit_incidence)<-c("low","median","high","time")
    
    r_fit<-data.frame(r_low,r_median,r_high,xout)
    names(r_fit)<-c("low","median","high","time")
    
    diag_fit <- data.frame(diag_low, diag_median, diag_high, xout)
    names(diag_fit) <- c("low","median","high","time")
    # Plot the synthetic data with the model predictions
    # Median and 95% Credible Interval
    
    
    
    return(list(df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
                r_fit_df=r_fit,sigma_pen_values=sigma_df,iota_value=params_df,
                iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist, diag_df = diag_fit))
    
    
  }
  
  prev_df_tot <- NULL
  incidence_df_tot <- NULL
  kappa_df_tot<-NULL
  diagnoses_df_tot <- NULL
  iota_values_tot<-NULL
  tot_fit_summary <- NULL
  sample_start<-samples_data_frame$sample_time_hiv[1] - 1970
  tot_rhat_check <- cbind(c(1:100),rep(0,100))
  
  for( i in 1:iteration_number){
    
    
    diags_data<-samples_data_frame[i,]
    
    
    stan_data_discrete<-list(
      n_obs = length(data_about_sampling$sample_years),
      y = as.array(diags_data),
      time_steps_euler = 501,
      penalty_order = data_about_sampling$penalty_order,
      estimate_period = 5,
      time_steps_year = 51,
      X_design = spline_matrix,
      D_penalty = penalty_matrix,
      mu = data_about_sampling$mu,
      sigma = data_about_sampling$sigma,
      mu_i = data_about_sampling$mu_i,
      delta = data_about_sampling$delta,
      dt = 1,
      dt_2 = 0.1,
      rows_to_interpret = as.array(data_about_sampling$rows_to_evaluate)
    )
    
    params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  
    
    
    mod_hiv_prev <- rstan::stan("simpleepp/stan_files/chunks/SID_model/SID_model.stan", data = stan_data_discrete,
                                pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                                control = list(adapt_delta = 0.99))
    
    if(data_about_sampling$check_rhat == T){
      fit_summary <- summary(mod_hiv_prev, probs = c(0.5))$summary
      N <- dim(fit_summary)[[1]]
      for (n in 1:N) {
        rhat <- fit_summary[,6][n]
        if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
          a<-sprintf('Rhat for parameter %s is %s!',rownames(fit_summary)[n], rhat)
          no_warning <- FALSE
        }
      }
      tot_rhat_check[i,2] <- a
    }
    fit_summary <- rstan::summary(mod_hiv_prev)$summary
    fit_summary <- data.frame(fit_summary)
    fit_summary$iter <- rep(i , nrow(fit_summary))
    
    xout<-seq(1970,2020,0.1)
    
    
    stan_output_random_walk_second_n_100<-plot_stan_model_fit(model_output = mod_hiv_prev,sim_sample = diags_data,
                                                              plot_name = "Random walk second order, n = 100",xout = xout,
                                                              sim_output = simulated_true_df)
    
    prev_df<-stan_output_random_walk_second_n_100$df_output
    prev_df$iteration<-rep(i,nrow(prev_df))
    
    incidence_df<-stan_output_random_walk_second_n_100$incidence_df
    incidence_df$iteration<-rep(i,nrow(incidence_df))
    
    
    kappa_df<-stan_output_random_walk_second_n_100$r_fit_df
    kappa_df$iteration<-rep(i,nrow(kappa_df))
    
    
    iota_value<-stan_output_random_walk_second_n_100$iota_value
    iota_value$iteration<-rep(i,nrow(iota_value))
    
    diag_df <- stan_output_random_walk_second_n_100$diag_df
    diag_df$iteration <- rep(i, nrow(diag_df))
    
    
    prev_df_tot<-rbind(prev_df_tot,prev_df)
    incidence_df_tot <- rbind(incidence_df_tot,incidence_df)
    kappa_df_tot <- rbind(kappa_df_tot,kappa_df)
    diagnoses_df_tot <- rbind(diagnoses_df_tot, diag_df)
    iota_values_tot <- rbind(iota_values_tot,iota_value)
    tot_fit_summary <- rbind.data.frame(tot_fit_summary, fit_summary)
    
    print((i/iteration_number)*100)
    
  }
  
  data_about_sampling$simul_type<-"Random_Walk"
  
  return(list(prev=prev_df_tot,incidence=incidence_df_tot,kappa=kappa_df_tot,iota_values=iota_values_tot,
              data_about_run = data_about_sampling, rhat_check = tot_rhat_check, fit_sumo = tot_fit_summary,
              diagnoses_df = diagnoses_df_tot))
  
}
