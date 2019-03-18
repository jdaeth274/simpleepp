###############################################################################
## Running the simple SID HIV model from STAN #################################
###############################################################################
require(reshape2)
require(splines)
require(ggplot2)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


expose_stan_functions("./simpleepp/stan_files/chunks/SID_models/SID_model.stan")

###############################################################################
## Now lets run the simulated model from the STAN functions ###################
###############################################################################



run_simulated_model<-function(params,times){
  
  mu <- params$mu                               # Non HIV mortality / exit from population
  sigma <- params$sigma                       # Progression from stages of infection
  mu_i <- params$mu_i                                   #c(0.003, 0.008, 0.035, 0.27)     # Mortality by stage, no ART
  iota<-params$iota
  delta <- params$delta
  
  
  
  dt <- times$dt                                     # time step
  nsteps <- as.integer(times$years/dt)                   # number of steps
  xstart <- times$start                                # the start of the epidemic
  step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
  xout <- c(xstart, step_vector)    
  
  rlogistic <- function(t, p) {
    p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
  }
  
  kappa_params<-c(log(params$kappa[1]),log(params$kappa[2]),params$kappa[3],params$kappa[4])
  kappa<-exp(rlogistic(step_vector,kappa_params))
  
  mod<-simpleepp_no_art(kappa = kappa, iota, mu, sigma, delta, mu_i, dt)
  sim_mod<-data.frame(mod)
  names(sim_mod)<-c("kappa","lambda","prevalence","diagnoses")
  sim_mod$prev_percent<-sim_mod$prevalence * 100
  sim_mod$time<-c(xstart,step_vector)
  
  return(list(sim_df=sim_mod,kappa_values=kappa))
  
  
}

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20, 3.1)           # Progression from stages of infection
mu_i <- c(0.003, 0.008,0.01, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
delta <- c(0.001,0.002,0.005,0.009)
iota<-0.0001

dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total



params_sim<-list(mu=mu,sigma=sigma,mu_i=mu_i,kappa=kappa,iota=iota, delta = delta)
times_sim<-list(dt=dt,years=50,start=xstart)

sim_model_output<-run_simulated_model(params_sim,times_sim)

sim_plot<-function(sim_df){
  
  prev_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=prev_percent),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Prevalence %",title="Simulated model prevalence over time")
  
  incidence_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=lambda),colour="midnightblue",size=1.05)+
    labs(x="Time",y="Incidence",title="Simulated model incidence over time")
  
  kappa_plot<-ggplot(data = sim_df)+geom_line(aes(x=time,y=kappa),colour="midnightblue",size=1.05)+
    labs(x="Time",y="kappa",title="Kappa parameter over time in simulated data")
  
  diagnoses_plot <- ggplot(data = sim_df) + 
    geom_line(aes(x=time, y = diagnoses), colour = "midnightblue", size =1.05) + 
    labs(x= "Time", y= "Diagnoses", title = "Diagnoses through time")
  
  melt_df<-sim_df[,-c(5)]
  melted_df_funky<-melt(melt_df,id="time")
  
  four_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot, incidence=incidence_plot, kappa=kappa_plot, diagnoses = diagnoses_plot,
              whole=four_variable_plot))
  
  
}

plotted_sim<-sim_plot(sim_model_output$sim_df)
plot(plotted_sim$whole)

###############################################################################
## Now that we've run the model once we can look to try and take the output ###
## of number of diagnoses and use these to fit our future models ##############
###############################################################################

year_range <- 1970:2015

data_collector <- function(sim_res, year_range, xstart, number_reps = 1){
  
  start_val <- (year_range[1] - xstart)*10 + 1
  end_val <- ((year_range[1] - xstart)*10) + ((length(year_range) - 1) * 10) + 1
  
  diagnoses_tot <- sim_res$diagnoses
  
  diag_to_keep <- seq(from  = start_val, to = end_val, by = 10)
  
  diag_to_keep <- round(diagnoses_tot[diag_to_keep])
  
  random_diags <- matrix(data = NA, nrow = number_reps,
                         ncol = length(diag_to_keep))
  
  for (i in 1:number_reps)
  random_diags[i,] <- rpois(length(diag_to_keep), diag_to_keep)
  
  return(random_diags)
}

diags_data <- data_collector(sim_model_output$sim_df, year_range, xstart)

###############################################################################
## Lets just create a large random datasets of the diagnoses to fit our #######
## cluster run models with ####################################################
###############################################################################

diags_data_to_be_saved <- data_collector(sim_model_output$sim_df,
                                         year_range, xstart,
                                         number_reps = 100)

write.table(diags_data_to_be_saved, file = "./simpleepp/analysis/SID_models_diag_10_3_2019.tsv",
            row.names = FA)

###############################################################################
## Now we've got our data from the simulated model, we can also generate the ##
## RW matrices for the input into the dataframe. ##############################
###############################################################################



xord<-seq(1970.1,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xord,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=2)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-3:45*10+1

diags_data <- diags_data[,-c(1:3)]

stan_data_discrete<-list(
  n_obs = length(diags_data),
  y = as.vector(diags_data),
  time_steps_euler = 501,
  penalty_order = 2,
  estimate_period = 5,
  time_steps_year = 51,
  X_design = spline_matrix,
  D_penalty = penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  delta = delta,
  dt = 1,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate)
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen", "phi_pen")  

test_stan_hiv<- stan("~/Dropbox/jeff_hiv_work/simpleepp/stan_files/chunks/SID_models/SID_model.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  

Sys.time()
mod_hiv_prev_pen_2 <- stan("~/Dropbox/jeff_hiv_work/simpleepp/stan_files/chunks/SID_models/SID_model.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.9))
Sys.time()

diags_df <- cbind.data.frame(diags_data, xout[rows_to_evaluate])
names(diags_df) <- c("diags","time")

plot_stan_model_fit<-function(model_output,sim_output,plot_name,xout, diags_dat){
  
  posts_hiv <- rstan::extract(model_output)
  
  
  iota_dist<-posts_hiv$iota
  params<-median(posts_hiv$iota)
  params_low<-quantile(posts_hiv$iota,c(0.025))
  params_high<-quantile(posts_hiv$iota,c(0.975))
  
  params_df<-rbind.data.frame(params_low,params,params_high)
  names(params_df)<-c("iota")
  
  sigma_pen_dist<-posts_hiv$sigma_pen
  sigma_values<-median(posts_hiv$sigma_pen)
  sigma_low<-quantile(posts_hiv$sigma_pen,c(0.025))
  sigma_high<-quantile(posts_hiv$sigma_pen,probs=c(0.975))
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
  
  
  
  prev_median<-(apply(posts_hiv$fitted_output[,,3],2,median))*100
  prev_low<-(apply(posts_hiv$fitted_output[,,3],2,quantile,probs=c(0.025)))*100
  prev_high<-(apply(posts_hiv$fitted_output[,,3],2,quantile,probs=c(0.975)))*100
  
  
  incidence_median<-apply(posts_hiv$fitted_output[,,2],2,median)
  incidence_low<-apply(posts_hiv$fitted_output[,,2],2,quantile,probs=c(0.025))
  incidence_high<-apply(posts_hiv$fitted_output[,,2],2,quantile,probs=c(0.975))
  
  r_median<-apply(posts_hiv$fitted_output[,,1],2,median)
  r_low<-apply(posts_hiv$fitted_output[,,1],2,quantile,probs=c(0.025))
  r_high<-apply(posts_hiv$fitted_output[,,1],2,quantile,probs=c(0.975))
  

  diags_median <- apply(posts_hiv$fitted_output[,,4], 2, median)
  diags_low <- apply(posts_hiv$fitted_output[,,4], 2, quantile, probs=c(0.025), na.rm = TRUE)
  diags_high <- apply(posts_hiv$fitted_output[,,4], 2, quantile, probs=c(0.975), na.rm = TRUE)
  
  # Combine into two data frames for plotting
  #df_sample = data.frame(sample_prop, sample_time)
  df_fit_prevalence = data.frame(prev_median, prev_low, prev_high, xout )
  names(df_fit_prevalence)<-c("median","low","high","time")
  df_fit_prevalence$credible_skew<-(df_fit_prevalence$high - df_fit_prevalence$median) - (df_fit_prevalence$median - df_fit_prevalence$low)
  
  df_fit_incidence<-data.frame(incidence_low,incidence_median,incidence_high,xout)
  names(df_fit_incidence)<-c("low","median","high","time")
  
  r_fit<-data.frame(r_low,r_median,r_high,xout)
  names(r_fit)<-c("low","median","high","time")
  
  diags_fit <- data.frame(diags_low, diags_median, diags_high, xout)
  names(diags_fit) <- c("low","median","high","time")
  # Plot the synthetic data with the model predictions
  # Median and 95% Credible Interval
  
  
  plotter<- ggplot(data = df_fit_prevalence) + 
    geom_line(data = df_fit_prevalence, aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(data = df_fit_prevalence,aes(x=time,ymin=low,ymax=high),
                colour="midnightblue",alpha=0.2,fill="midnightblue")+
    geom_line(data = sim_output,aes(x=time,y=prev_percent),colour="yellow",size=1)+
    coord_cartesian(xlim=c(1965,2025))+labs(x="Time",y="Prevalence (%)", title=plot_name)
  
  incidence_plot<-ggplot(data=df_fit_incidence)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="midnightblue",alpha=0.2,colour="midnightblue")+
    geom_line(data = sim_output,aes(x=time,y=lambda),colour="yellow",size=1)+
    labs(x="time",y="incidence",title="incidence_plot")
  
  r_plot<- ggplot(data = r_fit)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="midnightblue",colour="midnightblue",alpha=0.2)+
    geom_line(data = sim_output,aes(x=time,y=kappa),colour="yellow",size=1)+
  labs(x="Time",y="r value through time",title="R logistic through time")
  
  diags_plot <- ggplot(data = diags_fit) + geom_line(aes(x = time, y = median), colour = "midnightblue", size = 1)+
    geom_ribbon(aes(x = time, ymin = low, ymax = high), fill = "midnightblue", colour = "midnightblue", alpha =0.2)+
    geom_line(data = sim_output, aes(x= time, y = diagnoses), colour = "yellow", size = 1)+
    geom_point(data = diags_dat, aes(x= time, y= diags), colour = "red", size = 1)+
    labs(x = "Time", y = "diagnoses", title = "New diagnoses through time" )
  
  return(list(prevalence_plot=(plotter),df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
              r_fit_df=r_fit,incidence_plot=incidence_plot,r_plot=r_plot,sigma_pen_values=sigma_df,iota_value=params_df,
              iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist, diags_plot = diags_plot, diags_deef = diags_fit))
  
  
}


xout<-seq(1970,2020,0.1)


test_peno_2<-plot_stan_model_fit(model_output = mod_hiv_prev_pen_2,
                               plot_name = "Random walk second order",xout = xout,
                               sim_output = sim_model_output$sim_df, diags_dat = diags_df)
test_peno_2$prevalence_plot
test_peno_2$incidence_plot
test_peno_2$r_plot
test_peno_2$diags_plot

###############################################################################
## Lets do one for the RW first order methods #################################
###############################################################################

xout<-seq(1970,2019.9,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=1)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-0:45*10+1

stan_data_discrete<-list(
  n_obs = length(year_range),
  y = as.vector(diags_data),
  time_steps_euler = 501,
  penalty_order = 1,
  estimate_period = 5,
  time_steps_year = 51,
  X_design = spline_matrix,
  D_penalty = penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  delta = delta,
  dt = 1,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate)
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

test_stan_hiv<- stan("~/Dropbox/jeff_hiv_work/simpleepp/stan_files/chunks/SID_models/SID_model.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  


mod_hiv_prev_pen_1 <- stan("~/Dropbox/jeff_hiv_work/simpleepp/stan_files/chunks/SID_models/SID_model.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.99))


xout<-seq(1970,2020,0.1)


test_peno_pen_1<-plot_stan_model_fit(model_output = mod_hiv_prev_pen_1,
                               plot_name = "Random walk first order",xout = xout,
                               sim_output = sim_model_output$sim_df)
test_peno_pen_1$prevalence_plot
test_peno_pen_1$incidence_plot
test_peno_pen_1$r_plot

