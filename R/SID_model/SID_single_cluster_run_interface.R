###############################################################################
## Setting up the SID model for the DIDE cluster ##############################
###############################################################################
require(reshape2)
require(splines)
require(ggplot2)
require(rstan)
require(ggpubr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/HOMES_drive/")
options(didehpc.username = "jd2117",didehpc.home = "~/HOMES_drive/simpleepp",didehpc.cluster = "fi--didemrchnb")

didehpc::didehpc_config(cores = 3,parallel = FALSE)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## !!!!!!!!!!!!!!!!!!!!!!!! Remember to turn on pulse secure at this point !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
didehpc::didehpc_config()

context::context_log_start()

root <- "contexts_4"

ctx<- context::context_save(root,packages = c("rstan","ggplot2","splines"),
                            sources = "./simpleepp/R/SID_model/SID_model_single_run_delta_rw.R")
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)

################################################################################################################################
## So the above lines of code give me access to the cluster, with the obj object giving me queue functionalaity, I've also #####
## asked for three cores for stan to use for each of the chains ################################################################
################################################################################################################################
obj$cluster_load()
obj$task_status()

###################################################################################################################################
## Now we will run through the data from begininng ################################################################################
###################################################################################################################################

expose_stan_functions("~/HOMES_drive/simpleepp/stan_files/chunks/SID_model/SID_model_RW_delta.stan")

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
  
  delta_ps <- c(params$delta_parms[1],params$delta_parms[2], params$delta_parms[3],
                params$delta_parms[4])
  delta_step <-seq(xstart, by = dt, length.out = nsteps + 1)
  delta_single <- rlogistic(delta_step, delta_ps)
  
  delta <- matrix(data = 0, nrow = 4, ncol = length(kappa) + 1)
  delta[1,] <- delta_single
  delta[2,] <- delta_single * 1.02
  delta[3,] <- delta_single * 1.04
  delta[4,] <- delta_single * 1.1
  
  plot(delta[4,], type="l",col="red")
  lines(delta[3,], col = "blue")
  lines(delta[2,], col="green")
  lines(delta[1,], col="cyan")
  
  
  
  mod<-simpleepp_no_art(kappa = kappa, iota, mu, sigma, delta, mu_i, dt)
  sim_mod<-data.frame(mod)
  names(sim_mod)<-c("kappa","lambda","prevalence","diagnoses","undiagnosed","aids_deaths")
  sim_mod$prev_percent<-sim_mod$prevalence * 100
  sim_mod$time<-c(xstart,step_vector)
  
  return(list(sim_df=sim_mod,kappa_values=kappa, delta_vals = delta))
  
  
}

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20, 3.1)           # Progression from stages of infection
mu_i <- c(0.003, 0.008,0.01, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
delta_parms <- c(0.05,0.9,0.4,1985)
iota<-0.0001

dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total



params_sim<-list(mu=mu,sigma=sigma,mu_i=mu_i,kappa=kappa,iota=iota, delta_parms = delta_parms)
times_sim<-list(dt=dt,years=50,start=xstart)

sim_model_output<-run_simulated_model(params_sim,times_sim)

delta <- sim_model_output$delta_vals

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
  
  undiag_plot <- ggplot(data = sim_df) +
    geom_line(aes(x=time,y=undiagnosed), colour = "midnightblue", size = 1.05) +
    labs(x = "Time", y = "Number undiagnosed", title = "Undiagnosed through time")
  
  aids_deaths_plot <- ggplot(data = sim_df) +
    geom_line(aes(x = time, y = aids_deaths), colour = "midnightblue", size = 1.05) + 
    labs(x = "Time", y = "AIDS deaths", title = "AIDS deaths through time")
  
  melt_df<-sim_df[,-c(4,5,6,7)]
  melted_df_funky<-melt(melt_df,id="time")
  
  four_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot, incidence=incidence_plot, kappa=kappa_plot, diagnoses = diagnoses_plot,
              whole=four_variable_plot, undiagnosed = undiag_plot, aids_deaths = aids_deaths_plot))
  
  
}

plotted_sim<-sim_plot(sim_model_output$sim_df)
plot(plotted_sim$whole)
plot(plotted_sim$undiagnosed)
plot(plotted_sim$aids_deaths)

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
xord_diag <- seq(1970,2015,0.1)
spline_matrix<- splineDesign(1969:2016,xord_diag,ord = 2)


###############################################################################
## I think if we take our logistic curve, then use the rnorm function to ######
## generate some random noise around this shape for us then to form the log ###
## curve with #################################################################
###############################################################################
xstart = 1970
dt = 0.1
nsteps = 451

rlogistic <- function(t, p) {
  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}

delta_ps <- c(0.05,0.9,0.4,1985)

delta_step <-seq(xstart, by = dt, length.out = nsteps)
delta_single <- rlogistic(delta_step, delta_ps)

for (i in 1:length(delta_single)){
  delta_single[i] <- rnorm(1, mean = delta_single[i], sd = (delta_single[i]/8))
  if( delta_single[i] > 1){
    delta_single[i] <- 1
  }
}
plot(delta_single, type = "l")

yearly_delts <- delta_single[(0:45)*10 + 1]

rw_random <- spline_matrix %*% yearly_delts
plot(rw_random, type = "l")
###############################################################################
## Now lets see what the epidemics look like with this more random value ######
## for delta ##################################################################
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
  
  delta_single <- params$delta_curve
  
  delta <- matrix(data = 0, nrow = 4, ncol = length(kappa) + 1)
  delta[1,] <- delta_single
  for (j in 2:4){
    for(k in 1:length(delta_single)){
      delta[j,k] <- delta_single[k] * params$multi_vals[j-1]
      if (delta[j,k] > 1){
        delta[j, k] <- 1
      }
    }  
  }
  
  plot(delta[4,], type="l",col="red")
  lines(delta[3,], col = "blue")
  lines(delta[2,], col="green")
  lines(delta[1,], col="cyan")
  
  
  
  mod<-simpleepp_no_art(kappa = kappa, iota, mu, sigma, delta, mu_i, dt)
  sim_mod<-data.frame(mod)
  names(sim_mod)<-c("kappa","lambda","prevalence","diagnoses","undiagnosed","aids_deaths")
  sim_mod$prev_percent<-sim_mod$prevalence * 100
  sim_mod$time<-c(xstart,step_vector)
  
  return(list(sim_df=sim_mod,kappa_values=kappa, delta_vals = delta))
  
  
}

mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20, 3.1)           # Progression from stages of infection
mu_i <- c(0.003, 0.008,0.01, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
delta_parms <- c(0.05,0.9,0.4,1985)
iota<-0.0001
multiplicative_vals = c(1.02,1.04,1.1)
dt <- 0.1                                     # time step
nsteps <- as.integer(45/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total



params_sim<-list(mu=mu,sigma=sigma,mu_i=mu_i,kappa=kappa,iota=iota, delta_curve = rw_random,
                 multi_vals = multiplicative_vals)
times_sim<-list(dt=dt,years=45,start=xstart)

sim_model_output<-run_simulated_model(params_sim,times_sim)

delta <- sim_model_output$delta_vals

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
  
  max_diag <- max(sim_df$diagnoses)
  sim_df$diag_shape <- sim_df$diagnoses/max_diag
  
  
  
  undiag_plot <- ggplot(data = sim_df) +
    geom_line(aes(x=time,y=undiagnosed), colour = "midnightblue", size = 1.05) +
    labs(x = "Time", y = "Number undiagnosed", title = "Undiagnosed through time")
  
  aids_deaths_plot <- ggplot(data = sim_df) +
    geom_line(aes(x = time, y = aids_deaths), colour = "midnightblue", size = 1.05) + 
    labs(x = "Time", y = "AIDS deaths", title = "AIDS deaths through time")
  
  melt_df<-sim_df[,-c(4,5,6,7)]
  melted_df_funky<-melt(melt_df,id="time")
  
  four_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot, incidence=incidence_plot, kappa=kappa_plot, diagnoses = diagnoses_plot,
              whole=four_variable_plot, undiagnosed = undiag_plot, aids_deaths = aids_deaths_plot))
  
  
}

plotted_sim<-sim_plot(sim_model_output$sim_df)
plot(plotted_sim$whole)
plot(plotted_sim$diagnoses)
plot(plotted_sim$undiagnosed)
plot(plotted_sim$aids_deaths)

###############################################################################
## Now lets get our diagnoses data from this simulation #######################
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
## Now we've got our data from the simulated model, we can also generate the ##
## RW matrices for the input into the dataframe. ##############################
###############################################################################



xord<-seq(1970.1,2015,0.1)
xord_diag <- seq(1970,2015,0.1)
spline_matrix<-splineDesign(1969:2016,xord,ord = 2)            ## This matrix is the spline design one 
spline_mat_diag <- splineDesign(1969:2016,xord_diag,ord = 2)
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=2)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-3:45*10+1

diags_data <- diags_data[,-c(1:3)]

stan_data_discrete<-list(
  n_obs = length(diags_data),
  y = as.vector(diags_data),
  time_steps_euler = 451,
  penalty_order = 2,
  estimate_period = 0,
  time_steps_year = 46,
  X_design = spline_matrix,
  X_design_diag = spline_mat_diag,
  D_penalty = penalty_matrix,
  mu = mu,
  sigma = sigma,
  mu_i = mu_i,
  delta = delta,
  dt = 1,
  dt_2 = 0.1,
  rows_to_interpret = as.array(rows_to_evaluate),
  poisson_or_negbin = 1,
  multiplicative_vals = c(1.02,1.04,1.1),
  delta_priors = c(0,5,0.1,3),
  kappa_priors = c(0,5,0,2.5)
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen", "phi_pen", "delta_beta")  

test_stan_hiv<- stan("~/HOMES_drive/simpleepp/stan_files/chunks/SID_model/SID_model_RW_delta.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  

mod_hiv_prev_pen_2 <- obj$enqueue(single_run_rw_delta(stan_data = stan_data_discrete,
                           parameters_to_mon = params_monitor_hiv,
                           warmup_length = 500,iter_length = 1500,
                           adapt_delta = 0.95), name = "single_run_rw")
mod_hiv_prev_pen_2$status()
mod_hiv_prev_pen_2$log()

single_run_res <- mod_hiv_prev_pen_2$result()

###############################################################################
## Now lets analyse the results ###############################################
###############################################################################

diags_df <- cbind.data.frame(diags_data, xout[rows_to_evaluate])
names(diags_df) <- c("diags","time")

plot_stan_model_fit<-function(model_output,sim_output,plot_name,xout, diags_dat, X_Design){
  
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
  
  phi_vals_median <- median(posts_hiv$phi_pen)
  phi_low<-quantile(posts_hiv$phi_pen,c(0.025))
  phi_high<-quantile(posts_hiv$phi_pen,probs=c(0.975))
  phi_df<-rbind.data.frame(phi_low,phi_vals_median,phi_high)
  names(phi_df)<-c("phi_pen")
  
  beta_vals_median <- apply(posts_hiv$beta, 2, median)
  beta_vals_low <- apply(posts_hiv$beta,2, quantile, probs=(0.025))
  beta_vals_high <- apply(posts_hiv$beta, 2, quantile, probs=(0.975))
  beta_df<-cbind.data.frame(beta_vals_low,beta_vals_median,beta_vals_high)
  names(beta_df)<-c("low","median","high")
  
  diag_beta_median <-  X_Design %*%  exp(apply(posts_hiv$delta_beta, 2, median)) 
  diag_beta_low <- X_Design %*% exp(apply(posts_hiv$delta_beta, 2, quantile, probs=(0.025))) 
  diag_beta_high <- X_Design %*% exp(apply(posts_hiv$delta_beta, 2, quantile, probs=(0.975))) 
  diag_beta_df <- cbind.data.frame(diag_beta_low, diag_beta_median, diag_beta_high)
  names(diag_beta_df) <- c("low","median","high")
  
  
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
  
  undiag_median <- apply(posts_hiv$fitted_output[,,5], 2, median)
  undiag_low <- apply(posts_hiv$fitted_output[,,5], 2, quantile, probs=c(0.025), na.rm = TRUE)
  undiag_high <- apply(posts_hiv$fitted_output[,,5], 2, quantile, probs=c(0.975), na.rm = TRUE)
  
  aids_median <- apply(posts_hiv$fitted_output[,,6], 2, median)
  aids_low <- apply(posts_hiv$fitted_output[,,6], 2, quantile, probs=c(0.025), na.rm = TRUE)
  aids_high <- apply(posts_hiv$fitted_output[,,6], 2, quantile, probs=c(0.975), na.rm = TRUE)
  
  
  # Combine into two data frames for plotting
  #df_sample = data.frame(sample_prop, sample_time)
   df_fit_prevalence = data.frame(prev_median, prev_low, prev_high, xout )
   names(df_fit_prevalence)<-c("median","low","high","time")
   df_fit_prevalence$credible_skew<-(df_fit_prevalence$high - df_fit_prevalence$median) - (df_fit_prevalence$median - df_fit_prevalence$low)
  # 
  df_fit_incidence<-data.frame(incidence_low,incidence_median,incidence_high,xout)
  names(df_fit_incidence)<-c("low","median","high","time")
  
  r_fit<-data.frame(r_low,r_median,r_high,xout)
  names(r_fit)<-c("low","median","high","time")
  
  diags_fit <- data.frame(diags_low, diags_median, diags_high, xout)
  names(diags_fit) <- c("low","median","high","time")
  
  undiag_fit <- data.frame(undiag_low, undiag_median, undiag_high, xout)
  names(undiag_fit) <- c("low","median","high","time")
  
  aids_fit <- data.frame(aids_low, aids_median, aids_high, xout)
  names(aids_fit) <- c("low","median","high","time")
  
  # 
  diag_beta_df$time <- xout
  # 
  real_diag_rates <- data.frame(t(sim_output$delta_vals))
  names(real_diag_rates) <- "orig_curve"
  real_diag_rates$time <- xout
  
  
   plotter<- ggplot(data = df_fit_prevalence) + 
     geom_line(data = df_fit_prevalence, aes(x=time,y=median),colour="midnightblue",size=1)+
     geom_ribbon(data = df_fit_prevalence,aes(x=time,ymin=low,ymax=high),
                 colour="midnightblue",alpha=0.2,fill="midnightblue")+
     geom_line(data = sim_output$sim_df,aes(x=time,y=prev_percent),colour="yellow",size=1)+
     coord_cartesian(xlim=c(1965,2025))+labs(x="Time",y="Prevalence (%)", title=plot_name)
   
  incidence_plot<-ggplot(data=df_fit_incidence)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="midnightblue",alpha=0.2,colour="midnightblue")+
    geom_line(data = sim_output$sim_df,aes(x=time,y=lambda),colour="yellow",size=1)+
    labs(x="time",y="incidence",title="incidence_plot")
  
  r_plot<- ggplot(data = r_fit)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="midnightblue",colour="midnightblue",alpha=0.2)+
    geom_line(data = sim_output$sim_df,aes(x=time,y=kappa),colour="yellow",size=1)+
    labs(x="Time",y="r value through time",title="R logistic through time")
  
  diags_plot <- ggplot(data = diags_fit) + geom_line(aes(x = time, y = median), colour = "midnightblue", size = 1)+
    geom_ribbon(aes(x = time, ymin = low, ymax = high), fill = "midnightblue", colour = "midnightblue", alpha =0.2)+
    geom_line(data = sim_output$sim_df, aes(x= time, y = diagnoses), colour = "yellow", size = 1)+
    geom_point(data = diags_dat, aes(x= time, y= diags), colour = "red", size = 1)+
    labs(x = "Time", y = "diagnoses", title = "New diagnoses through time" )
  
  diag_rate_plot <- ggplot(data = diag_beta_df) + geom_line(aes(x = time, y = median), colour = "midnightblue", size = 1) +
    geom_ribbon(aes(x = time, ymin = low, ymax = high), fill = "midnightblue", colour = "midnightblue", alpha = 0.2)+
    geom_line(data = real_diag_rates, aes(x = time, y = orig_curve), colour = "yellow", size = 1) 
  
  undiag_plot <- ggplot(data = undiag_fit)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="midnightblue",colour="midnightblue",alpha=0.2)+
    geom_line(data = sim_output$sim_df,aes(x=time,y=undiagnosed),colour="yellow",size=1)+
    labs(x="Time",y="Number undiagnosed",title="Undiagnosed through time")
  
  aids_deaths_plot <- ggplot(data = aids_fit)+geom_line(aes(x=time,y=median),colour="midnightblue",size=1)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="midnightblue",colour="midnightblue",alpha=0.2)+
    geom_line(data = sim_output$sim_df,aes(x=time,y=aids_deaths),colour="yellow",size=1)+
    labs(x="Time",y="AIDS deaths",title="AIDS deaths through time")
  
  
  combed_plot <- ggarrange(plotter, incidence_plot, r_plot, diags_plot, diag_rate_plot,
                           undiag_plot, aids_deaths_plot, ncol = 2, nrow = 4)
  
  
  return(list(prevalence_plot=(plotter),df_output=df_fit_prevalence,incidence_df=df_fit_incidence,
              r_fit_df=r_fit,incidence_plot=incidence_plot,r_plot=r_plot,sigma_pen_values=sigma_df,iota_value=params_df,
              iota_dist=iota_dist,sigma_pen_dist=sigma_pen_dist, diags_plot = diags_plot, diags_deef = diags_fit,
              phi_vals = phi_df, beta_df = beta_df, combed_plot = combed_plot, diag_rate_plot = diag_rate_plot,
              undiag_plot = undiag_plot, aids_deaths_plot = aids_deaths_plot))
  
  
}


xout<-seq(1970,2015,0.1)


test_peno_8<-plot_stan_model_fit(model_output = single_run_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)

test_peno_5$diag_rate_plot
test_peno_5$prevalence_plot
test_peno_5$incidence_plot
test_peno_5$r_plot
test_peno_5$diags_plot
test_peno_5$diags_deef

test_peno_6$diag_rate_plot
test_peno_6$prevalence_plot
test_peno_6$incidence_plot
test_peno_6$r_plot
test_peno_6$diags_plot
test_peno_6$diags_deef

test_peno_7$diag_rate_plot
test_peno_7$prevalence_plot
test_peno_7$incidence_plot
test_peno_7$r_plot
test_peno_7$diags_plot
test_peno_7$diags_deef
test_peno_7$undiag_plot
test_peno_7$aids_deaths_plot

test_peno_8$diag_rate_plot
test_peno_8$prevalence_plot
test_peno_8$incidence_plot
test_peno_8$r_plot
test_peno_8$diags_plot
test_peno_8$diags_deef
test_peno_8$undiag_plot
test_peno_8$aids_deaths_plot


###############################################################################
## Now we've added in the functionality to name the priors directly, lets run #
## through the priors iteratively #############################################
###############################################################################
obj <- didehpc::queue_didehpc(ctx,config)


run_ids_2 <- NULL
prior_vals <- seq(from = 2, to = 4, length.out = 5)

for(k in 1:5){
  j <- prior_vals[k]
  delta_priors <- c(0,(j*2),0.1,j)
  kappa_priors <- c(0, (j*2), 0, j)
  
  stan_data_discrete<-list(
    n_obs = length(diags_data),
    y = as.vector(diags_data),
    time_steps_euler = 451,
    penalty_order = 2,
    estimate_period = 0,
    time_steps_year = 46,
    X_design = spline_matrix,
    X_design_diag = spline_mat_diag,
    D_penalty = penalty_matrix,
    mu = mu,
    sigma = sigma,
    mu_i = mu_i,
    delta = delta,
    dt = 1,
    dt_2 = 0.1,
    rows_to_interpret = as.array(rows_to_evaluate),
    poisson_or_negbin = 1,
    multiplicative_vals = c(1.02,1.04,1.1),
    delta_priors = delta_priors,
    kappa_priors = kappa_priors
  )
  
  name_run <- paste("Prior_combo:",as.character(j))
  
  params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen", "phi_pen", "delta_beta")
  
  cluster_run <- obj$enqueue(single_run_rw_delta(stan_data = stan_data_discrete,
                                                        parameters_to_mon = params_monitor_hiv,
                                                        warmup_length = 500,iter_length = 1500,
                                                        adapt_delta = 0.95), name = name_run)
  
  print(cluster_run$status())
  run_ids_2 <- c(run_ids_2, cluster_run$id)
  run_name <- paste("cluster_run",j,sep = "_")
  
  

}

prior_1_combo <- obj$task_get(run_ids[1])
prior_2_combo <- obj$task_get(run_ids[2])
prior_3_combo <- obj$task_get(run_ids[3])
prior_4_combo <- obj$task_get(run_ids[4])
prior_5_combo <- obj$task_get(run_ids[5])

prior_1_combo$status()
prior_2_combo$status()
prior_3_combo$status()
prior_4_combo$status()
prior_5_combo$status()

prior_1_combo$log()
prior_2_combo$log()
prior_3_combo$log()
prior_4_combo$log()
prior_5_combo$log()

prior_1_combo_res <- prior_1_combo$result()
prior_2_combo_res <- prior_2_combo$result()
prior_3_combo_res <- prior_3_combo$result()
prior_4_combo_res <- prior_4_combo$result()
prior_5_combo_res <- prior_5_combo$result()

prior_res_2<-plot_stan_model_fit(model_output = prior_2_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)

prior_res_3<-plot_stan_model_fit(model_output = prior_3_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)
prior_res_4<-plot_stan_model_fit(model_output = prior_4_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)
prior_res_5<-plot_stan_model_fit(model_output = prior_5_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)

prior_res_2$combed_plot
prior_res_3$combed_plot
prior_res_4$combed_plot
prior_res_5$combed_plot


###############################################################################
## So between 3 and 4 looks the best bet, lets look in more detail at these ###
## results ####################################################################
###############################################################################

prior_2_combo <- obj$task_get(run_ids_2[1])
prior_2.5_combo <- obj$task_get(run_ids_2[2])
prior_3_combo <- obj$task_get(run_ids_2[3])
prior_3.5_combo <- obj$task_get(run_ids_2[4])
prior_4_combo <- obj$task_get(run_ids_2[5])

prior_2_combo$status()
prior_2.5_combo$status()
prior_3_combo$status()
prior_3.5_combo$status()
prior_4_combo$status()

prior_2_combo$log()
prior_2.5_combo$log()
prior_3_combo$log()
prior_3.5_combo$log()
prior_4_combo$log()

prior_2_combo_res <- prior_2_combo$result()
prior_2.5_combo_res <- prior_2.5_combo$result()
prior_3_combo_res <- prior_3_combo$result()
prior_3.5_combo_res <- prior_3.5_combo$result()
prior_4_combo_res <- prior_4_combo$result()

prior_res_2<-plot_stan_model_fit(model_output = prior_2_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)

prior_res_2.5<-plot_stan_model_fit(model_output = prior_2.5_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)
prior_res_3<-plot_stan_model_fit(model_output = prior_3_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)
prior_res_3.5<-plot_stan_model_fit(model_output = prior_3.5_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)

prior_res_4<-plot_stan_model_fit(model_output = prior_4_combo_res,
                                   plot_name = "Random walk second order",xout = xout,
                                   sim_output = sim_model_output, diags_dat = diags_df,
                                   X_Design = stan_data_discrete$X_design_diag)


prior_res_2$combed_plot
prior_res_2.5$combed_plot
prior_res_3$combed_plot
prior_res_5$combed_plot


###############################################################################
## So the most promising of those results were between 2 and 3, lets try ######
## out the difference between these two values ################################
###############################################################################

obj <- didehpc::queue_didehpc(ctx,config)

testing_priors <- function(stan_data, prior_start, prior_end,
                           length_o_priors){

  run_ids_3 <- NULL
  prior_vals <- seq(from = prior_start, to = prior_end,
                    length.out = length_o_priors)
  
  for(k in 1:5){
    j <- prior_vals[k]
    delta_priors <- c(0,(j*2),0.1,j)
    kappa_priors <- c(0, (j*2), 0, j)
    
    stan_data_discrete <- stan_data
    stan_data_discrete$delta_priors <- delta_priors
    stan_data_discrete$kappa_priors <- kappa_priors
    
    name_run <- paste("Prior_combo:",as.character(j))
    
    params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen", "phi_pen", "delta_beta")
    
    cluster_run <- obj$enqueue(single_run_rw_delta(stan_data = stan_data_discrete,
                                                   parameters_to_mon = params_monitor_hiv,
                                                   warmup_length = 500,iter_length = 1500,
                                                   adapt_delta = 0.95), name = name_run)
    
    print(cluster_run$status())
    run_ids_3 <- c(run_ids_3, cluster_run$id)
    
    
    
  }
  
  return(run_ids_3)
}

run_2_3_outs <- testing_priors(stan_data_discrete, prior_start = 2,
                               prior_end = 3, length_o_priors = 5)

prior_2_combo <- obj$task_get(run_2_3_outs[1])
prior_2.25_combo <- obj$task_get(run_2_3_outs[2])
prior_2.5_combo <- obj$task_get(run_2_3_outs[3])
prior_2.75_combo <- obj$task_get(run_2_3_outs[4])
prior_3_combo <- obj$task_get(run_2_3_outs[5])

prior_2_combo$status()
prior_2.25_combo$status()
prior_2.5_combo$status()
prior_2.75_combo$status()
prior_3_combo$status()

prior_2_combo$log()
prior_2.25_combo$log()
prior_2.5_combo$log()
prior_2.75_combo$log()
prior_3_combo$log()

prior_2_combo_res <- prior_2_combo$result()
prior_2.25_combo_res <- prior_2.25_combo$result()
prior_2.5_combo_res <- prior_2.5_combo$result()
prior_2.75_combo_res <- prior_2.75_combo$result()
prior_3_combo_res <- prior_3_combo$result()

prior_res_2<-plot_stan_model_fit(model_output = prior_2_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)

prior_res_2.25<-plot_stan_model_fit(model_output = prior_2.25_combo_res,
                                   plot_name = "Random walk second order",xout = xout,
                                   sim_output = sim_model_output, diags_dat = diags_df,
                                   X_Design = stan_data_discrete$X_design_diag)
prior_res_2.5<-plot_stan_model_fit(model_output = prior_2.5_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)
prior_res_2.75<-plot_stan_model_fit(model_output = prior_2.75_combo_res,
                                   plot_name = "Random walk second order",xout = xout,
                                   sim_output = sim_model_output, diags_dat = diags_df,
                                   X_Design = stan_data_discrete$X_design_diag)

prior_res_3<-plot_stan_model_fit(model_output = prior_3_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)


prior_res_2$combed_plot
prior_res_2.25$combed_plot
prior_res_2.5$combed_plot
prior_res_2.75$combed_plot
prior_res_3$combed_plot


###############################################################################
## right lets run it with the poisson distribution instead ####################
###############################################################################

testing_priors_poisson <- function(stan_data, prior_start, prior_end,
                           length_o_priors){
  
  run_ids_3 <- NULL
  prior_vals <- seq(from = prior_start, to = prior_end,
                    length.out = length_o_priors)
  
  for(k in 1:5){
    j <- prior_vals[k]
    delta_priors <- c(0,(j*2),0.1,j)
    kappa_priors <- c(0, (j*2), 0, j)
    
    stan_data_discrete <- stan_data
    stan_data_discrete$delta_priors <- delta_priors
    stan_data_discrete$kappa_priors <- kappa_priors
    stan_data_discrete$poisson_or_negbin <- 0
    
    name_run <- paste("Prior_combo_Poisson:",as.character(j))
    
    params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen", "phi_pen", "delta_beta")
    
    cluster_run <- obj$enqueue(single_run_rw_delta(stan_data = stan_data_discrete,
                                                   parameters_to_mon = params_monitor_hiv,
                                                   warmup_length = 500,iter_length = 1500,
                                                   adapt_delta = 0.95), name = name_run)
    
    print(cluster_run$status())
    run_ids_3 <- c(run_ids_3, cluster_run$id)
    
    
    
  }
  
  return(run_ids_3)
}

run_1_5_outs <- testing_priors_poisson(stan_data_discrete, prior_start = 1,
                               prior_end = 5, length_o_priors = 5)

prior_poiss_1_combo <- obj$task_get(run_1_5_outs[1])
prior_poiss_2_combo <- obj$task_get(run_1_5_outs[2])
prior_poiss_3_combo <- obj$task_get(run_1_5_outs[3])
prior_poiss_4_combo <- obj$task_get(run_1_5_outs[4])
prior_poiss_5_combo <- obj$task_get(run_1_5_outs[5])

prior_poiss_1_combo$status()
prior_poiss_2_combo$status()
prior_poiss_3_combo$status()
prior_poiss_4_combo$status()
prior_poiss_5_combo$status()

prior_poiss_1_combo$log()
prior_poiss_2_combo$log()
prior_poiss_3_combo$log()
prior_poiss_4_combo$log()
prior_poiss_5_combo$log()

prior_poiss_1_combo_res <- prior_poiss_1_combo$result()
prior_poiss_2_combo_res <- prior_poiss_2_combo$result()
prior_poiss_3_combo_res <- prior_poiss_3_combo$result()
prior_poiss_4_combo_res <- prior_poiss_4_combo$result()
prior_poiss_5_combo_res <- prior_poiss_5_combo$result()

prior_poiss_res_1<-plot_stan_model_fit(model_output = prior_poiss_1_combo_res,
                                 plot_name = "Random walk second order",xout = xout,
                                 sim_output = sim_model_output, diags_dat = diags_df,
                                 X_Design = stan_data_discrete$X_design_diag)

prior_poiss_res_2<-plot_stan_model_fit(model_output = prior_poiss_2_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_3<-plot_stan_model_fit(model_output = prior_poiss_3_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_4<-plot_stan_model_fit(model_output = prior_poiss_4_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_5<-plot_stan_model_fit(model_output = prior_poiss_5_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)


prior_poiss_res_1$combed_plot
prior_poiss_res_2$combed_plot
prior_poiss_res_3$combed_plot
prior_poiss_res_4$combed_plot
prior_poiss_res_5$combed_plot

###############################################################################
## So the wide priors really aren't working for this, lets try them out #######
## between 0.5 and 1.5 ########################################################
###############################################################################

run_0.5_1.5_outs <- testing_priors_poisson(stan_data_discrete, prior_start = 0.5,
                                       prior_end = 1.5, length_o_priors = 5)

prior_poiss_0.5_combo <- obj$task_get(run_0.5_1.5_outs[1])
prior_poiss_0.75_combo <- obj$task_get(run_0.5_1.5_outs[2])
prior_poiss_1_combo <- obj$task_get(run_0.5_1.5_outs[3])
prior_poiss_1.25_combo <- obj$task_get(run_0.5_1.5_outs[4])
prior_poiss_1.5_combo <- obj$task_get(run_0.5_1.5_outs[5])


prior_poiss_0.5_combo$status()
prior_poiss_0.75_combo$status()
prior_poiss_1_combo$status()
prior_poiss_1.25_combo$status()
prior_poiss_1.5_combo$status()

prior_poiss_0.5_combo$log()
prior_poiss_0.75_combo$log()
prior_poiss_1_combo$log()
prior_poiss_1.25_combo$log()
prior_poiss_1.5_combo$log()

prior_poiss_0.5_combo_res <- prior_poiss_0.5_combo$result()
prior_poiss_0.75_combo_res <- prior_poiss_0.75_combo$result()
prior_poiss_1_combo_res <- prior_poiss_1_combo$result()
prior_poiss_1.25_combo_res <- prior_poiss_1.25_combo$result()
prior_poiss_1.5_combo_res <- prior_poiss_1.5_combo$result()

prior_poiss_res_0.5<-plot_stan_model_fit(model_output = prior_poiss_0.5_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)

prior_poiss_res_0.75<-plot_stan_model_fit(model_output = prior_poiss_0.75_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_1<-plot_stan_model_fit(model_output = prior_poiss_1_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_1.25<-plot_stan_model_fit(model_output = prior_poiss_1.25_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_1.5<-plot_stan_model_fit(model_output = prior_poiss_1.5_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)


prior_poiss_res_0.5$combed_plot
prior_poiss_res_0.75$combed_plot
prior_poiss_res_1$combed_plot
prior_poiss_res_1.25$combed_plot
prior_poiss_res_1.5$combed_plot

###############################################################################
## None of those runs were any good, lets try one with some really wide #######
## priors now #################################################################
###############################################################################

run_2_12_outs <- testing_priors_poisson(stan_data_discrete, prior_start = 2,
                                           prior_end = 12, length_o_priors = 5)

prior_poiss_2_combo <- obj$task_get(run_2_12_outs[1])
prior_poiss_4.5_combo <- obj$task_get(run_2_12_outs[2])
prior_poiss_7_combo <- obj$task_get(run_2_12_outs[3])
prior_poiss_9.5_combo <- obj$task_get(run_2_12_outs[4])
prior_poiss_12_combo <- obj$task_get(run_2_12_outs[5])


prior_poiss_2_combo$status()
prior_poiss_4.5_combo$status()
prior_poiss_7_combo$status()
prior_poiss_9.5_combo$status()
prior_poiss_12_combo$status()

prior_poiss_2_combo$log()
prior_poiss_4.5_combo$log()
prior_poiss_7_combo$log()
prior_poiss_9.5_combo$log()
prior_poiss_12_combo$log()

prior_poiss_2_combo_res <- prior_poiss_2_combo$result()
prior_poiss_4.5_combo_res <- prior_poiss_4.5_combo$result()
prior_poiss_7_combo_res <- prior_poiss_7_combo$result()
prior_poiss_9.5_combo_res <- prior_poiss_9.5_combo$result()
prior_poiss_12_combo_res <- prior_poiss_12_combo$result()

prior_poiss_res_2<-plot_stan_model_fit(model_output = prior_poiss_2_combo_res,
                                         plot_name = "Random walk second order",xout = xout,
                                         sim_output = sim_model_output, diags_dat = diags_df,
                                         X_Design = stan_data_discrete$X_design_diag)

prior_poiss_res_4.5<-plot_stan_model_fit(model_output = prior_poiss_4.5_combo_res,
                                          plot_name = "Random walk second order",xout = xout,
                                          sim_output = sim_model_output, diags_dat = diags_df,
                                          X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_7<-plot_stan_model_fit(model_output = prior_poiss_7_combo_res,
                                       plot_name = "Random walk second order",xout = xout,
                                       sim_output = sim_model_output, diags_dat = diags_df,
                                       X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_9.5<-plot_stan_model_fit(model_output = prior_poiss_9.5_combo_res,
                                          plot_name = "Random walk second order",xout = xout,
                                          sim_output = sim_model_output, diags_dat = diags_df,
                                          X_Design = stan_data_discrete$X_design_diag)
prior_poiss_res_12<-plot_stan_model_fit(model_output = prior_poiss_12_combo_res,
                                         plot_name = "Random walk second order",xout = xout,
                                         sim_output = sim_model_output, diags_dat = diags_df,
                                         X_Design = stan_data_discrete$X_design_diag)


prior_poiss_res_2$combed_plot
prior_poiss_res_4.5$combed_plot
prior_poiss_res_7$combed_plot
prior_poiss_res_9.5$combed_plot
prior_poiss_res_12$combed_plot



