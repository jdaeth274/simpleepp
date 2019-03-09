###############################################################################
## Running the simple SID HIV model from STAN #################################
###############################################################################
require(reshape2)
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
  sigma <- params$sigma           # Progression from stages of infection
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

data_collector <- function(sim_res, year_range, xstart){
  
  start_val <- (year_range[1] - xstart)*10 + 1
  end_val <- ((year_range[1] - xstart)*10) + ((length(year_range) - 1) * 10) + 1
  
  diagnoses_tot <- sim_res$diagnoses
  
  diag_to_keep <- seq(from  = start_val, to = end_val, by = 10)
  
  diag_to_keep <- round(diagnoses_tot[diag_to_keep])
  
  return(diag_to_keep)
}

diags_data <- data_collector(sim_model_output$sim_df, year_range, xstart)
###############################################################################
## Now we've got our data from the simulated model, we can also generate the ##
## RW matrices for the input into the dataframe. ##############################
###############################################################################

xout<-seq(1970,2019.9,0.1)
spline_matrix<-splineDesign(1969:2021,xout,ord = 2)            ## This matrix is the spline design one 
penalty_matrix<-diff(diag(ncol(spline_matrix)), diff=2)        ## This matrix creates the differences between your kappa values 
rows_to_evaluate<-0:45*10+1

stan_data_discrete<-list(
  n_obs = length(year_range),
  y = as.array(diags_data),
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

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen")  

test_stan_hiv<- stan("~/Dropbox/jeff_hiv_work/simpleepp/stan_files/chunks/SID_models/SID_model.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  


mod_hiv_prev <- stan("hiv_project/simpleepp/stan_files/chunks/cd4_matrix_random_walk.stan", data = stan_data_discrete,
                     pars = params_monitor_hiv,chains = 3,warmup = 500,iter = 1500,
                     control = list(adapt_delta = 0.99))

  
  
  


