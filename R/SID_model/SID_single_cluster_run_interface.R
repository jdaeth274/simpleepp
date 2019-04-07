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

setwd("~/homes_drive/")
options(didehpc.username = "jd2117",didehpc.home = "~/homes_drive/simpleepp",didehpc.cluster = "fi--didemrchnb")

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

expose_stan_functions("~/Dropbox/jeff_hiv_work/simpleepp/stan_files/chunks/SID_models/SID_model_RW_delta.stan")

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
  names(sim_mod)<-c("kappa","lambda","prevalence","diagnoses")
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
  
  melt_df<-sim_df[,-c(5)]
  melted_df_funky<-melt(melt_df,id="time")
  
  four_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot, incidence=incidence_plot, kappa=kappa_plot, diagnoses = diagnoses_plot,
              whole=four_variable_plot))
  
  
}

plotted_sim<-sim_plot(sim_model_output$sim_df)
plot(plotted_sim$whole)
plot(plotted_sim$prevalence)


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


###############################################################################
## Now we've got our data from the simulated model, we can also generate the ##
## RW matrices for the input into the dataframe. ##############################
###############################################################################



xord<-seq(1970.1,2020,0.1)
xord_diag <- seq(1970,2020,0.1)
spline_matrix<-splineDesign(1969:2021,xord,ord = 2)            ## This matrix is the spline design one 
spline_mat_diag <- splineDesign(1969:2021,xord_diag,ord = 2)
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
  multiplicative_vals = c(1.02,1.04,1.1)
)

params_monitor_hiv<-c("y_hat","iota","fitted_output","beta","sigma_pen", "phi_pen", "delta_beta")  

test_stan_hiv<- stan("~/Dropbox/jeff_hiv_work/simpleepp/stan_files/chunks/SID_models/SID_model_RW_delta.stan",
                     data = stan_data_discrete,
                     pars = params_monitor_hiv,
                     chains = 1, iter = 10)  

mod_hiv_prev_pen_2 <- obj$enqueue(single_run_rw_delta(stan_data = stan_data_discrete,
                           parameters_to_mon = params_monitor_hiv,
                           warmup_length = 500,iter_length = 1500,
                           adapt_delta = 0.95), name = "single_run_rw")
mod_hiv_prev_pen_2$status()
mod_hiv_prev_pen_2$log()
