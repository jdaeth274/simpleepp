###############################################################################
## Creating a more stochastic diagnosis rate curve ############################
###############################################################################

require(splines)
require(ggplot2)
require(rstan)

###############################################################################
## Lets form the RW matrix first for our runs and they try and add in some ####
## stochasticity ##############################################################
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

expose_stan_functions("~/Dropbox/jeff_hiv_work/simpleepp/stan_files/chunks/SID_models/SID_model_RW_delta.stan")

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
  
  
  melt_df<-sim_df[,-c(4,5)]
  melted_df_funky<-melt(melt_df,id="time")
  
  four_variable_plot<-ggplot(data = melted_df_funky)+geom_line(aes(x=time,y=value,colour=variable),size=1.05)+
    labs(x="Time",y="Value",title="Simulated model output")
  
  return(list(prevalence=prev_plot, incidence=incidence_plot, kappa=kappa_plot, diagnoses = diagnoses_plot,
              whole=four_variable_plot))
  
  
}

plotted_sim<-sim_plot(sim_model_output$sim_df)
plot(plotted_sim$whole)
plot(plotted_sim$diagnoses)

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



