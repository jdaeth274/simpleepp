###############################################################################
## Setting up the SID model for the DIDE cluster ##############################
###############################################################################

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
                            sources = "simpleepp/R/SID_model/SID_model_cluster_function.R")
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)

################################################################################################################################
## So the above lines of code give me access to the cluster, with the obj object giving me queue functionalaity, I've also #####
## asked for three cores for stan to use for each of the chains ################################################################
################################################################################################################################
obj$cluster_load()
obj$task_status()


diag_samples <- read.table("~/Dropbox/jeff_hiv_work/simpleepp/analysis/SID_models_diag_10_3_2019.tsv")
diag_samples <- diag_samples[-1,]


mu <- 1/35                               # Non HIV mortality / exit from population
sigma <- 1/c(3.16, 2.13, 3.20, 3.1)           # Progression from stages of infection
mu_i <- c(0.003, 0.008,0.01, 0.035, 0.27)     # Mortality by stage, no ART
kappa<-c(0.5,0.1,0.3,1995)
delta <- c(0.001,0.002,0.005,0.009)
iota<-0.0001

params <- NULL

dt <- 0.1                                     # time step
nsteps <- as.integer(50/dt)                   # number of steps
xstart <- 1970                                # the start of the epidemic
step_vector <- seq(xstart+dt, by=dt, length.out=nsteps)  # steps
xout <- c(xstart, step_vector)                #the vector of steps in total

###################################################################################################################################
## Now we will run through the data from begininng ################################################################################
###################################################################################################################################

sample_range<-1974:2015
                            ##### !!!!!!!!!!!!!!!!!!!!!! Remember to change this for when you sample
penalty_order<-1
rows_to_evaluate<-4:45*10+1
diag_samples <- diag_samples[,-c(1:4)]

data_about_sampling<-list(penalty_order=penalty_order,sample_years=length(sample_range),
                          rows_to_evaluate=rows_to_evaluate, mu = mu,
                          mu_i = mu_i, delta = delta, iota = iota, dt = dt, check_rhat = FALSE)

diag_SID_model_test <- obj$enqueue(fitting_SID_model(samples_data_frame = diag_samples, 
                                                                  data_about_sampling = data_about_sampling,
                                                                  iteration_number = 100,params = params,
                                                                  simulated_true_df = sim_model_output$sim_df), 
                                   name = "SID model Neg Binom")
                       


diag_SID_model_test$status()
diag_SID_model_test$context_id()

diag_sid_mod_id <- diag_SID_model_test$id
diag_SID_model_test$id

n_100_RW_first_order_loop_id<-n_100_RW_first_order_loop$id
save(diag_sid_mod_id,
     file = "~/Dropbox/jeff_hiv_work/simpleepp/data/SID_mod_cluster_id_23_3_19")
load("~/Dropbox/jeff_hiv_work/simpleepp/data/SID_mod_cluster_id_23_3_19", verbose = T)
diag_SID_model_test <- obj$task_get(diag_sid_mod_id)
