## Load the simulation functions
source("silent_missing_simulation_fns.R")

# Set number of simualtion replications
nsims = 1000

#### Run simulation using parallel processing  ####

# Create combinations of parameters that we vary
parameters = as_tibble(expand_grid(p_ed_param1b = c(0, -2) , # missingness <-> omega
                                   ed2_beta2 = c(0, -1.5), # ED2 <-> omega
                                   ed2_nz_beta2 = c(0, -0.25), # ED2 <-> omega
                                   low_concern_beta1 = c(0, 4), # low-concern <-> omega 
                                   illness_param2 = c(0, -1), # health mean <-> omega
                                   p_ed_param1a = c(-1,-1.5,-3), # overall missingness
                                   low_concern_beta0 = c(-1.39,0,1.39)) # overall prop. low concern
                       ) %>% 
  # Params that don't vary
  mutate(plan_n_mean = 100, #
         alpha = 0.25,
         j = 200, # 
         illness_param1 = -1,
         illness_var = .25,
         omega_param1 = 3,
         omega_param2 = 1, 
         omega_hat_sd = .1,
         ed1_nz_beta0 = -1,
         ed1_nz_beta1 = .5,
         ed1_beta0 = 0.5,
         ed1_beta1 = 1, 
         ed2_nz_beta0 = -1,
         ed2_nz_beta1 = 1,
         ed2_beta0 = 0,
         ed2_beta1 = 2,
         p_ed_var = 0.25) %>%
  # Get rid of parameter combinations with "mixed" values of ed_beta2 and ed2_nz_beta2
  filter(((ed2_beta2 < 0)&(ed2_nz_beta2 < 0))|((ed2_beta2 == 0)&(ed2_nz_beta2 ==0)))

## Set up the parallelization infrastructure 
cl <- parallel::makeCluster(10) 
doParallel::registerDoParallel(cl) 

## Run the simulation
results <- foreach(i = 1:nsims,
                   .packages=c('dplyr','tidyr','purrr'),
                   .combine = rbind) %dopar% { 
                     parameters %>% 
                       mutate(results = pmap( parameters, run_simulation ) ) %>%
                       unnest(results) %>% 
                       mutate(runID = i)
                   } # end of dopar

stopCluster(cl) #end the parallelization 

saveRDS(results, file = paste0(Sys.Date(),"_silent_missing_sim_results.rds"))
