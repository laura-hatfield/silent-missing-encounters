### DGP with parallel processing ###

#### 0. Load libraries ####
library(tidyverse)
#library(dplyr)
#library(tidyr)
#library(boot)
#library(tibble)
#library(purrr)
#library(splitstackshape)
#library(data.table)
library(doParallel)
library(foreach)


#### 1. Function to generate data  ####

generate_data <- function(plan_n_mean, alpha, j,
                          illness_param1, illness_param2, illness_var,  
                          omega_param1, omega_param2, omega_hat_sd,
                          low_concern_beta0, low_concern_beta1,
                          ed1_nz_beta0, ed1_nz_beta1, 
                          ed1_beta0, ed1_beta1, 
                          ed2_nz_beta0, ed2_nz_beta1, ed2_nz_beta2,
                          ed2_beta0, ed2_beta1, ed2_beta2,
                          p_ed_param1a, p_ed_param1b,p_ed_var) {
  
  # Need to add more of these checks:
  if (illness_var >= 1 | illness_var <= 0) stop("illness_var must be in (0,1)")
  if (p_ed_var >= 1 | p_ed_var <= 0) stop("p_ed_var must be in (0,1)")
  
  # Generate plan sizes
  n_min = round(plan_n_mean * (1 - alpha) )
  n_max = round(plan_n_mean * (1 + alpha) )

  # Generate plan-level data
  simdata_plans <- tibble(planID = 1:j,
                          # uniform distribution on plan sizes
                          nj=sample (n_min:n_max, j, replace=TRUE),
                          # plan level mean illness
                          omega=rbeta(j, omega_param1, omega_param2),
                          # noise added to plan-level mean illness
                          omega.err=rnorm(j,0,omega_hat_sd)) %>% 
    mutate(noisy_omega=omega+omega.err,  
           # Mean of the illness parameters across plans
           center_omega=omega-mean(omega), 
           center_noisy_omega=noisy_omega-mean(noisy_omega),
           # Low-concern indicator depends on centered omega
           low_concern=rbinom(j,1,inv.logit(low_concern_beta0+low_concern_beta1*center_omega)), 
           # Mean of illness distribution depends on centered omega
           theta_mean = inv.logit(illness_param1 + illness_param2*center_omega),
           # Mean of missingness proportion depends on centered omega
           p_ed_mean = ifelse(low_concern==1,0,inv.logit(p_ed_param1a + p_ed_param1b*center_omega))) %>%
    select(-omega.err)

  total_n = sum(simdata_plans$nj)

  # Generate person-level identifiers, merge to plan, and generate 
  # * individual illness level theta
  # * true emergency ED visits ed1 (as a "two-part model")
  # * non-emergency ED visits ed2 (as a "two-part model")
  
  simdata_taf <- tibble(personID=as.factor(c(1:total_n)),
                          planID = rep(1:j, simdata_plans$nj)) %>%
    left_join(simdata_plans,by="planID") %>%
    # Distribution of health depends on plan characteristics
    # (higher quality plans have healthier enrollees;  note that theta is illness level, theta=0 perfect health)
    mutate(theta = rbeta(total_n, theta_mean*(1-illness_var)/illness_var, (1-theta_mean)*(1-illness_var)/illness_var), 
           # Injury emergency visits depend on health only
           ed1_nz_prop = inv.logit(ed1_nz_beta0 + ed1_nz_beta1*(theta-theta_mean)), 
           ed1_lambda = exp(ed1_beta0 + ed1_beta1*(theta-theta_mean)),
           ed1_0 = rbinom(total_n,1,ed1_nz_prop),
           ed1 = ifelse( ed1_0 == 0, 0, rpois(total_n, lambda = ed1_lambda)),
           # Other (non-)emergency visits depend on health and plan quality
           ed2_nz_prop = inv.logit(ed2_nz_beta0 + ed2_nz_beta1*(theta-theta_mean) + ed2_nz_beta2*center_omega),
           ed2_lambda= exp(ed2_beta0 + ed2_beta1*(theta-theta_mean) + ed2_beta2*center_omega),
           ed2_0 = rbinom(total_n,1,ed2_nz_prop),
           ed2 = ifelse( ed2_0 == 0, 0, rpois(total_n, lambda = ed2_lambda)),
           ed = ed1 + ed2,
           # Used later in the approaches
           health_adj1=ed1_lambda*ed1_nz_prop,
           health_adj2=ed1_lambda*ed2_nz_prop,
           taf_ed1=ifelse(ed1==0,0,rbinom(total_n,ed1,1-p_ed_mean)), 
           taf_ed2=ifelse(ed2==0,0,rbinom(total_n,ed2,1-p_ed_mean)),
           taf_ed=taf_ed1+taf_ed2)

  return(simdata_taf)
}
# end of function to generate the data

#### 2. Function for Approach 1  ####

approach1 <- function(data) {
  
  # Step 1. Estimate relationship between omega and utilization, 
  # using only low_concern states (where true = observed by assumption)
  # and using only ED2 (where there is a relationship)
  
  # 1.3. Regress ED2 on omega hat in low-concern plans
  n_low_concern <- dim(filter(data,low_concern==1) %>% 
                         group_by(planID) %>% summarize())[1]
  
  if (n_low_concern > 1){
    # 1.1. Casemix adjust each type of ED visit
    simdata_apr1 <- data %>% 
      mutate(taf_ed2_adj=taf_ed2-health_adj2, # This is "excess" utilization beyond health
             true_ed2_adj=ed2-health_adj2)

    # In low-concern plans, regress residual visits on plan proxy for utilization
    a1_ed_omega_reg <- lm(formula = taf_ed2_adj ~ center_noisy_omega,
                          data = subset(simdata_apr1, low_concern == 1))
    # Predict the quality portion of health-adjusted ED utilization for *all* plans
    data_with_predict <- simdata_apr1 %>% 
      bind_cols(ed2_qual = predict(a1_ed_omega_reg,newdata=simdata_apr1,type="response")) %>%
      # Add together the quality and health predictions for ed
      # mutate(ed_hat = ifelse(low_concern == 1, taf_ed, 
      #                         ifelse((ed2_qual + health_adj2 + health_adj1)<1,0,
      #                               ed2_qual + health_adj2 + health_adj1)),
      #        ed_hat = ifelse(ed_hat > taf_ed, ed_hat, taf_ed))
      mutate(ed_hat_1=ed2_qual + health_adj2 + health_adj1, # calculate our estimated ED visits
             ed_hat_2=ifelse(ed_hat_1<1,0,ed_hat_1), # quick and dirty 2-part model fix
             # sum up to plan here and calculate 
             ###
             ed_hat_3=ifelse(low_concern == 1, taf_ed, ed_hat_2), # replace with TAF in low concern plans (no missingness)
             ed_hat=ifelse(ed_hat_3 > taf_ed, ed_hat_3, taf_ed)) # replace with TAF if estimated is less than TAF (implement connection to TAF) 
    
    # Create tibble for the results data
    return(data_with_predict %>% select(personID,planID,omega,noisy_omega,center_noisy_omega,low_concern,nj,p_ed_mean,
                                        ed1,ed2,ed,health_adj1,health_adj2,
                                        taf_ed1,taf_ed2,taf_ed,
                                        ed2_qual,ed_hat_1,ed_hat_2,ed_hat_3,ed_hat))
    
  } else {
    print(paste0("there are only 0 or 1 low concern plans"))
    return(simdata_apr1 %>%
             mutate(ed_hat=NA) %>%
             select(personID,planID,omega,noisy_omega,center_noisy_omega,low_concern,nj,p_ed_mean,
                    ed1,ed2,ed,health_adj1,health_adj2,
                    taf_ed1,taf_ed2,taf_ed,
                    ed2_qual,ed_hat_1,ed_hat_2,ed_hat_3,ed_hat))
  }
}
# end of function for Approach 1

#### 3. Function for Approach 2  ####

approach2 <- function( data ) {
  # Health-adjusted utilization
  simdata_apr2 <- data %>% 
    mutate(taf_ed1_adj=taf_ed1-health_adj1,
           taf_ed2_adj=taf_ed2-health_adj2)
  
  # compute plan-level ratios of non-emergency to emergency visits
  plan_ED_mean <- simdata_apr2 %>% group_by(planID) %>% 
    summarise(taf_ed1_planmean=mean(taf_ed1), 
              taf_ed2_planmean=mean(taf_ed2),
              taf_ed2toed1_multiplier=taf_ed2_planmean/taf_ed1_planmean)
  
  simdata_apr2 <- simdata_apr2 %>% left_join(plan_ED_mean,by="planID")
  
  # Logic: after health adjustment, 
  #   any variation in injury ED visits is due to data quality
  # Predict expected injury ED visits (by plan quality)
  a2_ed_omega_reg <- lm(formula = taf_ed1_adj ~ center_noisy_omega,data = simdata_apr2)
  data_with_predict <- simdata_apr2 %>%
    bind_cols(a2_ed1_delta = predict(a2_ed_omega_reg,type="response")) %>%
    # Apply plan-level adjustments to non-injury ED visits 
    # need to adjust this delta subtraction to ED1 and ED2 distribution
    mutate(a2_ed2_hat = taf_ed2 - a2_ed1_delta*taf_ed2toed1_multiplier,
           a2_ed1_hat = taf_ed1 - a2_ed1_delta,
           ed_hat_1 = ifelse(a2_ed1_hat + a2_ed2_hat > 0,a2_ed1_hat + a2_ed2_hat,0), # calculate our estimated ED visits
           ed_hat_2 = ifelse(ed_hat_1<1,0,ed_hat_1), # quick and dirty 2-part model fix
           ed_hat = ifelse(ed_hat_2 > taf_ed, ed_hat_2, taf_ed)) # replace with TAF if estimated is less than TAF (implement connection to TAF)
  
  
  # Create tibble for the results data
  tibble(data_with_predict %>% select(personID,planID,omega,noisy_omega,center_noisy_omega,low_concern,nj,p_ed_mean,
                                      ed1,ed2,ed,health_adj1,health_adj2,
                                      taf_ed1,taf_ed2,taf_ed,a2_ed1_delta,
                                      a2_ed2_hat,a2_ed1_hat,ed_hat_1,ed_hat_2,ed_hat,taf_ed2toed1_multiplier))
}
# end of function for Approach 2

#### 4. Function to analyze data with the two approaches  ####

analyze_data <- function( data) {
  appr1 <- approach1( data )
  appr2 <- approach2( data )
  long_dat <- bind_rows(appr1 = appr1, appr2 = appr2, .id="approach") %>%
    # Create 2 individual-level metrics: 
    #   1. difference of estimated and true, 
    #   2. difference of estimated and taf
    mutate(est_ed_diff = ed_hat - ed,
           est_taf_diff = ed_hat - taf_ed) %>%
    group_by(approach,planID) %>%
    # create plan-level metrics (mean differences and total sum differences)
    summarise(mean_est_ed_diff = mean( est_ed_diff ),
              mean_est_taf_diff = mean( est_taf_diff ), 
              sum_est_ed_diff = sum( est_ed_diff),
              sum_est_taf_diff = sum( est_taf_diff),
              omega = mean(omega),
              nj=mean(nj), 
              noisy_omega = mean(noisy_omega),
              center_noisy_omega = mean(center_noisy_omega),
              mean_proportion_missing = mean(p_ed_mean), 
              sum_ed = sum(ed),
              sum_taf_ed = sum(taf_ed),
              sum_ed_hat = sum(ed_hat),
              mean_ed = mean(ed),
              mean_taf_ed = mean(taf_ed),
              mean_ed_hat = mean(ed_hat),
              low_concern = mean(low_concern),
              mean_health_adj1 = mean(health_adj1),
              mean_health_adj2 = mean(health_adj2),.groups="keep")
  return(long_dat)
}
# end of function to analyze data

#### 5. Function to run the simulation ####

run_simulation <- function(plan_n_mean, alpha, j,
                           illness_param1, illness_param2, illness_var,
                           omega_param1, omega_param2, omega_hat_sd,
                           low_concern_beta0,low_concern_beta1,
                           ed1_nz_beta0, ed1_nz_beta1, 
                           ed1_beta0, ed1_beta1, 
                           ed2_nz_beta0, ed2_nz_beta1, ed2_nz_beta2,
                           ed2_beta0, ed2_beta1, ed2_beta2,
                           p_ed_param1a, p_ed_param1b,p_ed_var) {
  # Generate the data
  mydata <- generate_data(plan_n_mean, alpha, j,
                          illness_param1, illness_param2, illness_var,
                          omega_param1, omega_param2, omega_hat_sd,
                          low_concern_beta0,low_concern_beta1,
                          ed1_nz_beta0, ed1_nz_beta1, 
                          ed1_beta0, ed1_beta1, 
                          ed2_nz_beta0, ed2_nz_beta1, ed2_nz_beta2,
                          ed2_beta0, ed2_beta1, ed2_beta2,
                          p_ed_param1a, p_ed_param1b,p_ed_var)
  # Analyze the data
  results <- analyze_data(mydata)

  return(results)
}
# end of function to clean up and group the data

