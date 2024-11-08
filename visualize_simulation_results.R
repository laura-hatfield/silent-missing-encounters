library(tidyverse)
theme_set(theme_minimal()+theme(legend.position = "bottom"))

#### Load and process full-scale simulation results ####
raw <- readRDS("silent_missing_sim_results.rds")

## Setting up for nice labels on the parameters that varied across simulation reps
p_ed_param1b.labs <- c("Strong plan effect on\nmissingness", "No plan effect on\nmissingness")
names(p_ed_param1b.labs) <- c("-2", "0")

ed2_beta2.labs <- c("Strong", "None")
names(ed2_beta2.labs) <- c('-1.5', '0')

low_concern_beta1.labs <- c('Strong plan effect on\nhigh-quality indicator','No plan effect on\nhigh-quality indicator')
names(low_concern_beta1.labs) <- c('4','0')

low_concern_beta0.labs <- c('20% low-concern','50% low-concern','80% low-concern')
names(low_concern_beta0.labs) <- c('-1.39','0','1.39')

p_ed_param1a.labs <- c("High missingness", "Medium missingness", "Low missingness")
names(p_ed_param1a.labs) <- c("-1", "-1.5", "-3")

low_concern.labs <- c("High concern", "Low concern")
names(low_concern.labs) <- c("0", "1")

#### Summarize results across simulation replications ####

results_pretty <- raw %>%
  # Compute ED visit rate in each plan
  mutate(true_rate=sum_ed/nj,
         appr1_est_rate=sum_appr1_ed_hat/nj,
         appr2_est_rate=sum_appr2_ed_hat/nj,
         appr3_est_rate=ifelse(low_concern==1,sum_taf_ed/nj,NA),
         appr4_est_rate=sum_taf_ed/nj) %>%
  # In each simulation run for each parameter combination, compute mean ED visit rate (across plans)
  group_by(runID, p_ed_param1b, ed2_beta2, ed2_nz_beta2,
           low_concern_beta1, illness_param2, p_ed_param1a, low_concern_beta0) %>% 
  summarize(across(c(true_rate,appr1_est_rate,appr2_est_rate,appr3_est_rate,appr4_est_rate),~mean(.x,na.rm=TRUE)),.groups="keep") %>%
  # Compute error compared to true ED visit rate
  mutate(across(c(appr1_est_rate,appr2_est_rate,appr3_est_rate,appr4_est_rate),~(.x-true_rate)/true_rate*100,.names="err_{.col}")) %>%
  # Code scenarios by realized missingness percent (which is the TAF as-is error rate)
  mutate(realized_miss=cut(abs(err_appr4_est_rate),breaks=c(0,10,20,100),
                           labels=c('Missingness <10%','Missingness 10-20%','Missingness >20%')),
         # Give the parameters nicer names
         low_concern_beta1=factor(low_concern_beta1,labels=low_concern_beta1.labs,levels=names(low_concern_beta1.labs)),
         ed2_beta2=factor(ed2_beta2,labels=ed2_beta2.labs,levels=names(ed2_beta2.labs)),
         p_ed_param1a=factor(p_ed_param1a,labels=p_ed_param1a.labs,levels=names(p_ed_param1a.labs)),
         p_ed_param1b=factor(p_ed_param1b,labels=p_ed_param1b.labs,levels=names(p_ed_param1b.labs)),
         low_concern_beta0=factor(low_concern_beta0,labels=low_concern_beta0.labs,levels=names(low_concern_beta0.labs)))

results_long <- results_pretty %>% select(runID, p_ed_param1b, ed2_beta2, ed2_nz_beta2,
                                          low_concern_beta1, illness_param2, p_ed_param1a, low_concern_beta0,
                                          realized_miss,err_appr1_est_rate,err_appr2_est_rate,err_appr3_est_rate,err_appr4_est_rate) %>%
  pivot_longer(cols=c(err_appr1_est_rate,err_appr2_est_rate,err_appr3_est_rate,err_appr4_est_rate),values_to="err_pct") %>%
  mutate(approach=factor(name,levels=paste0("err_appr",1:4,"_est_rate"),labels=c('Learn Utilization','Learn Missingness','High-quality only','Uncorrected')))

# Final Figure: 
ggplot(results_long,aes(x=ed2_beta2,y=err_pct)) + geom_hline(yintercept=0) + 
  geom_boxplot(aes(col=approach),outlier.size=0.5) + 
  facet_grid(low_concern_beta1~realized_miss) +
  scale_color_manual("",labels=c("Learn Utilization","Learn Missingness","High-quality only","Uncorrected"),values=c("#E69F00", "#00BFC4", "#7CAE00",  "#C77CFF")) +
  scale_y_continuous("% difference from true visit rate") + xlab("Plan effect on ED visits")
