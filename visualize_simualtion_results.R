
raw <- readRDS("silent_missing_sim_results.rds")


## Give the parameters better names
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


## Instead of showing plan-level distributions across simulation reps,
## suppose we were trying to estimate the per-person # ED visits
## Compare the distribution (across simulation reps) of the estimated ED visits per person
## in our 2 approaches plus just using low-concern states
appr12 <- raw %>%
  # within plans, compute ED visit rates
  mutate(true_rate=sum_ed/nj,
         est_rate=sum_ed_hat/nj) %>%
  # across plans, in each approach and simulation run, estimate the mean per person ED visit rate:
  group_by(runID, approach, p_ed_param1b, ed2_beta2, ed2_nz_beta2, 
                                low_concern_beta1, illness_param2, p_ed_param1a, 
                                low_concern_beta0) %>% 
  summarize(true_rate = mean(true_rate),
            est_rate = mean(est_rate)) %>%
  mutate(error_pct=(est_rate-true_rate)/true_rate*100)

## Use low-concern states' observed (TAF) data to estimate, but compare to truth for all plans
appr3 <- raw %>% filter(approach=="appr1") %>% # approach doesn't matter; just picking one so as not to double count
  mutate(true_rate=sum_ed/nj, 
         est_rate=ifelse(low_concern==1,sum_taf_ed/nj,NA),
         approach="appr3") %>%
  # across plans, in each approach and simulation run, estimate the mean per person ED visit rate:
  group_by(runID, approach, p_ed_param1b, ed2_beta2, ed2_nz_beta2, 
                                low_concern_beta1, illness_param2, p_ed_param1a, 
                                low_concern_beta0) %>% 
  summarize(true_rate = mean(true_rate),
            # na.rm=T ensures that the non-low-concern states' estimates are not included
            est_rate = mean(est_rate,na.rm=T)) %>%
  mutate(error_pct=(est_rate-true_rate)/true_rate*100)

## Use observed TAF data in all states
appr4 <- raw %>% filter(approach=="appr1") %>% # approach doesn't matter; just picking one so as not to double count
  mutate(true_rate=sum_ed/nj,
         est_rate=sum_taf_ed/nj,
         approach="appr4") %>%
  # across plans, in each approach and simulation run, estimate the mean per person ED visit rate:
  group_by(runID, approach, p_ed_param1b, ed2_beta2, ed2_nz_beta2, 
                                low_concern_beta1, illness_param2, p_ed_param1a, 
                                low_concern_beta0) %>% 
  summarize(true_rate = mean(true_rate),
            est_rate = mean(est_rate)) %>%
  mutate(error_pct=(est_rate-true_rate)/true_rate*100)

baseline <- appr4 %>% mutate(baseline_rate=est_rate) %>% select(-c("est_rate","approach","error_pct"))

true <- baseline %>% mutate(approach='true')

rate_compare <- bind_rows(appr12,appr3) %>% bind_rows(appr4)

rate_compare_2 <- left_join(rate_compare, baseline, by=c('runID','p_ed_param1b', 'ed2_beta2', 'ed2_nz_beta2', 
                                'low_concern_beta1', 'illness_param2', 'p_ed_param1a', 
                                'low_concern_beta0')) %>%
  select(-c("approach.y","true_rate.y")) %>% rename(approach=approach.x,true_rate=true_rate.x)

x <- rate_compare %>% filter(approach=='appr4') %>% mutate(error_pct_abs=abs(error_pct))
x$realized_miss <- cut(x$error_pct_abs,
                       breaks=c(0,10,20,30),
                       labels=c('Missingness <10%','Missingness 10-20%','Missingness >20%'))
summary(x$error_pct)


realized <- x %>% subset(select=-c(approach,true_rate,est_rate,error_pct,error_pct_abs))

rate_compare_3 <- bind_rows(rate_compare_2,true)

realized_allapproaches <- rate_compare_3 %>% left_join(realized,by=c("runID",'p_ed_param1b', 'ed2_beta2', 'ed2_nz_beta2', 
                                'low_concern_beta1', 'illness_param2', 'p_ed_param1a', 
                                'low_concern_beta0'))

realized_allapproaches <- realized_allapproaches %>%
mutate(
        low_concern_beta1=factor(low_concern_beta1,labels=low_concern_beta1.labs,levels=names(low_concern_beta1.labs)),
         ed2_beta2=factor(ed2_beta2,labels=ed2_beta2.labs,levels=names(ed2_beta2.labs)),
         p_ed_param1a=factor(p_ed_param1a,labels=p_ed_param1a.labs,levels=names(p_ed_param1a.labs)),
         p_ed_param1b=factor(p_ed_param1b,labels=p_ed_param1b.labs,levels=names(p_ed_param1b.labs)),
         low_concern_beta0=factor(low_concern_beta0,labels=low_concern_beta0.labs,levels=names(low_concern_beta0.labs))) %>% mutate(rate_to_TAF=est_rate/baseline_rate) 

true_2 <- realized_allapproaches %>% filter(approach=="true") %>% mutate(rate_to_TAF=true_rate/baseline_rate)

fourapproaches <- realized_allapproaches %>% filter(approach!="true") 

final_data <- bind_rows(fourapproaches,true_2)

# Final Figure: 
ggplot(filter(final_data,approach!="true"),
       aes(x=ed2_beta2,y=error_pct)) + geom_hline(yintercept=0) + 
  geom_boxplot(aes(col=approach),outlier.size=0.5) + 
  facet_grid(low_concern_beta1~realized_miss) +
  scale_color_manual("",labels=c("Learn Utilization","Learn Missingness","High-quality only","Uncorrected"),values=c("#E69F00", "#00BFC4", "#7CAE00",  "#C77CFF")) +
  scale_y_continuous("% difference from true visit rate") + xlab("Plan effect on ED visits")
