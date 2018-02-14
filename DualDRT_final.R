# load libraries
library(readr)
library(dplyr)
library(forcats)
library(ggplot2)

# Read in relevant datsets --------------
visual <- read_csv("visual.csv") %>% 
  group_by(sub) %>% 
  mutate(n_cars = n_distinct(model)) %>%
  filter(n_cars > 3)

cognitive <- read_csv("cognitive.csv") %>% 
  group_by(sub) %>% 
  mutate(n_cars = n_distinct(model)) %>%
  filter(n_cars > 3)



# Demographics --------------
demographics <- read_csv("5AAA_MASTER_Demographics.csv") %>%
  mutate(sub = parse_integer(sub)) %>%
  filter(Validity == 1) %>%
  group_by(sub, gender) %>%
  summarise(age = mean(age))

# makes a dataframe showing how many times each participant
# has repeated
n_p <- visual %>%
  group_by(sub) %>%
  summarise(repeats = mean(n_cars)) 
summary(n_p)





# Linear Mixed Effects Models--------------
library(nlme)
library(multcomp)

# Visual Load

# makes a dataframe with the mean hit rate to each stimulus
# for each condition (Single, SuRT, N-back) at the participant level
d <- visual %>%
  group_by(sub, task, stim) %>%
  summarise(hr_m = mean(hit_rate)) %>%
  filter(!is.na(hr_m)) %>%
  ungroup() %>%
  mutate(task = as.factor(task),
         stim = as.factor(stim),
         sub = as.factor(sub))

# linear mixed effects model 
model1 <- lme(hr_m ~ task * stim, random = ~1|sub/task/stim, data = d)
anova(model1)
summary(glht(model1, linfct = mcp(task = "Tukey")), test = adjusted(type="bonferroni"))

# boxplot of Hit Rate 
d %>%
  mutate(Task = fct_recode(task, `N-back` = "nBack"),
         Stimulus = fct_recode(stim, `Remote LED` = "LED_dash",
                               `Vibrotactor` = "tactor")) %>%
  ggplot(mapping = aes(x = task, y = hr_m, fill = Stimulus)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#f0f0f0", "#999999")) +
  theme_bw() +
  ylab("Mean Hit Rate") +
  xlab("Task") +
  ggsave("Hit_rate.png", width = 6, height = 3.5)


#Cognitive Load

# makes a dataframe with the mean reaction time to each stimulus
# for each condition (Single, SuRT, N-back) at the participant level
d2 <- cognitive %>%
  group_by(sub, task, stim) %>%
  summarise(rt_m = mean(mean_rt)) %>%
  filter(!is.na(rt_m)) %>%
  ungroup() %>%
  mutate(task = as.factor(task),
         stim = as.factor(stim),
         sub = as.factor(sub))

# linear mixed effects model
model2 <- lme(rt_m ~ task * stim, random = ~1|sub/task/stim, data = d2)
anova(model2)
summary(glht(model2, linfct = mcp(task = "Tukey")), test = adjusted(type="bonferroni"))

# boxplot of Reaction Time
d2 %>%
  mutate(Task = fct_recode(task, `N-back` = "nBack"),
         Stimulus = fct_recode(stim, `Remote LED` = "LED_dash",
                               `Vibrotactor` = "tactor")) %>%
  ggplot(mapping = aes(x = Task, y = rt_m, fill = Stimulus)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#f0f0f0", "#999999")) +
  theme_bw() +
  ylab("Mean Reaction Time") +
  xlab("Task") +
  ggsave("Reaction_time_test.png", width = 6, height = 3.5)



# Effect size analysis at the participant level -----------------
library(tidyr)

# Calculates Cohen's d for RT to the tactor
cohen.d_RTtactor <- cognitive %>%
  filter(task != "SuRT", stim != "LED_dash", n_cars >= 3) %>%
  spread(task, mean_rt) %>%
  filter(!is.na(nBack), !is.na(Single)) %>%
  group_by(sub, stim) %>%
  summarise(Single_mean = mean(Single), 
            nBack_mean = mean(nBack),
            Single_SD2 = sd(Single)^2,
            nBack_SD2 = sd(nBack)^2) %>%
  mutate(sd.pool = sqrt((Single_SD2 + nBack_SD2) / 2)) %>%
  mutate(cohen.d = (nBack_mean - Single_mean) / sd.pool) %>%
  # this filters out participants who only participated once
  filter(!is.na(sd.pool))

# Calculates Cohen's d for RT to the LED
cohen.d_RTLED <- cognitive %>%
  filter(task != "SuRT", stim != "tactor", n_cars >= 3) %>%
  spread(task, mean_rt) %>%
  filter(!is.na(nBack), !is.na(Single)) %>%
  group_by(sub, stim) %>%
  summarise(Single_mean = mean(Single), 
            nBack_mean = mean(nBack),
            Single_SD2 = sd(Single)^2,
            nBack_SD2 = sd(nBack)^2) %>%
  mutate(sd.pool = sqrt((Single_SD2 + nBack_SD2) / 2)) %>%
  mutate(cohen.d = (nBack_mean - Single_mean) / sd.pool) %>%
  # this filters out participants who only participated once
  filter(!is.na(sd.pool))

# Combine these two datasets and visualize them together in a boxplot
cohen.d_cognitive <- rbind(cohen.d_RTLED, cohen.d_RTtactor)

cohen.d_cognitive %>%
  mutate(Stimulus = fct_recode(stim, `Remote LED` = "LED_dash",
                               `Vibro-Tactor` = "tactor")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = Stimulus, y = cohen.d, fill = Stimulus), width = .3) +
  scale_fill_manual(values = c("#f0f0f0", "#999999")) +
  ylab("Effect Size (Cohen's d)") +
  theme_bw()

# Run a paired samples t test to see if the reaction time for each stimulus 
# is significantly different
cd_cognitive_ttest <- t.test(cohen.d ~ stim, cohen.d_cognitive)

      # supplies sd and mean for reporting alongside t test results
cohen.d_cognitive %>%
  dplyr::select(stim, cohen.d) %>%
  group_by(stim) %>%
  summarise(cohen.d_sd = sd(cohen.d),
            cohen.d_mean = mean(cohen.d))

cd_cognitive_ttest


# Calculates Cohen's d for Hit rate to the tactor
cohen.d_HRtactor <- visual %>%
  filter(task != "nBack", stim != "LED_dash", n_cars > 3) %>%
  spread(task, hit_rate) %>%
  filter(!is.na(SuRT), !is.na(Single)) %>%
  group_by(sub, stim) %>%
  summarise(Single_mean = mean(Single), 
            SuRT_mean = mean(SuRT),
            Single_SD2 = sd(Single)^2,
            SuRT_SD2 = sd(SuRT)^2) %>%
  mutate(sd.pool = sqrt((Single_SD2 + SuRT_SD2) / 2)) %>%
  mutate(cohen.d = abs((SuRT_mean - Single_mean) / sd.pool)) 

# Calculates Cohen's d for Hit rate to the LED
cohen.d_HRLED <- visual %>%
  filter(task != "nBack", stim != "tactor", n_cars > 3) %>%
  spread(task, hit_rate) %>%
  filter(!is.na(SuRT), !is.na(Single)) %>%
  group_by(sub, stim) %>%
  summarise(Single_mean = mean(Single), 
            SuRT_mean = mean(SuRT),
            Single_SD2 = sd(Single)^2,
            SuRT_SD2 = sd(SuRT)^2) %>%
  mutate(sd.pool = sqrt((Single_SD2 + SuRT_SD2) / 2)) %>%
  mutate(cohen.d = abs((SuRT_mean - Single_mean) / sd.pool)) 

# Combine these two datasets and visualize them together in a boxplot
cohen.d_visual <- rbind(cohen.d_HRLED, cohen.d_HRtactor)

cohen.d_visual %>%
  mutate(Stimulus = fct_recode(stim, `Remote LED` = "LED_dash",
                               `Vibro-Tactor` = "tactor")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = Stimulus, y = cohen.d, fill = Stimulus), width = .3) +
  scale_fill_manual(values = c("#f0f0f0", "#999999")) +
  ylab("Effect Size (Cohen's d)") +
  theme_bw()

# Run a paired samples t test to see if the reaction time for each stimulus 
# is significantly different
cd_visual_ttest <- t.test(cohen.d ~ stim, cohen.d_visual)

cohen.d_visual %>%
  group_by(stim) %>%
  dplyr::select(stim, cohen.d) %>%
  summarise(cohen.d_sd = sd(cohen.d, na.rm = T),
            cohen.d_mean = mean(cohen.d, na.rm = T))

cd_visual_ttest
