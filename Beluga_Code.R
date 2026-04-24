# Beluga Whales (Delphinapterus leucas) Discriminate Individual Callers Across Multiple Call Types
# 2026
#
# Author: 
# Software: RStudio Version 2025.09.2+418


# DESCRIPTION
# This script reproduces the statistical analyses for the figures, table, and the in-text results for the number of looks 
# toward the underwater speaker after the playback had occurred. For each beluga call type, a generalized linear mixed model (GLMM) 
# is fit with session_order as a fixed effect and target identity, caller identity, and trial number as random effects. 
# Poisson is used by default; if check_overdispersion() returns as over or underdispersed, the model is refit as negative binomial (nbinom2) 
# and the better fitting model (lower AIC, less over/underdispersion) is retained as the final model. Pairwise session contrasts 
# are extracted with emmeans on the response scale (back-transformed rate ratios) with Tukey adjustment. Cohen's d and 95% CIs are 
# computed for the three key contrasts (0 vs 1, 1 vs 3, 3 vs 4) using the effsize package with paired = TRUE on the raw look counts.


# REPRODUCIBILITY NOTE
# All seven call types use the same model structure. 
# The choice between Poisson and negative binomial is data-driven and recorded in the `final_family` object so the user can verify which family was used per call type.

#Load in library ----------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(emmeans)
library(performance)
library(effsize)
library(writexl)


# Load in dataset ---------------------------------------------------------

# Update the read excel path to the location of the dataset on your device

beluga <- read_excel("Downloads/Beluga_Dataset_Schoenhoft_et_al_2026.xlsx")

View(beluga)

# Create Factors ----------------------------------------------------------

beluga$session_order <- factor(beluga$session_order)
beluga$session_order <- relevel(beluga$session_order, ref = "0")
beluga$target <- factor(beluga$target)
beluga$sod <- factor(beluga$sod)
beluga$tod <- factor(beluga$tod)
beluga$year <- factor(beluga$year)
beluga$caller_ID <- factor(beluga$caller_ID)
beluga$trial_number <- factor(beluga$trial_number)
beluga$playback_type <- factor(beluga$playback_type)

# Sort data by Call Type --------------------------------------------------

callerccc   <- beluga %>% filter(playback_type == "CCC")
callerscc   <- beluga %>% filter(playback_type == "SCC")
callerbipt  <- beluga %>% filter(playback_type == "BiPT")
callersegw  <- beluga %>% filter(playback_type == "SegW")
callerbirws <- beluga %>% filter(playback_type == "BiRWS")
callerfw    <- beluga %>% filter(playback_type == "FW")
callertbc   <- beluga %>% filter(playback_type == "TBC")

# Playback Analysis --------------------------------------------------

# Each call type follows the same structure:
#   Step 1: Plot
#   Step 2: Fit poisson GLMM
#   Step 3: Check overdispersion
#   Step 4: Fit negative binomial in case poisson does not fit well
#   Step 5: Compare AIC; select final model
#   Step 6: Pairwise session contrasts on the response scale (Tukey adjusted)
#   Step 7: Cohen's d and 95% CI for the three key contrasts

# Complex Contact Call (CCC) ----------------------------------------------
 
# Step 1: Plot
ggplot(callerccc, aes(x = session_order, y = num_look)) +
  geom_boxplot(fill = "powderblue", color = "black") +
  labs(title = "CCC: Number of Looks by Session Order",
       x = "Session Order",
       y = "Number of Looks") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 
# Step 2: Fit poisson GLMM
ccc_pois <- glmer(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data    = callerccc,
  family  = poisson(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
 
summary(ccc_pois)
 
# Step 3: Check overdispersion
check_overdispersion(ccc_pois)

#0.96
#No overdispersion, keep poisson model
 
# Step 4: Fit negative binomial in case Poisson does not fit well
ccc_nb <- glmmTMB(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data   = callerccc,
  family = nbinom2)
 
summary(ccc_nb)
check_overdispersion(ccc_nb)
 
# Step 5: Compare AIC
AIC(ccc_pois, ccc_nb)
 
# Final model
ccc_final <- ccc_pois
 
# Step 6: Pairwise session contrasts on the response scale (Tukey adjusted)
ccc_emm   <- emmeans(ccc_final, ~ session_order, type = "response")

ccc_pairs <- as.data.frame(summary(pairs(ccc_emm, adjust = "tukey")))
if ("ratio" %in% names(ccc_pairs)) ccc_pairs <- rename(ccc_pairs, estimate = ratio)
ccc_pairs$playback_type <- "CCC"
ccc_pairs
 
# Step 7: Cohen's d and 95% CI for the three key contrasts (paired by trial)

# Contrast 0 - 1
ccc_wide_01 <- callerccc %>%
  filter(session_order %in% c("0", "1")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

ccc_d01_obj <- effsize::cohen.d(ccc_wide_01$s1, ccc_wide_01$s0, paired = TRUE)
 
# Contrast 1 - 3
ccc_wide_13 <- callerccc %>%
  filter(session_order %in% c("1", "3")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

ccc_d13_obj <- effsize::cohen.d(ccc_wide_13$s3, ccc_wide_13$s1, paired = TRUE)
 
# Contrast 3 - 4
ccc_wide_34 <- callerccc %>%
  filter(session_order %in% c("3", "4")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

ccc_d34_obj <- effsize::cohen.d(ccc_wide_34$s4, ccc_wide_34$s3, paired = TRUE)
 
ccc_d <- data.frame(
  playback_type = "CCC",
  contrast = c("session_order0 - session_order1",
               "session_order1 - session_order3",
               "session_order3 - session_order4"),
  cohens_d = c(as.numeric(ccc_d01_obj$estimate),
               as.numeric(ccc_d13_obj$estimate),
               as.numeric(ccc_d34_obj$estimate)),
  ci_low = c(ccc_d01_obj$conf.int[["lower"]],
             ccc_d13_obj$conf.int[["lower"]],
             ccc_d34_obj$conf.int[["lower"]]),
  ci_high = c(ccc_d01_obj$conf.int[["upper"]],
              ccc_d13_obj$conf.int[["upper"]],
              ccc_d34_obj$conf.int[["upper"]]),
  stringsAsFactors = FALSE)

ccc_d
 

# Simple Contact Call (SCC) -----------------------------------------------
 
# Step 1: Plot
ggplot(callerscc, aes(x = session_order, y = num_look)) +
  geom_boxplot(fill = "thistle", color = "black") +
  labs(title = "SCC: Number of Looks by Session Order",
       x = "Session Order",
       y = "Number of Looks") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 
# Step 2: Poisson GLMM
scc_pois <- glmer(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data    = callerscc,
  family  = poisson(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
 
summary(scc_pois)
 
# Step 3: Check overdispersion
check_overdispersion(scc_pois)

#0.93
#No overdispersion, keep poisson
 
# Step 4: Negative binomial
scc_nb <- glmmTMB(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data   = callerscc,
  family = nbinom2)
 
summary(scc_nb)
check_overdispersion(scc_nb)
 
# Step 5: Compare AIC
AIC(scc_pois, scc_nb)
 
# Final model
scc_final <- scc_pois
 
# Step 6: Pairwise session contrasts on the response scale
scc_emm   <- emmeans(scc_final, ~ session_order, type = "response")

scc_pairs <- as.data.frame(summary(pairs(scc_emm, adjust = "tukey")))
if ("ratio" %in% names(scc_pairs)) scc_pairs <- rename(scc_pairs, estimate = ratio)
scc_pairs$playback_type <- "SCC"
scc_pairs
 
# Step 7: Cohen's d and 95% CI
scc_wide_01 <- callerscc %>%
  filter(session_order %in% c("0", "1")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()
scc_d01_obj <- effsize::cohen.d(scc_wide_01$s1, scc_wide_01$s0, paired = TRUE)
 
scc_wide_13 <- callerscc %>%
  filter(session_order %in% c("1", "3")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()
scc_d13_obj <- effsize::cohen.d(scc_wide_13$s3, scc_wide_13$s1, paired = TRUE)
 
scc_wide_34 <- callerscc %>%
  filter(session_order %in% c("3", "4")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()
scc_d34_obj <- effsize::cohen.d(scc_wide_34$s4, scc_wide_34$s3, paired = TRUE)
 
scc_d <- data.frame(
  playback_type = "SCC",
  contrast = c("session_order0 - session_order1",
               "session_order1 - session_order3",
               "session_order3 - session_order4"),
  cohens_d = c(as.numeric(scc_d01_obj$estimate),
               as.numeric(scc_d13_obj$estimate),
               as.numeric(scc_d34_obj$estimate)),
  ci_low = c(scc_d01_obj$conf.int[["lower"]],
             scc_d13_obj$conf.int[["lower"]],
             scc_d34_obj$conf.int[["lower"]]),
  ci_high = c(scc_d01_obj$conf.int[["upper"]],
              scc_d13_obj$conf.int[["upper"]],
              scc_d34_obj$conf.int[["upper"]]),
  stringsAsFactors = FALSE)

scc_d
 

# Biphonated Pulsed Tonal (BiPT) ------------------------------------------
 
# Step 1: Plot
ggplot(callerbipt, aes(x = session_order, y = num_look)) +
  geom_boxplot(fill = "honeydew2", color = "black") +
  labs(title = "BiPT: Number of Looks by Session Order",
       x = "Session Order",
       y = "Number of Looks") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 
# Step 2: Poisson GLMM
bipt_pois <- glmer(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data    = callerbipt,
  family  = poisson(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
 
summary(bipt_pois)
 
# Step 3: Check overdispersion
check_overdispersion(bipt_pois)

#1.05
#Little to no overdispersion ≈1, keep poisson
 
# Step 4: Negative binomial
bipt_nb <- glmmTMB(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data   = callerbipt,
  family = nbinom2)
 
summary(bipt_nb)
check_overdispersion(bipt_nb)
 
# Step 5: Compare AIC
AIC(bipt_pois, bipt_nb)
 
# Final model
bipt_final <- bipt_pois
 
# Step 6: Pairwise session contrasts on the response scale
bipt_emm   <- emmeans(bipt_final, ~ session_order, type = "response")

bipt_pairs <- as.data.frame(summary(pairs(bipt_emm, adjust = "tukey")))
if ("ratio" %in% names(bipt_pairs)) bipt_pairs <- rename(bipt_pairs, estimate = ratio)
bipt_pairs$playback_type <- "BiPT"
bipt_pairs
 
# Step 7: Cohen's d and 95% CI
bipt_wide_01 <- callerbipt %>%
  filter(session_order %in% c("0", "1")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

bipt_d01_obj <- effsize::cohen.d(bipt_wide_01$s1, bipt_wide_01$s0, paired = TRUE)
 
bipt_wide_13 <- callerbipt %>%
  filter(session_order %in% c("1", "3")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

bipt_d13_obj <- effsize::cohen.d(bipt_wide_13$s3, bipt_wide_13$s1, paired = TRUE)
 
bipt_wide_34 <- callerbipt %>%
  filter(session_order %in% c("3", "4")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

bipt_d34_obj <- effsize::cohen.d(bipt_wide_34$s4, bipt_wide_34$s3, paired = TRUE)
 
bipt_d <- data.frame(
  playback_type = "BiPT",
  contrast = c("session_order0 - session_order1",
               "session_order1 - session_order3",
               "session_order3 - session_order4"),
  cohens_d = c(as.numeric(bipt_d01_obj$estimate),
               as.numeric(bipt_d13_obj$estimate),
               as.numeric(bipt_d34_obj$estimate)),
  ci_low = c(bipt_d01_obj$conf.int[["lower"]],
             bipt_d13_obj$conf.int[["lower"]],
             bipt_d34_obj$conf.int[["lower"]]),
  ci_high = c(bipt_d01_obj$conf.int[["upper"]],
              bipt_d13_obj$conf.int[["upper"]],
              bipt_d34_obj$conf.int[["upper"]]),
  stringsAsFactors = FALSE)

bipt_d
 
# Segmented Whistle (SegW) ------------------------------------------------
 
# Step 1: Plot
ggplot(callersegw, aes(x = session_order, y = num_look)) +
  geom_boxplot(fill = "palevioletred", color = "black") +
  labs(title = "SegW: Number of Looks by Session Order",
       x = "Session Order",
       y = "Number of Looks") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 
# Step 2: Poisson GLMM
segw_pois <- glmer(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data    = callersegw,
  family  = poisson(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
 
summary(segw_pois)
 
# Step 3: Check overdispersion
check_overdispersion(segw_pois)

#0.85
#mild underdispersion, still within normal limits, keep poisson
 
# Step 4: Negative binomial
segw_nb <- glmmTMB(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data   = callersegw,
  family = nbinom2)
 
summary(segw_nb)
check_overdispersion(segw_nb)
 
# Step 5: Compare AIC
AIC(segw_pois, segw_nb)
 
# Final model
segw_final <- segw_pois
 
# Step 6: Pairwise session contrasts on the response scale
segw_emm   <- emmeans(segw_final, ~ session_order, type = "response")

segw_pairs <- as.data.frame(summary(pairs(segw_emm, adjust = "tukey")))
if ("ratio" %in% names(segw_pairs)) segw_pairs <- rename(segw_pairs, estimate = ratio)
segw_pairs$playback_type <- "SegW"
segw_pairs
 
# Step 7: Cohen's d and 95% CI
segw_wide_01 <- callersegw %>%
  filter(session_order %in% c("0", "1")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

segw_d01_obj <- effsize::cohen.d(segw_wide_01$s1, segw_wide_01$s0, paired = TRUE)
 
segw_wide_13 <- callersegw %>%
  filter(session_order %in% c("1", "3")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

segw_d13_obj <- effsize::cohen.d(segw_wide_13$s3, segw_wide_13$s1, paired = TRUE)
 
segw_wide_34 <- callersegw %>%
  filter(session_order %in% c("3", "4")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

segw_d34_obj <- effsize::cohen.d(segw_wide_34$s4, segw_wide_34$s3, paired = TRUE)
 
segw_d <- data.frame(
  playback_type = "SegW",
  contrast = c("session_order0 - session_order1",
               "session_order1 - session_order3",
               "session_order3 - session_order4"),
  cohens_d = c(as.numeric(segw_d01_obj$estimate),
               as.numeric(segw_d13_obj$estimate),
               as.numeric(segw_d34_obj$estimate)),
  ci_low = c(segw_d01_obj$conf.int[["lower"]],
             segw_d13_obj$conf.int[["lower"]],
             segw_d34_obj$conf.int[["lower"]]),
  ci_high = c(segw_d01_obj$conf.int[["upper"]],
              segw_d13_obj$conf.int[["upper"]],
              segw_d34_obj$conf.int[["upper"]]),
  stringsAsFactors = FALSE)

segw_d
 
# Biphonated r-Shaped Whistle (BiRWS) -------------------------------------
 
# Step 1: Plot
ggplot(callerbirws, aes(x = session_order, y = num_look)) +
  geom_boxplot(fill = "goldenrod1", color = "black") +
  labs(title = "BiRWS: Number of Looks by Session Order",
       x = "Session Order",
       y = "Number of Looks") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 
# Step 2: Poisson GLMM
birws_pois <- glmer(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data    = callerbirws,
  family  = poisson(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
 
summary(birws_pois)
 
# Step 3: Check overdispersion
check_overdispersion(birws_pois)

#0.94
#No overdispersion, keep poisson
 
# Step 4: Negative binomial
birws_nb <- glmmTMB(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data   = callerbirws,
  family = nbinom2)
 
summary(birws_nb)
check_overdispersion(birws_nb)
 
# Step 5: Compare AIC
AIC(birws_pois, birws_nb)
 
# Final model
birws_final <- birws_pois
 
# Step 6: Pairwise session contrasts on the response scale
birws_emm   <- emmeans(birws_final, ~ session_order, type = "response")

birws_pairs <- as.data.frame(summary(pairs(birws_emm, adjust = "tukey")))
if ("ratio" %in% names(birws_pairs)) birws_pairs <- rename(birws_pairs, estimate = ratio)
birws_pairs$playback_type <- "BiRWS"
birws_pairs
 
# Step 7: Cohen's d and 95% CI
birws_wide_01 <- callerbirws %>%
  filter(session_order %in% c("0", "1")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

birws_d01_obj <- effsize::cohen.d(birws_wide_01$s1, birws_wide_01$s0, paired = TRUE)
 
birws_wide_13 <- callerbirws %>%
  filter(session_order %in% c("1", "3")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

birws_d13_obj <- effsize::cohen.d(birws_wide_13$s3, birws_wide_13$s1, paired = TRUE)
 
birws_wide_34 <- callerbirws %>%
  filter(session_order %in% c("3", "4")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

birws_d34_obj <- effsize::cohen.d(birws_wide_34$s4, birws_wide_34$s3, paired = TRUE)
 
birws_d <- data.frame(
  playback_type = "BiRWS",
  contrast = c("session_order0 - session_order1",
               "session_order1 - session_order3",
               "session_order3 - session_order4"),
  cohens_d = c(as.numeric(birws_d01_obj$estimate),
               as.numeric(birws_d13_obj$estimate),
               as.numeric(birws_d34_obj$estimate)),
  ci_low = c(birws_d01_obj$conf.int[["lower"]],
             birws_d13_obj$conf.int[["lower"]],
             birws_d34_obj$conf.int[["lower"]]),
  ci_high = c(birws_d01_obj$conf.int[["upper"]],
              birws_d13_obj$conf.int[["upper"]],
              birws_d34_obj$conf.int[["upper"]]),
  stringsAsFactors = FALSE)

birws_d
 
# Flat Whistle (FW) -------------------------------------------------------
 
# Step 1: Plot
ggplot(callerfw, aes(x = session_order, y = num_look)) +
  geom_boxplot(fill = "palegreen3", color = "black") +
  labs(title = "FW: Number of Looks by Session Order",
       x = "Session Order",
       y = "Number of Looks") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 
# Step 2: Poisson GLMM
fw_pois <- glmer(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data    = callerfw,
  family  = poisson(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
 
summary(fw_pois)
 
# Step 3: Check overdispersion
check_overdispersion(fw_pois)

#Mild underdispersion, keep poisson
 
# Step 4: Negative binomial
fw_nb <- glmmTMB(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data   = callerfw,
  family = nbinom2)
 
summary(fw_nb)
check_overdispersion(fw_nb)
 
# Step 5: Compare AIC
AIC(fw_pois, fw_nb)
 
# Final model
fw_final <- fw_pois
 
# Step 6: Pairwise session contrasts on the response scale
fw_emm   <- emmeans(fw_final, ~ session_order, type = "response")

fw_pairs <- as.data.frame(summary(pairs(fw_emm, adjust = "tukey")))
if ("ratio" %in% names(fw_pairs)) fw_pairs <- rename(fw_pairs, estimate = ratio)
fw_pairs$playback_type <- "FW"
fw_pairs
 
# Step 7: Cohen's d and 95% CI
fw_wide_01 <- callerfw %>%
  filter(session_order %in% c("0", "1")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

fw_d01_obj <- effsize::cohen.d(fw_wide_01$s1, fw_wide_01$s0, paired = TRUE)
 
fw_wide_13 <- callerfw %>%
  filter(session_order %in% c("1", "3")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

fw_d13_obj <- effsize::cohen.d(fw_wide_13$s3, fw_wide_13$s1, paired = TRUE)
 
fw_wide_34 <- callerfw %>%
  filter(session_order %in% c("3", "4")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

fw_d34_obj <- effsize::cohen.d(fw_wide_34$s4, fw_wide_34$s3, paired = TRUE)
 
fw_d <- data.frame(
  playback_type = "FW",
  contrast = c("session_order0 - session_order1",
               "session_order1 - session_order3",
               "session_order3 - session_order4"),
  cohens_d = c(as.numeric(fw_d01_obj$estimate),
               as.numeric(fw_d13_obj$estimate),
               as.numeric(fw_d34_obj$estimate)),
  ci_low = c(fw_d01_obj$conf.int[["lower"]],
             fw_d13_obj$conf.int[["lower"]],
             fw_d34_obj$conf.int[["lower"]]),
  ci_high = c(fw_d01_obj$conf.int[["upper"]],
              fw_d13_obj$conf.int[["upper"]],
              fw_d34_obj$conf.int[["upper"]]),
  stringsAsFactors = FALSE)

fw_d
 

# Tonal Band Change (TBC) -------------------------------------------------
 
# Step 1: Plot
ggplot(callertbc, aes(x = session_order, y = num_look)) +
  geom_boxplot(fill = "darkslategray3", color = "black") +
  labs(title = "TBC: Number of Looks by Session Order",
       x = "Session Order",
       y = "Number of Looks") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 
# Step 2: Poisson GLMM
tbc_pois <- glmer(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data    = callertbc,
  family  = poisson(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
 
summary(tbc_pois)
 
# Step 3: Check overdispersion
check_overdispersion(tbc_pois)

#1.6
#Overdispersed, switch to NB
 
# Step 4: Negative binomial
tbc_nb <- glmmTMB(
  num_look ~ session_order + (1 | target) + (1 | caller_ID) + (1 | trial_number),
  data   = callertbc,
  family = nbinom2)
 
summary(tbc_nb)
check_overdispersion(tbc_nb)

#1.36
#Mild overdisperson, keep NB
 
# Step 5: Compare AIC
AIC(tbc_pois, tbc_nb)
 
# Final model
tbc_final <- tbc_nb
 
# Step 6: Pairwise session contrasts on the response scale
tbc_emm   <- emmeans(tbc_final, ~ session_order, type = "response")

tbc_pairs <- as.data.frame(summary(pairs(tbc_emm, adjust = "tukey")))
if ("ratio" %in% names(tbc_pairs)) tbc_pairs <- rename(tbc_pairs, estimate = ratio)
tbc_pairs$playback_type <- "TBC"
tbc_pairs
 
# Step 7: Cohen's d and 95% CI
tbc_wide_01 <- callertbc %>%
  filter(session_order %in% c("0", "1")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

tbc_d01_obj <- effsize::cohen.d(tbc_wide_01$s1, tbc_wide_01$s0, paired = TRUE)
 
tbc_wide_13 <- callertbc %>%
  filter(session_order %in% c("1", "3")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

tbc_d13_obj <- effsize::cohen.d(tbc_wide_13$s3, tbc_wide_13$s1, paired = TRUE)
 
tbc_wide_34 <- callertbc %>%
  filter(session_order %in% c("3", "4")) %>%
  dplyr::select(trial_number, session_order, num_look) %>%
  pivot_wider(names_from = session_order, values_from = num_look,
              names_prefix = "s") %>%
  drop_na()

tbc_d34_obj <- effsize::cohen.d(tbc_wide_34$s4, tbc_wide_34$s3, paired = TRUE)
 
tbc_d <- data.frame(
  playback_type = "TBC",
  contrast = c("session_order0 - session_order1",
               "session_order1 - session_order3",
               "session_order3 - session_order4"),
  cohens_d = c(as.numeric(tbc_d01_obj$estimate),
               as.numeric(tbc_d13_obj$estimate),
               as.numeric(tbc_d34_obj$estimate)),
  ci_low = c(tbc_d01_obj$conf.int[["lower"]],
             tbc_d13_obj$conf.int[["lower"]],
             tbc_d34_obj$conf.int[["lower"]]),
  ci_high = c(tbc_d01_obj$conf.int[["upper"]],
              tbc_d13_obj$conf.int[["upper"]],
              tbc_d34_obj$conf.int[["upper"]]),
  stringsAsFactors = FALSE)

tbc_d
 
# Combine All Session Contrasts -------------------------------------------
summary_contrasts <- bind_rows(
  ccc_pairs, scc_pairs, bipt_pairs, segw_pairs,
  birws_pairs, fw_pairs, tbc_pairs) %>%
  dplyr::select(playback_type, contrast, estimate, SE, df, p.value)
 
summary_contrasts
 
# Format every contrast as "estimate ± SE, p = p value"
results <- summary_contrasts %>%
  mutate(formatted = paste0(
    sprintf("%.2f", estimate), " ± ", sprintf("%.2f", SE),
    ", p=",
    ifelse(
      p.value < 0.0001,
      format(p.value, scientific = TRUE, digits = 2),
      sprintf("%.4f", p.value))))

#See results
results

results <- results %>%
 mutate(
    contrast = dplyr::recode(
      contrast,
      "session_order0 - session_order1" = "0 to 1",
      "session_order0 - session_order2" = "0 to 2",
      "session_order0 - session_order3" = "0 to 3",
      "session_order0 - session_order4" = "0 to 4",
      "session_order1 - session_order2" = "1 to 2",
      "session_order1 - session_order3" = "1 to 3",
      "session_order1 - session_order4" = "1 to 4",
      "session_order2 - session_order3" = "2 to 3",
      "session_order2 - session_order4" = "2 to 4",
      "session_order3 - session_order4" = "3 to 4"))

 
# Split results into one data frame per call type
results_by_playtype <- split(results, results$playback_type)
 
results_by_playtype$CCC
results_by_playtype$SCC
results_by_playtype$BiPT
results_by_playtype$SegW
results_by_playtype$BiRWS
results_by_playtype$FW
results_by_playtype$TBC
 
# Reshape into table layout: rows = contrasts, columns = call types
table_ready <- results %>%
  dplyr::select(playback_type, contrast, formatted) %>%
  pivot_wider(names_from = playback_type, values_from = formatted) %>%
  dplyr::select(contrast, any_of(c("CCC", "SCC", "BiPT", "SegW", "BiRWS", "TBC", "FW")))
 
#Relabel for ease of reading
colnames(table_ready)[1] <- "Session Contrasts"

table_ready
 
# Key Contrasts with Effect Sizes -----------------------------------------

# Combine model estimates with Cohen's d and 95% CI for the main three contrasts that drive our interpretations:
# 0 - 1  (control to first stimulus)        habituation onset
# 1 - 3  (across repeated stimuli)          habituation
# 3 - 4  (caller change)                    dishabituation / discrimination
 
key_contrasts <- c(
  "session_order0 - session_order1",
  "session_order1 - session_order3",
  "session_order3 - session_order4")
 
key_emmeans_results <- summary_contrasts %>%
  filter(contrast %in% key_contrasts)
 
key_effect_sizes <- bind_rows(
  ccc_d, scc_d, bipt_d, segw_d, birws_d, fw_d, tbc_d)
 
key_effect_sizes
 
key_results_final <- key_emmeans_results %>%
  left_join(key_effect_sizes, by = c("playback_type", "contrast")) %>%
  mutate(formatted_result = paste0(
    round(estimate, 2), " ± ", round(SE, 2),
    ", p = ",
    ifelse(p.value < 0.0001,
           format(p.value, scientific = TRUE, digits = 3),
           signif(p.value, 3)),
    ", d = ", round(cohens_d, 2),
    ", 95% CI [", round(ci_low, 2), ", ", round(ci_high, 2), "]"))
 
key_results_final
 
 
# Export to Excel ---------------------------------------------------------

write_xlsx(
  list(
    all_session_contrasts = summary_contrasts,
    table2                = table_ready,
    key_results_main_text = key_results_final),
  "Beluga_Caller_Results.xlsx")

#See where R exported excel file
getwd()
