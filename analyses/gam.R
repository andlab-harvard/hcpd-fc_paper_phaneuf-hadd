##### Conduct HCPD FC GAM Analyses #####

# Written by: Camille Phaneuf-Hadd (cphaneuf@g.harvard.edu)
# Last updated: 12/10/25
# Question: how do different measures of FC change across development?

##### Set up Script ##### 

# Load needed libraries
require(pacman) # for p_load()
p_load(tidyverse, # for df manipulation
       dplyr, # for %>% and other operators
       ggplot2, # for plotting
       sjPlot, # for plot_model()
       mgcv, # for GAMs
       gratia, # for appraise()
       itsadug, # for plot_diff()
       ggeffects, # for ggpredict()
       gratia) # for derivatives()

# Load shared HCPD FC functions
source("utilities.R")

# Set data directory paths
getwd()
processed_data_dir <- "../data/"
analyzed_data_dir <- "../results/gam/"
factor_data_dir <- "../results/psychom/"

# Read in demographic and behavioral data
uncorr_toolbox <- read.csv(paste0(processed_data_dir, "uncorr_toolbox.csv"))

# Read in factor score data
cfa_2_scores <- read.csv(paste0(factor_data_dir, "cfa_2factor_scores.csv"))
cfa_2_scores$X <- NULL # Remove extraneous column

# Convert uncorrected toolbox data into long format
fluid_measures_uncorr <- uncorr_toolbox[, c("Subject", "Age", toolbox_z)]
names(fluid_measures_uncorr) <- c("Subject", "Age", long_task_labels)
fluid_measures_uncorr_long <- fluid_measures_uncorr %>% pivot_longer(cols = long_task_labels, names_to = "Measure", "values_to" = "Score")
fluid_measures_uncorr_long$Measure <- factor(fluid_measures_uncorr_long$Measure, levels = long_task_labels)

# Compare uncorrected toolbox scores to the mean of every other
fluid_measures_uncorr_means <- fluid_measures_uncorr
fluid_measures_uncorr_means$NoFlanker <- c(-999)
fluid_measures_uncorr_means$NoDCCS <- c(-999)
fluid_measures_uncorr_means$NoPattern <- c(-999)
fluid_measures_uncorr_means$NoPSMT <- c(-999)
fluid_measures_uncorr_means$NoList <- c(-999)
for (i in 1:length(fluid_measures_uncorr_means$Subject)) {
  flanker <- fluid_measures_uncorr_means$Flanker[i]
  dccs <- fluid_measures_uncorr_means$DCCS[i]
  pattern <- fluid_measures_uncorr_means$PatternComparison[i]
  psmt <- fluid_measures_uncorr_means$PSMT[i]
  list <- fluid_measures_uncorr_means$ListSorting[i]
  
  # No Flanker
  no_flanker_sum <- sum(c(dccs, pattern, psmt, list), na.rm = TRUE)
  no_flanker_n <- sum(!is.na(c(dccs, pattern, psmt, list)))
  fluid_measures_uncorr_means$NoFlanker[i] <- no_flanker_sum / no_flanker_n
  
  # No DCCS
  no_dccs_sum <- sum(c(flanker, pattern, psmt, list), na.rm = TRUE)
  no_dccs_n <- sum(!is.na(c(flanker, pattern, psmt, list)))
  fluid_measures_uncorr_means$NoDCCS[i] <- no_dccs_sum / no_dccs_n
  
  # No Pattern Comparison
  no_pattern_sum <- sum(c(flanker, dccs, psmt, list), na.rm = TRUE)
  no_pattern_n <- sum(!is.na(c(flanker, dccs, psmt, list)))
  fluid_measures_uncorr_means$NoPattern[i] <- no_pattern_sum / no_pattern_n
  
  # No PSMT
  no_psmt_sum <- sum(c(flanker, dccs, pattern, list), na.rm = TRUE)
  no_psmt_n <- sum(!is.na(c(flanker, dccs, pattern, list)))
  fluid_measures_uncorr_means$NoPSMT[i] <- no_psmt_sum / no_psmt_n
  
  # No List Sorting
  no_list_sum <- sum(c(flanker, dccs, pattern, psmt), na.rm = TRUE)
  no_list_n <- sum(!is.na(c(flanker, dccs, pattern, psmt)))
  fluid_measures_uncorr_means$NoList[i] <- no_list_sum / no_list_n
}
sum(fluid_measures_uncorr_means$NoFlanker == -999) # Should be 0
sum(fluid_measures_uncorr_means$NoDCCS == -999) # Should be 0
sum(fluid_measures_uncorr_means$NoPattern == -999) # Should be 0
sum(fluid_measures_uncorr_means$NoPSMT == -999) # Should be 0
sum(fluid_measures_uncorr_means$NoList == -999) # Should be 0

# Convert pairwise mean data frames to long format
no_flanker <- fluid_measures_uncorr_means[, c("Subject", "Age", "Flanker", "NoFlanker")]
no_flanker_long <- no_flanker %>% pivot_longer(cols = c("Flanker", "NoFlanker"), names_to = "Measure", "values_to" = "Score")
no_flanker_long$Measure <- factor(no_flanker_long$Measure, levels = c("Flanker", "NoFlanker"))
no_flanker_long$Subject <- as.factor(no_flanker_long$Subject)
no_dccs <- fluid_measures_uncorr_means[, c("Subject", "Age", "DCCS", "NoDCCS")]
no_dccs_long <- no_dccs %>% pivot_longer(cols = c("DCCS", "NoDCCS"), names_to = "Measure", "values_to" = "Score")
no_dccs_long$Measure <- factor(no_dccs_long$Measure, levels = c("DCCS", "NoDCCS"))
no_dccs_long$Subject <- as.factor(no_dccs_long$Subject)
no_pattern <- fluid_measures_uncorr_means[, c("Subject", "Age", "PatternComparison", "NoPattern")]
no_pattern_long <- no_pattern %>% pivot_longer(cols = c("PatternComparison", "NoPattern"), names_to = "Measure", "values_to" = "Score")
no_pattern_long$Measure <- factor(no_pattern_long$Measure, levels = c("PatternComparison", "NoPattern"))
no_pattern_long$Subject <- as.factor(no_pattern_long$Subject)
no_psmt <- fluid_measures_uncorr_means[, c("Subject", "Age", "PSMT", "NoPSMT")]
no_psmt_long <- no_psmt %>% pivot_longer(cols = c("PSMT", "NoPSMT"), names_to = "Measure", "values_to" = "Score")
no_psmt_long$Measure <- factor(no_psmt_long$Measure, levels = c("PSMT", "NoPSMT"))
no_psmt_long$Subject <- as.factor(no_psmt_long$Subject)
no_list <- fluid_measures_uncorr_means[, c("Subject", "Age", "ListSorting", "NoList")]
no_list_long <- no_list %>% pivot_longer(cols = c("ListSorting", "NoList"), names_to = "Measure", "values_to" = "Score")
no_list_long$Measure <- factor(no_list_long$Measure, levels = c("ListSorting", "NoList"))
no_list_long$Subject <- as.factor(no_list_long$Subject)

##### Notes About GAMM Function and Parameter Decisions #####

# gam() syntax ---

# * gamm() or gam() are both appropriate because Score is normally distributed for every 
#   task and factor --> going with gam() to be consistent with Heffer et al. (2024)
# * method = REML mitigates over-fitting in wiggly splines and handles random effects well
# * random slope models fit with bam() (vs. gam()) because of convergence issues

# derivatives() syntax ---

# Use simultaneous CIs (interval = "simultaneous", unconditional = TRUE) for repeated tests against m = 0

# plot_diff() syntax ---

# * Use simultaneous CIs (sim.ci = TRUE) for consistency

##### Examine NIH Toolbox Data Across Age - Flanker ##### 

# Model Flanker scores across age
hist(fluid_measures_uncorr$Flanker) 
flanker_age_gam <- mgcv::gam(Flanker ~ s(Age), data = fluid_measures_uncorr, method = "REML") 
gam.check(flanker_age_gam) # k-index should be larger than 1 (no), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(flanker_age_gam) # Residuals should be normally distributed (yes)
summary(flanker_age_gam) # H0 = horizontal line is REJECTED; adjusted R2 is fine

# Identify the region of significant age-related change by computing the first derivative
# Credit (including for structure of all GAM plateau analyses to follow): Katherine Grisanzio 
# (code https://osf.io/fvy8d for paper https://andl.wjh.harvard.edu/files/2023/11/Grisanzio-et-al-2023-1.pdf)
flanker_deriv <- derivatives(flanker_age_gam, n = 100, eps = 1e-07, level = .95, interval = "simultaneous", unconditional = TRUE, n_sim = 10000) 
flanker_deriv_sig <- clip_on_siggratia(flanker_deriv)
flanker_deriv_sig$sig # 1s and 0s are continuous
flanker_deriv_plot_range <- range(flanker_deriv_sig[which(flanker_deriv_sig$sig == 1), "Age"])

# Plot flanker_age_gam effects
plot_model(flanker_age_gam, type = "pred", terms = "Age") # Flanker scores increase with age during childhood and adolescence, but then plateau into adulthood
flanker_age_gam_eff <- ggpredict(flanker_age_gam, terms = c("Age"))
ggplot(flanker_age_gam_eff, aes(x, predicted)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  scale_y_continuous(breaks = seq(-1.5, 1.5, by = .75)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = .75), colour = 'white') +
  annotate("rect", xmin = flanker_deriv_plot_range[1], xmax = flanker_deriv_plot_range[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = flanker_blue) +
  geom_line(size = 1.5, colour = flanker_blue) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .5, fill = flanker_blue) +
  labs(x = "Age (Years)", y = "Fitted Flanker Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_age_flanker.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - Flanker Contrasts #####

# Model uncorrected toolbox scores across age - Flanker compared to every other task
hist(no_flanker_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "flanker_gam.rda"))) {
  flanker_gam <- mgcv::gam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re"), data = no_flanker_long, method = "REML") 
  save(flanker_gam, file = paste0(analyzed_data_dir, "flanker_gam.rda"))
} else {
  load(paste0(analyzed_data_dir, "flanker_gam.rda"))
}
gam.check(flanker_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(flanker_gam) # Residuals should be normally distributed (yes)
summary(flanker_gam) # H0 = horizontal line is REJECTED for both levels of Measure; adjusted R2 is fine

# Plot flanker_gam effects
plot_model(flanker_gam, type = "pred", terms = "Age") # FC scores increase with age, especially during childhood
plot_model(flanker_gam, type = "pred", terms = "Measure") # Non-Flanker scores are higher than Flanker scores
plot_model(flanker_gam, type = "pred", terms = c("Age", "Measure"), show.data = TRUE) # Flanker and non-Flanker scores increase at similar rates with age
flanker_gam_eff <- ggpredict(flanker_gam, terms = c("Age", "Measure"))
ggplot(flanker_gam_eff, aes(x, predicted)) + 
  scale_y_continuous(breaks = seq(-2, 2, by = 1)) +
  geom_hline(yintercept = seq(-2, 2, by = 1), colour = 'white') +
  geom_line(aes(color = group), size = 1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .5) +
  scale_color_manual(values = c(flanker_blue, comp_yellow), labels = no_flanker_labels) +
  scale_fill_manual(values = c(flanker_blue, comp_yellow), labels = no_flanker_labels) +
  labs(x = "Age (Years)", y = "Fitted Scores", color = "Measure", fill = "Measure") +
  plot_theme + theme(legend.position = "bottom")
ggsave(paste0(analyzed_data_dir, 'plot_flanker_eff.png'), width = 7, height = 7)

# Plot significant region of flanker_gam
if (!file.exists(paste0(analyzed_data_dir, "flanker_gam_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  flanker_gam_plot <- plot_diff(flanker_gam, view = "Age", comp = list(Measure = c("Flanker", "NoFlanker")),
                                sim.ci = TRUE) # Flanker and non-Flanker scores never differ
  save(flanker_gam_plot, file = paste0(analyzed_data_dir, "flanker_gam_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "flanker_gam_plot.rda"))
}
ggplot(flanker_gam_plot, aes(x = Age, y = est)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +  
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = flanker_blue, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = flanker_blue) +
  labs(x = "Age (Years)", y = "Estimated Difference in\nFlanker - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_flanker.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - Flanker Contrasts with Random Slope #####

# Model uncorrected toolbox scores across age *with random slope* - Flanker compared to every other task
if (!file.exists(paste0(analyzed_data_dir, "flanker_gam_slopes.rda"))) {
  flanker_gam_slopes <- mgcv::bam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re") + s(Subject, Measure, bs = "re"),
                                  data = no_flanker_long,
                                  method = "fREML", # bam uses 'fREML' instead of 'REML'
                                  discrete = TRUE, # Speeds up computation
                                  nthreads = 4, # Use multiple cores if available
                                  verbose = TRUE)
  save(flanker_gam_slopes, file = paste0(analyzed_data_dir, "flanker_gam_slopes.rda"))
} else {
  load(paste0(analyzed_data_dir, "flanker_gam_slopes.rda"))
}
gam.check(flanker_gam_slopes) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(flanker_gam_slopes) # Residuals should be normally distributed (yes)
# summary(flanker_gam_slopes) # Causes R to hang
mgcv::gam.vcomp(flanker_gam_slopes) # Examine random effects directly

# Plot significant region of flanker_gam_slopes
if (!file.exists(paste0(analyzed_data_dir, "flanker_gam_slopes_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  flanker_gam_slopes_plot <- plot_diff(flanker_gam_slopes, view = "Age", comp = list(Measure = c("Flanker", "NoFlanker")),
                                       sim.ci = TRUE) # Flanker and non-Flanker scores never differ
  save(flanker_gam_slopes_plot, file = paste0(analyzed_data_dir, "flanker_gam_slopes_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "flanker_gam_slopes_plot.rda"))
}
ggplot(flanker_gam_slopes_plot, aes(x = Age, y = est)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +  
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = flanker_blue, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = flanker_blue) +
  labs(x = "Age (Years)", y = "Estimated Difference in\nFlanker - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_flanker_slopes.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - DCCS ##### 

# Model DCCS scores across age
hist(fluid_measures_uncorr$DCCS) 
dccs_age_gam <- mgcv::gam(DCCS ~ s(Age), data = fluid_measures_uncorr, method = "REML") 
gam.check(dccs_age_gam) # k-index should be larger than 1 (no), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(dccs_age_gam) # Residuals should be normally distributed (yes)
summary(dccs_age_gam) # H0 = horizontal line is REJECTED; adjusted R2 is fine

# Identify the region of significant age-related change by computing the first derivative
dccs_deriv <- derivatives(dccs_age_gam, n = 100, eps = 1e-07, level = .95, interval = "simultaneous", unconditional = TRUE, n_sim = 10000) 
dccs_deriv_sig <- clip_on_siggratia(dccs_deriv)
dccs_deriv_sig$sig # 1s and 0s are not continuous
dccs_deriv_plot_vals <- with(rle(dccs_deriv_sig$sig), 
                             mapply(function(s,e) range(dccs_deriv_sig$Age[s:e]),
                                    cumsum(lengths)[values == 1] - lengths[values == 1] + 1,
                                    cumsum(lengths)[values == 1]))
dccs_deriv_plot_range <- c(dccs_deriv_plot_vals[1,1], dccs_deriv_plot_vals[2,1], dccs_deriv_plot_vals[1,2], dccs_deriv_plot_vals[2,2])

# Plot dccs_age_gam effects
plot_model(dccs_age_gam, type = "pred", terms = "Age") # DCCS scores increase with age during childhood and early adolescence, but then plateau into late adolesence and adulthood
dccs_age_gam_eff <- ggpredict(dccs_age_gam, terms = c("Age"))
ggplot(dccs_age_gam_eff, aes(x, predicted)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  scale_y_continuous(breaks = seq(-1.5, 1.5, by = .75)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = .75), colour = 'white') +
  annotate("rect", xmin = dccs_deriv_plot_range[1], xmax = dccs_deriv_plot_range[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = dccs_blue) +
  annotate("rect", xmin = dccs_deriv_plot_range[3], xmax = dccs_deriv_plot_range[4], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = dccs_blue) +
  geom_line(size = 1.5, colour = dccs_blue) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .5, fill = dccs_blue) +
  labs(x = "Age (Years)", y = "Fitted DCCS Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_age_dccs.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - DCCS Contrasts ##### 

# Model uncorrected toolbox scores across age - DCCS compared to every other task
hist(no_dccs_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "dccs_gam.rda"))) {
  dccs_gam <- mgcv::gam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re"), data = no_dccs_long, method = "REML") 
  save(dccs_gam, file = paste0(analyzed_data_dir, "dccs_gam.rda"))
} else {
  load(paste0(analyzed_data_dir, "dccs_gam.rda"))
}
gam.check(dccs_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(dccs_gam) # Residuals should be normally distributed (yes)
summary(dccs_gam) # H0 = horizontal line is REJECTED for both levels of Measure; adjusted R2 is fine

# Plot dccs_gam effects
plot_model(dccs_gam, type = "pred", terms = "Age") # FC scores increase with age, especially during childhood
plot_model(dccs_gam, type = "pred", terms = "Measure") # DCCS scores are higher than non-DCCS scores
plot_model(dccs_gam, type = "pred", terms = c("Age", "Measure"), show.data = TRUE) # DCCS scores increase more dramatically with age than non-DCCS scores
dccs_gam_eff <- ggpredict(dccs_gam, terms = c("Age", "Measure"))
ggplot(dccs_gam_eff, aes(x, predicted)) + 
  scale_y_continuous(breaks = seq(-2, 2, by = 1)) +
  geom_hline(yintercept = seq(-2, 2, by = 1), colour = 'white') +
  geom_line(aes(color = group), size = 1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .5) +
  scale_color_manual(values = c(dccs_blue, comp_yellow), labels = no_dccs_labels) +
  scale_fill_manual(values = c(dccs_blue, comp_yellow), labels = no_dccs_labels) +
  labs(x = "Age (Years)", y = "Fitted Scores", color = "Measure", fill = "Measure") +
  plot_theme + theme(legend.position = "bottom")
ggsave(paste0(analyzed_data_dir, 'plot_dccs_eff.png'), width = 7, height = 7)

# Plot significant region of dccs_gam
if (!file.exists(paste0(analyzed_data_dir, "dccs_gam_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  dccs_gam_plot <- plot_diff(dccs_gam, view = "Age", comp = list(Measure = c("DCCS", "NoDCCS")),
                             sim.ci = TRUE) # DCCS and non-DCCS scores differ between ages 8.22 and 10.89, and again between 16.23 and 17.21
  save(dccs_gam_plot, file = paste0(analyzed_data_dir, "dccs_gam_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "dccs_gam_plot.rda"))
}
ggplot(dccs_gam_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 8.22, xmax = 10.89, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  annotate("rect", xmin = 16.23, xmax = 17.21, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = dccs_blue, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = dccs_blue) +
  labs(x = "Age (Years)", y = "Estimated Difference in\nDCCS - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_dccs.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - DCCS Contrasts with Random Slope ##### 

# Model uncorrected toolbox scores across age *with random slope* - DCCS compared to every other task
hist(no_dccs_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "dccs_gam_slopes.rda"))) {
  dccs_gam_slopes <- mgcv::bam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re") + s(Subject, Measure, bs = "re"),
                               data = no_dccs_long,
                               method = "fREML", # bam uses 'fREML' instead of 'REML'
                               discrete = TRUE, # Speeds up computation
                               nthreads = 4, # Use multiple cores if available
                               verbose = TRUE)
  save(dccs_gam_slopes, file = paste0(analyzed_data_dir, "dccs_gam_slopes.rda"))
} else {
  load(paste0(analyzed_data_dir, "dccs_gam_slopes.rda"))
}
gam.check(dccs_gam_slopes) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(dccs_gam_slopes) # Residuals should be normally distributed (yes)
mgcv::gam.vcomp(dccs_gam_slopes) # Examine random effects directly

# Plot significant region of dccs_gam_slopes
if (!file.exists(paste0(analyzed_data_dir, "dccs_gam_slopes_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  dccs_gam_slopes_plot <- plot_diff(dccs_gam_slopes, view = "Age", comp = list(Measure = c("DCCS", "NoDCCS")),
                                    sim.ci = TRUE) # DCCS and non-DCCS scores differ between ages 8.22 and 10.96, and again between 16.23 and 17.28
  save(dccs_gam_slopes_plot, file = paste0(analyzed_data_dir, "dccs_gam_slopes_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "dccs_gam_slopes_plot.rda"))
}
ggplot(dccs_gam_slopes_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 8.22, xmax = 10.96, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  annotate("rect", xmin = 16.23, xmax = 17.28, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +  
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = dccs_blue, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = dccs_blue) +
  labs(x = "Age (Years)", y = "Estimated Difference in\nDCCS - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_dccs_slopes.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - Pattern Comparison ##### 

# Model Pattern Comparison scores across age
hist(fluid_measures_uncorr$PatternComparison) 
pattern_age_gam <- mgcv::gam(PatternComparison ~ s(Age), data = fluid_measures_uncorr, method = "REML") 
gam.check(pattern_age_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(pattern_age_gam) # Residuals should be normally distributed (yes)
summary(pattern_age_gam) # H0 = horizontal line is REJECTED; adjusted R2 is fine

# Identify the region of significant age-related change by computing the first derivative
pattern_deriv <- derivatives(pattern_age_gam, n = 100, eps = 1e-07, level = .95, interval = "simultaneous", unconditional = TRUE, n_sim = 10000) 
pattern_deriv_sig <- clip_on_siggratia(pattern_deriv)
pattern_deriv_sig$sig # 1s and 0s are continuous
pattern_deriv_plot_range <- range(pattern_deriv_sig[which(pattern_deriv_sig$sig == 1), "Age"])

# Plot pattern_age_gam effects
plot_model(pattern_age_gam, type = "pred", terms = "Age") # Pattern Comparison scores increase with age during childhood and adolescence, but then plateau into adulthood
pattern_age_gam_eff <- ggpredict(pattern_age_gam, terms = c("Age"))
ggplot(pattern_age_gam_eff, aes(x, predicted)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  scale_y_continuous(breaks = seq(-1.5, 1.5, by = .75)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = .75), colour = 'white') +
  annotate("rect", xmin = pattern_deriv_plot_range[1], xmax = pattern_deriv_plot_range[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = pattern_blue) +
  geom_line(size = 1.5, colour = pattern_blue) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .5, fill = pattern_blue) +
  labs(x = "Age (Years)", y = "Fitted Pattern Comparison Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_age_pattern.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - Pattern Comparison Contrasts ##### 

# Model uncorrected toolbox scores across age - Pattern Comparison compared to every other task
hist(no_pattern_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "pattern_gam.rda"))) {
  pattern_gam <- mgcv::gam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re"), data = no_pattern_long, method = "REML") 
  save(pattern_gam, file = paste0(analyzed_data_dir, "pattern_gam.rda"))
} else {
  load(paste0(analyzed_data_dir, "pattern_gam.rda"))
}
gam.check(pattern_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(pattern_gam) # Residuals should be normally distributed (yes)
summary(pattern_gam) # H0 = horizontal line is REJECTED for both levels of Measure; adjusted R2 is fine

# Plot pattern_gam effects
plot_model(pattern_gam, type = "pred", terms = "Age") # FC scores increase with age, especially during childhood
plot_model(pattern_gam, type = "pred", terms = "Measure") # Non-Pattern scores are higher than Pattern scores
plot_model(pattern_gam, type = "pred", terms = c("Age", "Measure"), show.data = TRUE) # Pattern scores increase more dramatically with age than non-Pattern scores
pattern_gam_eff <- ggpredict(pattern_gam, terms = c("Age", "Measure"))
pattern_gam_eff$group <- as.character(pattern_gam_eff$group)
pattern_gam_eff$group <- replace(pattern_gam_eff$group, pattern_gam_eff$group == "PatternComparison", "Pattern")
pattern_gam_eff$group <- factor(pattern_gam_eff$group, levels = c("Pattern", "NoPattern"))
ggplot(pattern_gam_eff, aes(x, predicted)) + 
  scale_y_continuous(breaks = seq(-2, 2, by = 1)) +
  geom_hline(yintercept = seq(-2, 2, by = 1), colour = 'white') +
  geom_line(aes(color = group), size = 1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .5) +
  scale_color_manual(values = c(pattern_blue, comp_yellow), labels = no_pattern_labels) +
  scale_fill_manual(values = c(pattern_blue, comp_yellow), labels = no_pattern_labels) +
  labs(x = "Age (Years)", y = "Fitted Scores", color = "Measure", fill = "Measure") +
  plot_theme + theme(legend.position = "bottom")
ggsave(paste0(analyzed_data_dir, 'plot_pattern_eff.png'), width = 7, height = 7)

# Plot significant region of pattern_gam
if (!file.exists(paste0(analyzed_data_dir, "pattern_gam_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  pattern_gam_plot <- plot_diff(pattern_gam, view = "Age", comp = list(Measure = c("PatternComparison", "NoPattern")),
                                sim.ci = TRUE) # Pattern and non-Pattern scores differ between ages 8.22 and 12.65, and again between 16.86 and 21.99
  save(pattern_gam_plot, file = paste0(analyzed_data_dir, "pattern_gam_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "pattern_gam_plot.rda"))
}
ggplot(pattern_gam_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 8.22, xmax = 12.65, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  annotate("rect", xmin = 16.86, xmax = 21.99, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = pattern_blue, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = pattern_blue) +
  labs(x = "Age (Years)", y = "Estimated Difference in Pattern\nComparison - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_pattern.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - Pattern Comparison Contrasts with Random Slope ##### 

# Model uncorrected toolbox scores across age *with random slope* - Pattern Comparison compared to every other task
hist(no_pattern_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "pattern_gam_slopes.rda"))) {
  pattern_gam_slopes <- mgcv::bam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re") + s(Subject, Measure, bs = "re"),
                                  data = no_pattern_long,
                                  method = "fREML", # bam uses 'fREML' instead of 'REML'
                                  discrete = TRUE, # Speeds up computation
                                  nthreads = 4, # Use multiple cores if available
                                  verbose = TRUE)
  save(pattern_gam_slopes, file = paste0(analyzed_data_dir, "pattern_gam_slopes.rda"))
} else {
  load(paste0(analyzed_data_dir, "pattern_gam_slopes.rda"))
}
gam.check(pattern_gam_slopes) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(pattern_gam_slopes) # Residuals should be normally distributed (yes)
mgcv::gam.vcomp(pattern_gam_slopes) # Examine random effects directly

# Plot significant region of pattern_gam_slopes
if (!file.exists(paste0(analyzed_data_dir, "pattern_gam_slopes_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  pattern_gam_slopes_plot <- plot_diff(pattern_gam_slopes, view = "Age", comp = list(Measure = c("PatternComparison", "NoPattern")),
                                       sim.ci = TRUE) # Pattern and non-Pattern scores differ between ages 8.15 and 12.79, and again between 16.58 and 21.99
  save(pattern_gam_slopes_plot, file = paste0(analyzed_data_dir, "pattern_gam_slopes_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "pattern_gam_slopes_plot.rda"))
}
ggplot(pattern_gam_slopes_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 8.15, xmax = 12.79, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  annotate("rect", xmin = 16.58, xmax = 21.99, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +  
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = pattern_blue, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = pattern_blue) +
  labs(x = "Age (Years)", y = "Estimated Difference in Pattern\nComparison - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_pattern_slopes.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - PSMT ##### 

# Model PSMT scores across age
hist(fluid_measures_uncorr$PSMT) 
psmt_age_gam <- mgcv::gam(PSMT ~ s(Age), data = fluid_measures_uncorr, method = "REML") 
gam.check(psmt_age_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(psmt_age_gam) # Residuals should be normally distributed (yes)
summary(psmt_age_gam) # H0 = horizontal line is REJECTED; adjusted R2 is fine

# Identify the region of significant age-related change by computing the first derivative
psmt_deriv <- derivatives(psmt_age_gam, n = 100, eps = 1e-07, level = .95, interval = "simultaneous", unconditional = TRUE, n_sim = 10000) 
psmt_deriv_sig <- clip_on_siggratia(psmt_deriv)
psmt_deriv_sig$sig # 1s and 0s are not continuous 
psmt_deriv_plot_vals <- with(rle(psmt_deriv_sig$sig), 
                             mapply(function(s,e) range(psmt_deriv_sig$Age[s:e]),
                                    cumsum(lengths)[values == 1] - lengths[values == 1] + 1,
                                    cumsum(lengths)[values == 1]))
psmt_deriv_plot_range <- c(psmt_deriv_plot_vals[1,1], psmt_deriv_plot_vals[2,1], psmt_deriv_plot_vals[1,2], psmt_deriv_plot_vals[2,2])

# Plot psmt_age_gam effects
plot_model(psmt_age_gam, type = "pred", terms = "Age") # PSMT scores increase with age during childhood and adolescence, but then plateau into adulthood
psmt_age_gam_eff <- ggpredict(psmt_age_gam, terms = c("Age"))
ggplot(psmt_age_gam_eff, aes(x, predicted)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  scale_y_continuous(breaks = seq(-1.5, 1.5, by = .75)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = .75), colour = 'white') +
  annotate("rect", xmin = psmt_deriv_plot_range[1], xmax = psmt_deriv_plot_range[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = psmt_green) +
  annotate("rect", xmin = psmt_deriv_plot_range[3], xmax = psmt_deriv_plot_range[4], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = psmt_green) +
  geom_line(size = 1.5, colour = psmt_green) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .5, fill = psmt_green) +
  labs(x = "Age (Years)", y = "Fitted PSMT Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_age_psmt.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - PSMT Contrasts ##### 

# Model uncorrected toolbox scores across age - PSMT compared to every other task
hist(no_psmt_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "psmt_gam.rda"))) {
  psmt_gam <- mgcv::gam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re"), data = no_psmt_long, method = "REML") 
  save(psmt_gam, file = paste0(analyzed_data_dir, "psmt_gam.rda"))
} else {
  load(paste0(analyzed_data_dir, "psmt_gam.rda"))
}
gam.check(psmt_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(psmt_gam) # Residuals should be normally distributed (yes)
summary(psmt_gam) # H0 = horizontal line is REJECTED for both levels of Measure; adjusted R2 is fine

# Plot psmt_gam effects 
plot_model(psmt_gam, type = "pred", terms = "Age") # FC scores increase with age, especially during childhood
plot_model(psmt_gam, type = "pred", terms = "Measure") # Non-PSMT scores are higher than PSMT scores
plot_model(psmt_gam, type = "pred", terms = c("Age", "Measure"), show.data = TRUE) # Non-PSMT scores increase more dramatically with age than PSMT scores
psmt_gam_eff <- ggpredict(psmt_gam, terms = c("Age", "Measure"))
ggplot(psmt_gam_eff, aes(x, predicted)) + 
  scale_y_continuous(breaks = seq(-2, 2, by = 1)) +
  geom_hline(yintercept = seq(-2, 2, by = 1), colour = 'white') +
  geom_line(aes(color = group), size = 1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .5) +
  scale_color_manual(values = c(psmt_green, comp_yellow), labels = no_psmt_labels) +
  scale_fill_manual(values = c(psmt_green, comp_yellow), labels = no_psmt_labels) +
  labs(x = "Age (Years)", y = "Fitted Scores", color = "Measure", fill = "Measure") +
  plot_theme + theme(legend.position = "bottom")
ggsave(paste0(analyzed_data_dir, 'plot_psmt_eff.png'), width = 7, height = 7)

# Plot significant region of psmt_gam 
if (!file.exists(paste0(analyzed_data_dir, "psmt_gam_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  psmt_gam_plot <- plot_diff(psmt_gam, view = "Age", comp = list(Measure = c("PSMT", "NoPSMT")),
                             sim.ci = TRUE) # PSMT and non-PSMT scores differ between ages 8.01 and 9.70
  save(psmt_gam_plot, file = paste0(analyzed_data_dir, "psmt_gam_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "psmt_gam_plot.rda"))
}
ggplot(psmt_gam_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 8.01, xmax = 9.70, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = psmt_green, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = psmt_green) +
  labs(x = "Age (Years)", y = "Estimated Difference in\nPSMT - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_psmt.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - PSMT Contrasts with Random Slope ##### 

# Model uncorrected toolbox scores across age *with random slope* - PSMT compared to every other task
hist(no_psmt_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "psmt_gam_slopes.rda"))) {
  psmt_gam_slopes <- mgcv::bam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re") + s(Subject, Measure, bs = "re"),
                               data = no_psmt_long,
                               method = "fREML", # bam uses 'fREML' instead of 'REML'
                               discrete = TRUE, # Speeds up computation
                               nthreads = 4, # Use multiple cores if available
                               verbose = TRUE)
  save(psmt_gam_slopes, file = paste0(analyzed_data_dir, "psmt_gam_slopes.rda"))
} else {
  load(paste0(analyzed_data_dir, "psmt_gam_slopes.rda"))
}
gam.check(psmt_gam_slopes) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(psmt_gam_slopes) # Residuals should be normally distributed (yes)
mgcv::gam.vcomp(psmt_gam_slopes) # Examine random effects directly

# Plot significant region of psmt_gam_slopes
if (!file.exists(paste0(analyzed_data_dir, "psmt_gam_slopes_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  psmt_gam_slopes_plot <- plot_diff(psmt_gam_slopes, view = "Age", comp = list(Measure = c("PSMT", "NoPSMT")),
                                    sim.ci = TRUE) # PSMT and non-PSMT scores differ between ages 8.01 and 9.77
  save(psmt_gam_slopes_plot, file = paste0(analyzed_data_dir, "psmt_gam_slopes_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "psmt_gam_slopes_plot.rda"))
}
ggplot(psmt_gam_slopes_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 8.01, xmax = 9.77, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +  
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = psmt_green, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = psmt_green) +
  labs(x = "Age (Years)", y = "Estimated Difference in\nPSMT - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_psmt_slopes.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - List Sorting ##### 

# Model List Sorting scores across age
hist(fluid_measures_uncorr$ListSorting) 
list_age_gam <- mgcv::gam(ListSorting ~ s(Age), data = fluid_measures_uncorr, method = "REML") 
gam.check(list_age_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(list_age_gam) # Residuals should be normally distributed (yes)
summary(list_age_gam) # H0 = horizontal line is REJECTED; adjusted R2 is fine

# Identify the region of significant age-related change by computing the first derivative
list_deriv <- derivatives(list_age_gam, n = 100, eps = 1e-07, level = .95, interval = "simultaneous", unconditional = TRUE, n_sim = 10000) 
list_deriv_sig <- clip_on_siggratia(list_deriv)
list_deriv_sig$sig # 1s and 0s are continuous
list_deriv_plot_range <- range(list_deriv_sig[which(list_deriv_sig$sig == 1), "Age"])

# Plot list_age_gam effects
plot_model(list_age_gam, type = "pred", terms = "Age") # List Sorting scores increase with age during childhood and adolescence, but then plateau into adulthood
list_age_gam_eff <- ggpredict(list_age_gam, terms = c("Age"))
ggplot(list_age_gam_eff, aes(x, predicted)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  scale_y_continuous(breaks = seq(-1.5, 1.5, by = .75)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = .75), colour = 'white') +
  annotate("rect", xmin = list_deriv_plot_range[1], xmax = list_deriv_plot_range[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = list_green) +
  geom_line(size = 1.5, colour = list_green) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .5, fill = list_green) +
  labs(x = "Age (Years)", y = "Fitted List Sorting Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_age_list.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - List Sorting Contrasts ##### 

# Model uncorrected toolbox scores across age - List Sorting compared to every other task
hist(no_list_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "list_gam.rda"))) {
  list_gam <- mgcv::gam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re"), data = no_list_long, method = "REML") 
  save(list_gam, file = paste0(analyzed_data_dir, "list_gam.rda"))
} else {
  load(paste0(analyzed_data_dir, "list_gam.rda"))
}
gam.check(list_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(list_gam) # Residuals should be normally distributed (yes)
summary(list_gam) # H0 = horizontal line is REJECTED for both levels of Measure; adjusted R2 is fine

# Plot list_gam effects
plot_model(list_gam, type = "pred", terms = "Age") # FC scores increase with age, especially during childhood
plot_model(list_gam, type = "pred", terms = "Measure") # List scores are higher than non-List scores
plot_model(list_gam, type = "pred", terms = c("Age", "Measure"), show.data = TRUE) # Non-List scores increase more dramatically with age than List scores
list_gam_eff <- ggpredict(list_gam, terms = c("Age", "Measure"))
list_gam_eff$group <- as.character(list_gam_eff$group)
list_gam_eff$group <- replace(list_gam_eff$group, list_gam_eff$group == "ListSorting", "List")
list_gam_eff$group <- factor(list_gam_eff$group, levels = c("List", "NoList"))
ggplot(list_gam_eff, aes(x, predicted)) + 
  scale_y_continuous(breaks = seq(-2, 2, by = 1)) +
  geom_hline(yintercept = seq(-2, 2, by = 1), colour = 'white') +
  geom_line(aes(color = group), size = 1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .5) +
  scale_color_manual(values = c(list_green, comp_yellow), labels = no_list_labels) +
  scale_fill_manual(values = c(list_green, comp_yellow), labels = no_list_labels) +
  labs(x = "Age (Years)", y = "Fitted Scores", color = "Measure", fill = "Measure") +
  plot_theme + theme(legend.position = "bottom")
ggsave(paste0(analyzed_data_dir, 'plot_list_eff.png'), width = 7, height = 7)

# Plot significant region of list_gam
if (!file.exists(paste0(analyzed_data_dir, "list_gam_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  list_gam_plot <- plot_diff(list_gam, view = "Age", comp = list(Measure = c("ListSorting", "NoList")),
                             sim.ci = TRUE) # List and non-List scores differ between ages 19.32 and 21.29
  save(list_gam_plot, file = paste0(analyzed_data_dir, "list_gam_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "list_gam_plot.rda"))
}
ggplot(list_gam_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 19.32, xmax = 21.29, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = list_green, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = list_green) +
  labs(x = "Age (Years)", y = "Estimated Difference in List\nSorting - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_list.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - List Sorting Contrasts with Random Slope ##### 

# Model uncorrected toolbox scores across age *with random slope* - List Sorting compared to every other task
hist(no_list_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "list_gam_slopes.rda"))) {
  list_gam_slopes <- mgcv::bam(Score ~ s(Age, by = Measure) + Measure + s(Subject, bs = "re") + s(Subject, Measure, bs = "re"),
                               data = no_list_long,
                               method = "fREML", # bam uses 'fREML' instead of 'REML'
                               discrete = TRUE, # Speeds up computation
                               nthreads = 4, # Use multiple cores if available
                               verbose = TRUE)
  save(list_gam_slopes, file = paste0(analyzed_data_dir, "list_gam_slopes.rda"))
} else {
  load(paste0(analyzed_data_dir, "list_gam_slopes.rda"))
}
gam.check(list_gam_slopes) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(list_gam_slopes) # Residuals should be normally distributed (yes)
mgcv::gam.vcomp(list_gam_slopes) # Examine random effects directly

# Plot significant region of list_gam_slopes
if (!file.exists(paste0(analyzed_data_dir, "list_gam_slopes_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  list_gam_slopes_plot <- plot_diff(list_gam_slopes, view = "Age", comp = list(Measure = c("ListSorting", "NoList")),
                                    sim.ci = TRUE) # List and non-List scores differ between ages 19.18 and 21.29
  save(list_gam_slopes_plot, file = paste0(analyzed_data_dir, "list_gam_slopes_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "list_gam_slopes_plot.rda"))
}
ggplot(list_gam_slopes_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 19.18, xmax = 21.29, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.75, -.375, .375, .75), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.75, .75, by = .375)) +  
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = list_green, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = list_green) +
  labs(x = "Age (Years)", y = "Estimated Difference in List\nSorting - Remaining Task Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_list_slopes.png'), width = 7, height = 7)

##### Examine NIH Toolbox Data Across Age - Save Statistical Results ##### 

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'uncorr_toolbox_age_gams.txt'))

cat("Flanker\n")
cat("-------\n")
cat(gam.check(flanker_age_gam))
cat("\n")
summary(flanker_age_gam)

cat("\n\n")
cat("DCCS\n")
cat("----\n")
cat(gam.check(dccs_age_gam))
cat("\n")
summary(dccs_age_gam)

cat("\n\n")
cat("Pattern Comparison\n")
cat("------------------\n")
cat(gam.check(pattern_age_gam))
cat("\n")
summary(pattern_age_gam)

cat("\n\n")
cat("PSMT\n")
cat("----\n")
cat(gam.check(psmt_age_gam))
cat("\n")
summary(psmt_age_gam)

cat("\n\n")
cat("List Sorting\n")
cat("------------\n")
cat(gam.check(list_age_gam))
cat("\n")
summary(list_age_gam)

# Stop writing to file
sink()

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'uncorr_toolbox_gams.txt'))

cat("Flanker\n")
cat("-------\n")
cat(gam.check(flanker_gam))
cat("\n")
summary(flanker_gam)

cat("\n\n")
cat("DCCS\n")
cat("----\n")
cat(gam.check(dccs_gam))
cat("\n")
summary(dccs_gam)

cat("\n\n")
cat("Pattern Comparison\n")
cat("------------------\n")
cat(gam.check(pattern_gam))
cat("\n")
summary(pattern_gam)

cat("\n\n")
cat("PSMT\n")
cat("----\n")
cat(gam.check(psmt_gam))
cat("\n")
summary(psmt_gam)

cat("\n\n")
cat("List Sorting\n")
cat("------------\n")
cat(gam.check(list_gam))
cat("\n")
summary(list_gam)

# Stop writing to file
sink()

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'uncorr_toolbox_sig_ages.txt'))

cat("Flanker\n")
cat("-------\n")
cat(paste("Significant region of age-related change:", flanker_deriv_plot_range[1], flanker_deriv_plot_range[2], "\n"))
set.seed(123)
plot_diff(flanker_gam, view = "Age", comp = list(Measure = c("Flanker", "NoFlanker")), sim.ci = TRUE)

cat("\n\n")
cat("DCCS\n")
cat("----\n")
cat(paste("Significant region of age-related change:", dccs_deriv_plot_range[1], dccs_deriv_plot_range[2], dccs_deriv_plot_range[3], dccs_deriv_plot_range[4], "\n"))
set.seed(123)
plot_diff(dccs_gam, view = "Age", comp = list(Measure = c("DCCS", "NoDCCS")), sim.ci = TRUE)

cat("\n\n")
cat("Pattern Comparison\n")
cat("------------------\n")
cat(paste("Significant region of age-related change:", pattern_deriv_plot_range[1], pattern_deriv_plot_range[2], "\n"))
set.seed(123)
plot_diff(pattern_gam, view = "Age", comp = list(Measure = c("PatternComparison", "NoPattern")), sim.ci = TRUE)

cat("\n\n")
cat("PSMT\n")
cat("----\n")
cat(paste("Significant region of age-related change:", psmt_deriv_plot_range[1], psmt_deriv_plot_range[2], psmt_deriv_plot_range[3], psmt_deriv_plot_range[4], "\n"))
set.seed(123)
plot_diff(psmt_gam, view = "Age", comp = list(Measure = c("PSMT", "NoPSMT")), sim.ci = TRUE)

cat("\n\n")
cat("List Sorting\n")
cat("------------\n")
cat(paste("Significant region of age-related change:", list_deriv_plot_range[1], list_deriv_plot_range[2], "\n"))
set.seed(123)
plot_diff(list_gam, view = "Age", comp = list(Measure = c("ListSorting", "NoList")), sim.ci = TRUE)

# Stop writing to file
sink()
closeAllConnections()

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'uncorr_toolbox_gams_slopes.txt'))

cat("Flanker\n")
cat("-------\n")
mgcv::gam.vcomp(flanker_gam_slopes)

cat("\n\n")
cat("DCCS\n")
cat("----\n")
mgcv::gam.vcomp(dccs_gam_slopes)

cat("\n\n")
cat("Pattern Comparison\n")
cat("------------------\n")
mgcv::gam.vcomp(pattern_gam_slopes)

cat("\n\n")
cat("PSMT\n")
cat("----\n")
mgcv::gam.vcomp(psmt_gam_slopes)

cat("\n\n")
cat("List Sorting\n")
cat("------------\n")
mgcv::gam.vcomp(list_gam_slopes)

# Stop writing to file
sink()
closeAllConnections()

##### Examine Cognitive Control Factor Scores Across Age ##### 

# Model factor scores across age
hist(cfa_2_scores$CgCn) 
cgcn_age_gam <- mgcv::gam(CgCn ~ s(Age), data = cfa_2_scores, method = "REML") 
gam.check(cgcn_age_gam) # k-index should be larger than 1 (no), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(cgcn_age_gam) # Residuals should be normally distributed (yes)
summary(cgcn_age_gam) # H0 = horizontal line is REJECTED; adjusted R2 is fine

# Identify the region of significant age-related change by computing the first derivative
cgcn_deriv <- derivatives(cgcn_age_gam, n = 100, eps = 1e-07, level = .95, interval = "simultaneous", unconditional = TRUE, n_sim = 10000) 
cgcn_deriv_sig <- clip_on_siggratia(cgcn_deriv)
cgcn_deriv_sig$sig # 1s and 0s are continuous
cgcn_deriv_plot_range <- range(cgcn_deriv_sig[which(cgcn_deriv_sig$sig == 1), "Age"])

# Plot cgcn_age_gam effects
plot_model(cgcn_age_gam, type = "pred", terms = "Age") # CgCn scores increase with age during childhood and adolescence, but then plateau into adulthood
cgcn_age_gam_eff <- ggpredict(cgcn_age_gam, terms = c("Age"))
ggplot(cgcn_age_gam_eff, aes(x, predicted)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  scale_y_continuous(breaks = seq(-1, 1, by = .5)) +
  geom_hline(yintercept = seq(-1, 1, by = .5), colour = 'white') +
  annotate("rect", xmin = cgcn_deriv_plot_range[1], xmax = cgcn_deriv_plot_range[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = dccs_blue) +
  geom_line(size = 1.5, colour = flanker_blue) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .5, fill = flanker_blue) +
  labs(x = "Age (Years)", y = "Cognitive Control Factor Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_age_cgcn.png'), width = 7, height = 7)

##### Examine Memory Factor Scores Across Age ##### 

# Model factor scores across age
hist(cfa_2_scores$Mmry) 
mmry_age_gam <- mgcv::gam(Mmry ~ s(Age), data = cfa_2_scores, method = "REML") 
gam.check(mmry_age_gam) # k-index should be larger than 1 (no), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(mmry_age_gam) # Residuals should be normally distributed (yes)
summary(mmry_age_gam) # H0 = horizontal line is REJECTED; adjusted R2 is fine

# Identify the region of significant age-related change by computing the first derivative
mmry_deriv <- derivatives(mmry_age_gam, n = 100, eps = 1e-07, level = .95, interval = "simultaneous", unconditional = TRUE, n_sim = 10000) 
mmry_deriv_sig <- clip_on_siggratia(mmry_deriv)
mmry_deriv_sig$sig # 1s and 0s are continuous
mmry_deriv_plot_range <- range(mmry_deriv_sig[which(mmry_deriv_sig$sig == 1), "Age"])

# Plot mmry_age_gam effects
plot_model(mmry_age_gam, type = "pred", terms = "Age") # Mmry scores increase with age during childhood and adolescence, but then plateau into adulthood
mmry_age_gam_eff <- ggpredict(mmry_age_gam, terms = c("Age"))
ggplot(mmry_age_gam_eff, aes(x, predicted)) + 
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  scale_y_continuous(breaks = seq(-1, 1, by = .5)) +
  geom_hline(yintercept = seq(-1, 1, by = .5), colour = 'white') +
  annotate("rect", xmin = mmry_deriv_plot_range[1], xmax = mmry_deriv_plot_range[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = list_green) +
  geom_line(size = 1.5, colour = psmt_green) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .5, fill = psmt_green) +
  labs(x = "Age (Years)", y = "Memory Factor Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_age_mmry.png'), width = 7, height = 7)

##### Examine Factor Scores Across Age - Contrasts ##### 

# Convert factor score data frame to long format
factor_long <- cfa_2_scores %>% pivot_longer(cols = c("CgCn", "Mmry"), names_to = "Factor", "values_to" = "Score")
factor_long$Factor <- factor(factor_long$Factor, levels = c("CgCn", "Mmry"))
factor_long$Subject <- as.factor(factor_long$Subject)

# Plot factor scores across age - CgCn compared to Mmry
ggplot(data = factor_long, aes(x = Age, y = Score, color = Factor, fill = Factor)) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-3, -2, -1, 1, 2, 3), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-3, 3, by = 1)) +   
  geom_point(alpha = .2, size = 2.5, stroke = NA) +
  geom_smooth(method = "gam", alpha = .5, size = 1.5, se = FALSE) + 
  scale_color_manual(values = c(flanker_blue, psmt_green), labels = c("Cognitive Control", "Memory")) +
  scale_fill_manual(values = c(flanker_blue, psmt_green), labels = c("Cognitive Control", "Memory")) +
  labs(x = "Age (Years)", y = "Factor Score") +
  plot_theme + theme(legend.position = "bottom")
ggsave(paste0(analyzed_data_dir, 'factor_contrast_age.png'), width = 9, height = 7)

# Plot factor scores across age - CgCn compared to Mmry *as difference score*
cfa_2_scores$Diff <- cfa_2_scores$CgCn - cfa_2_scores$Mmry
ggplot(data = cfa_2_scores, aes(x = Age, y = Diff)) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-3, -2, -1, 1, 2, 3), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-3, 3, by = 1)) +   
  geom_point(alpha = .2, size = 2.5, stroke = NA, color = comp_yellow) +
  geom_smooth(method = "gam", alpha = .5, size = 1.5, se = FALSE, color = comp_yellow) + 
  labs(x = "Age (Years)", y = "Difference in Cognitive Control -\nMemory Factor Scores") +
  plot_theme 
ggsave(paste0(analyzed_data_dir, 'factor_contrast_diff_age.png'), width = 7, height = 7)

# Model factor scores across age - CgCn compared to Mmry 
hist(factor_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "factor_gam.rda"))) {
  factor_gam <- mgcv::gam(Score ~ s(Age, by = Factor) + Factor + s(Subject, bs = "re"), data = factor_long, method = "REML") 
  save(factor_gam, file = paste0(analyzed_data_dir, "factor_gam.rda"))
} else {
  load(paste0(analyzed_data_dir, "factor_gam.rda"))
}
gam.check(factor_gam) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(factor_gam) # Residuals should be normally distributed (yes)
summary(factor_gam) # H0 = horizontal line is REJECTED for both levels of Factor; adjusted R2 is fine

# Plot factor_gam effects
plot_model(factor_gam, type = "pred", terms = "Age") # FC scores increase with age, especially during childhood
plot_model(factor_gam, type = "pred", terms = "Factor") # CgCn factor scores are higher than Mmry factor scores
plot_model(factor_gam, type = "pred", terms = c("Age", "Factor"), show.data = TRUE) # CgCn scores increase more dramatically with age than Mmry scores
factor_gam_eff <- ggpredict(factor_gam, terms = c("Age", "Factor"))
ggplot(factor_gam_eff, aes(x, predicted)) + 
  scale_y_continuous(breaks = seq(-1.5, 1.5, by = .75)) +
  geom_hline(yintercept = seq(-1.5, 1.5, by = .75), colour = 'white') +
  geom_line(aes(color = group), size = 1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .5) +
  scale_color_manual(values = c(flanker_blue, psmt_green), labels = c("Cognitive Control", "Memory")) +
  scale_fill_manual(values = c(flanker_blue, psmt_green), labels = c("Cognitive Control", "Memory")) +
  labs(x = "Age (Years)", y = "Fitted Factor Scores", color = "Factor", fill = "Factor") +
  plot_theme + theme(legend.position = "bottom")
ggsave(paste0(analyzed_data_dir, 'plot_factor_eff.png'), width = 7, height = 7)

# Plot significant region of factor_gam
if (!file.exists(paste0(analyzed_data_dir, "factor_gam_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  factor_gam_plot <- plot_diff(factor_gam, view = "Age", comp = list(Factor = c("CgCn", "Mmry")),
                               sim.ci = TRUE) # Factor scores differ between ages 8.01 and 11.80, and again between 14.54 and 21.99
  save(factor_gam_plot, file = paste0(analyzed_data_dir, "factor_gam_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "factor_gam_plot.rda"))
}
ggplot(factor_gam_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 8.01, xmax = 11.80, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  annotate("rect", xmin = 14.54, xmax = 21.99, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.5, -.25, .25, .5), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.5, .5, by = .25)) +  
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = comp_yellow, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = comp_yellow) +
  labs(x = "Age (Years)", y = "Estimated Difference in Cognitive\nControl - Memory Factor Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_factor.png'), width = 7, height = 7)

##### Examine Factor Scores Across Age - Contrasts with Random Slope ##### 

# Model factor scores across age *with random slope* - CgCn compared to Mmry 
hist(factor_long$Score) 
if (!file.exists(paste0(analyzed_data_dir, "factor_gam_slopes.rda"))) {
  factor_gam_slopes <- mgcv::bam(Score ~ s(Age, by = Factor) + Factor + s(Subject, bs = "re") + s(Subject, Factor, bs = "re"),
                                 data = factor_long,
                                 method = "fREML", # bam uses 'fREML' instead of 'REML'
                                 discrete = TRUE, # Speeds up computation
                                 nthreads = 4, # Use multiple cores if available
                                 verbose = TRUE)
  save(factor_gam_slopes, file = paste0(analyzed_data_dir, "factor_gam_slopes.rda"))
} else {
  load(paste0(analyzed_data_dir, "factor_gam_slopes.rda"))
}
gam.check(factor_gam_slopes) # k-index should be larger than 1 (yes), p-value should be large (fine), edf should be clearly smaller than k' (yes)
appraise(factor_gam_slopes) # Residuals should be normally distributed (yes)
mgcv::gam.vcomp(factor_gam_slopes) # Examine random effects directly

# Plot significant region of factor_gam_slopes
if (!file.exists(paste0(analyzed_data_dir, "factor_gam_slopes_plot.rda"))) {
  set.seed(123) # For sim.ci consistency
  factor_gam_slopes_plot <- plot_diff(factor_gam_slopes, view = "Age", comp = list(Factor = c("CgCn", "Mmry")),
                                      sim.ci = TRUE) # Factor scores differ between ages 8.01 and 11.80, and again between 14.54 and 21.99
  save(factor_gam_slopes_plot, file = paste0(analyzed_data_dir, "factor_gam_slopes_plot.rda"))
} else {
  load(paste0(analyzed_data_dir, "factor_gam_slopes_plot.rda"))
}
ggplot(factor_gam_slopes_plot, aes(x = Age, y = est)) + 
  annotate("rect", xmin = 8.01, xmax = 11.80, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  annotate("rect", xmin = 14.54, xmax = 21.99, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = comp_yellow) +
  scale_x_continuous(breaks = seq(8, 22, by = 2)) +
  geom_hline(yintercept = c(-.5, -.25, .25, .5), colour = 'grey80') +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  scale_y_continuous(breaks = seq(-.5, .5, by = .25)) +  
  geom_ribbon(aes(ymin = est - sim.CI, ymax = est + sim.CI), fill = comp_yellow, alpha = .5, show.legend = FALSE) + 
  geom_line(size = 1.5, colour = comp_yellow) +
  labs(x = "Age (Years)", y = "Estimated Difference in Cognitive\nControl - Memory Factor Scores") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'plot_diff_factor_slopes.png'), width = 7, height = 7)

##### Examine Factor Scores Across Age - Save Statistical Results ##### 

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'factor_gams.txt'))

cat("CgCn Age GAM\n")
cat("------------\n")
cat(gam.check(cgcn_age_gam))
cat("\n")
summary(cgcn_age_gam)

cat("\n\n")
cat("Mmry Age GAM\n")
cat("------------\n")
cat(gam.check(mmry_age_gam))
cat("\n")
summary(mmry_age_gam)

cat("\n\n")
cat("Factor GAM\n")
cat("----------\n")
cat(gam.check(factor_gam))
cat("\n")
summary(factor_gam)

cat("\n\n")
cat("Significant Ages\n")
cat("----------------\n")
cat(paste("Significant region of CgCn age-related change:", cgcn_deriv_plot_range[1], cgcn_deriv_plot_range[2], "\n"))
cat(paste("Significant region of Mmry age-related change:", mmry_deriv_plot_range[1], mmry_deriv_plot_range[2], "\n"))
set.seed(123)
plot_diff(factor_gam, view = "Age", comp = list(Factor = c("CgCn", "Mmry")), sim.ci = TRUE)

cat("\n\n")
cat("Random Effects\n")
cat("--------------\n")
mgcv::gam.vcomp(factor_gam_slopes)

# Stop writing to file
sink()
