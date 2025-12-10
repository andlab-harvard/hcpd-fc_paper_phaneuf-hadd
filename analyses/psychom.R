##### Conduct HCPD FC Psychometric Analyses #####

# Written by: Camille Phaneuf-Hadd (cphaneuf@g.harvard.edu)
# Last updated: 12/4/25
# Question: how do different measures of FC change across development?

##### Set up Script ##### 

# Load needed libraries
require(pacman) # for p_load()
p_load(tidyverse, # for df manipulation
       dplyr, # for %>% and other operators
       plotrix, # for std.error()
       ggplot2, # for static plotting
       sjPlot, # for plot_model()
       lmerTest, # for mixed effects models
       Rmisc, # for summarySE()
       ggcorrplot, # for ggcorrplot()
       smacof, # for MDS functions
       psych, # for fa() and cortest.jenrich()
       lavaan, # for cfa()
       semPlot, # for semPaths()
       semTools, # for multigroup CFA comparisons
       reshape2, # for melt()
       qgraph, # for qgraph()
       NetworkComparisonTest, # for NCT()
       WebPower) # for multigroup CFA power analysis

# Load shared HCPD FC functions
source("utilities.R")

# Set data directory paths
getwd()
processed_data_dir <- "../data/"
analyzed_data_dir <- "../results/psychom/"

# Read in demographic and behavioral data
uncorr_toolbox <- read.csv(paste0(processed_data_dir, "uncorr_toolbox.csv"))

##### Examine Distribution of NIH Toolbox Data ##### 

# Plot Flanker score distribution
range(uncorr_toolbox$Flanker_z)
ggplot(data = uncorr_toolbox, aes(x = Flanker_z)) +
  scale_x_continuous(breaks = seq(-5, 3, by = 1)) +
  geom_hline(yintercept = seq(0, 250, by = 50), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 250, by = 50)) + 
  geom_histogram(binwidth = .5, alpha = .75, fill = flanker_blue, colour = flanker_blue) +
  geom_vline(xintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x = "Z-Scored Flanker", y = "Number of Participants") +
  plot_theme
ggsave(paste0(analyzed_data_dir, "flanker_z_hist.png"), plot = last_plot(), width = 9, height = 7)

# Plot DCCS score distribution
range(uncorr_toolbox$DCCS_z)
ggplot(data = uncorr_toolbox, aes(x = DCCS_z)) +
  scale_x_continuous(breaks = seq(-5, 3, by = 1)) +
  geom_hline(yintercept = seq(0, 250, by = 50), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 250, by = 50)) + 
  geom_histogram(binwidth = .5, alpha = .75, fill = dccs_blue, colour = dccs_blue) +
  geom_vline(xintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x = "Z-Scored DCCS", y = "Number of Participants") +
  plot_theme
ggsave(paste0(analyzed_data_dir, "dccs_z_hist.png"), plot = last_plot(), width = 9, height = 7)

# Plot Pattern Comparison score distribution
range(uncorr_toolbox$Pattern_z)
ggplot(data = uncorr_toolbox, aes(x = Pattern_z)) +
  scale_x_continuous(breaks = seq(-5, 3, by = 1)) +
  geom_hline(yintercept = seq(0, 250, by = 50), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 250, by = 50)) + 
  geom_histogram(binwidth = .5, alpha = .75, fill = pattern_blue, colour = pattern_blue) +
  geom_vline(xintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x = "Z-Scored Pattern Comparison", y = "Number of Participants") +
  plot_theme
ggsave(paste0(analyzed_data_dir, "pattern_z_hist.png"), plot = last_plot(), width = 9, height = 7)

# Plot PSMT score distribution
range(uncorr_toolbox$PSMT_z)
ggplot(data = uncorr_toolbox, aes(x = PSMT_z)) +
  scale_x_continuous(breaks = seq(-5, 3, by = 1)) +
  geom_hline(yintercept = seq(0, 250, by = 50), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 250, by = 50)) + 
  geom_histogram(binwidth = .5, alpha = .75, fill = psmt_green, colour = psmt_green) +
  geom_vline(xintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x = "Z-Scored PSMT", y = "Number of Participants") +
  plot_theme
ggsave(paste0(analyzed_data_dir, "psmt_z_hist.png"), plot = last_plot(), width = 9, height = 7)

# Plot List Sorting score distribution
range(uncorr_toolbox$List_z)
ggplot(data = uncorr_toolbox, aes(x = List_z)) +
  scale_x_continuous(breaks = seq(-5, 3, by = 1)) +
  geom_hline(yintercept = seq(0, 250, by = 50), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 250, by = 50)) + 
  geom_histogram(binwidth = .5, alpha = .75, fill = list_green, colour = list_green) +
  geom_vline(xintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x = "Z-Scored List Sorting", y = "Number of Participants") +
  plot_theme
ggsave(paste0(analyzed_data_dir, "list_z_hist.png"), plot = last_plot(), width = 9, height = 7)

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'toolbox_means_sds.txt'))

cat("FLANKER\n")
cat("-------\n")
cat("Mean:", mean(uncorr_toolbox$Flanker), "\n")
cat("Standard Deviation:", sd(uncorr_toolbox$Flanker), "\n")

cat("\nDCCS\n")
cat("----\n")
cat("Mean:", mean(uncorr_toolbox$DCCS), "\n")
cat("Standard Deviation:", sd(uncorr_toolbox$DCCS), "\n")

cat("\nPATTERN COMPARISON\n")
cat("------------------\n")
cat("Mean:", mean(uncorr_toolbox$Pattern), "\n")
cat("Standard Deviation:", sd(uncorr_toolbox$Pattern), "\n")

cat("\nPSMT\n")
cat("----\n")
cat("Mean:", mean(uncorr_toolbox$PSMT), "\n")
cat("Standard Deviation:", sd(uncorr_toolbox$PSMT), "\n")

cat("\nLIST SORTING\n")
cat("------------\n")
cat("Mean:", mean(uncorr_toolbox$List), "\n")
cat("Standard Deviation:", sd(uncorr_toolbox$List), "\n")

# Stop writing to file
sink()

##### Examine Relations Among NIH Toolbox Data ##### 

# Create correlation matrix among toolbox scores only
toolbox_only <- uncorr_toolbox[, toolbox_z]
names(toolbox_only) <- short_task_labels
toolbox_cormat <- cor(toolbox_only) # Compute Pearson correlations
toolbox_cormat

# Confirm that toolbox_only already has complete cases, so use = "pairwise.complete.obs" did not need to be set in cor()
sum(complete.cases(toolbox_only)) == dim(toolbox_only)[1] # TRUE

# Plot correlations among uncorrected toolbox scores
# Resources used:
# - (general) https://rpkgs.datanovia.com/ggcorrplot/reference/ggcorrplot.html
# - (axis colors) https://community.rstudio.com/t/axis-labels-with-individual-colors/77848
ggcorrplot(toolbox_cormat, digits = 2, method = "circle", lab = TRUE, lab_size = 1.25, 
           lab_col = 'white', outline.color = 'black') + 
  scale_fill_gradient2(low = 'white', high = reg_gray, breaks = c(0, 1), limit = c(0, 1)) +
  labs(fill = NULL) +
  xlab(label = NULL) + ylab(label = NULL) +
  plot_theme + 
  theme(axis.text.x = element_text(colour = cols),
        axis.text.y = element_text(colour = cols))
ggsave(paste0(analyzed_data_dir, 'uncorr_toolbox_cormat.png'), width = 3.25, height = 3.25, bg = 'white')

##### Examine Relations Among NIH Toolbox Data by Age Group ##### 

# Get floored age and age groups
uncorr_toolbox$Age_Floor <- floor(uncorr_toolbox$Age)
uncorr_toolbox <- uncorr_toolbox %>% mutate(Age_Group = case_when((Age_Floor <= 10) ~ "8-10 Years",
                                                                  (Age_Floor >= 11 & Age_Floor <= 13) ~ "11-13 Years",
                                                                  (Age_Floor >= 14 & Age_Floor <= 17) ~ "14-17 Years",
                                                                  (Age_Floor >= 18) ~ "18-21 Years"))

# Create data frames for each age group
age_8_10 <- uncorr_toolbox[uncorr_toolbox$Age_Group == "8-10 Years", ]
age_11_13 <- uncorr_toolbox[uncorr_toolbox$Age_Group == "11-13 Years", ]
age_14_17 <- uncorr_toolbox[uncorr_toolbox$Age_Group == "14-17 Years", ]
age_18_21 <- uncorr_toolbox[uncorr_toolbox$Age_Group == "18-21 Years", ]
age_8_10 <- age_8_10[, toolbox_z]
age_11_13 <- age_11_13[, toolbox_z]
age_14_17 <- age_14_17[, toolbox_z]
age_18_21 <- age_18_21[, toolbox_z]
names(age_8_10) <- short_task_labels
names(age_11_13) <- short_task_labels
names(age_14_17) <- short_task_labels
names(age_18_21) <- short_task_labels

# Confirm that age_X_Y already has complete cases, so use = "pairwise.complete.obs" does not need to be set in cor()
sum(complete.cases(age_8_10)) == dim(age_8_10)[1] # TRUE
sum(complete.cases(age_11_13)) == dim(age_11_13)[1] # TRUE
sum(complete.cases(age_14_17)) == dim(age_14_17)[1] # TRUE
sum(complete.cases(age_18_21)) == dim(age_18_21)[1] # TRUE

# Create correlation network of Partial Correlations for each age group
#   NOTE: a partial correlation (e.g., between PSMT and List) implies partialling out 
#         the effects of all other variables in the network (i.e., Flanker, DCCS, and Pattern)
pearson <- qgraph(toolbox_cormat[1:5, 1:5], layout = "spring", graph = "cor")
png(paste0(analyzed_data_dir, 'age_partial_cor.png'), width = 900, height = 1400)
op <- par(mfrow = c(2, 2))
par(family = "Avenir", mgp = c(2.5, 1, 0))
partial_8_10 <- qgraph(cor(age_8_10), layout = pearson$layout, graph = "pcor", 
                       labels = FALSE, edge.labels = TRUE, sampleSize = nrow(toolbox_only), edge.color = age_8_10_orange)  
text(partial_8_10$layout[, 1], partial_8_10$layout[, 2], labels = short_task_labels, cex = node_fontsize / par("cex"))
title("8-10 Years", font.main = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1, 
      center = TRUE)
partial_11_13 <- qgraph(cor(age_11_13), layout = pearson$layout, graph = "pcor", 
                        labels = FALSE, edge.labels = TRUE, sampleSize = nrow(toolbox_only), edge.color = age_11_13_red)  
text(partial_11_13$layout[, 1], partial_11_13$layout[, 2], labels = short_task_labels, cex = node_fontsize / par("cex"))
title("11-13 Years", font.main = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1, 
      center = TRUE)
partial_14_17 <- qgraph(cor(age_14_17), layout = pearson$layout, graph = "pcor", 
                        labels = FALSE, edge.labels = TRUE, sampleSize = nrow(toolbox_only), edge.color = age_14_17_pink)  
text(partial_14_17$layout[, 1], partial_14_17$layout[, 2], labels = short_task_labels, cex = node_fontsize / par("cex"))
title("14-17 Years", font.main = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1, 
      center = TRUE)
partial_18_21 <- qgraph(cor(age_18_21), layout = pearson$layout, graph = "pcor", 
                        labels = FALSE, edge.labels = TRUE, sampleSize = nrow(toolbox_only), edge.color = age_18_21_purple)  
text(partial_18_21$layout[, 1], partial_18_21$layout[, 2], labels = short_task_labels, cex = node_fontsize / par("cex"))
title("18-21 Years", font.main = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1, 
      center = TRUE)
par(op)
dev.off()

# Inspect range of Partial Correlations to manually add legend in annotated_figs.pptx
# max: 0.40
# max: 0.37
# max: 0.41
# max: 0.44

# Perform network comparison test between each age group
set.seed(123)
nct_child_vs_earlyado <- NCT(age_8_10, age_11_13, test.edges = TRUE, it = 250, p.adjust.methods = "holm")
set.seed(123)
nct_child_vs_lateado <- NCT(age_8_10, age_14_17, test.edges = TRUE, it = 250, p.adjust.methods = "holm")
set.seed(123)
nct_child_vs_adult <- NCT(age_8_10, age_18_21, test.edges = TRUE, it = 250, p.adjust.methods = "holm")
set.seed(123)
nct_earlyado_vs_lateado <- NCT(age_11_13, age_14_17, test.edges = TRUE, it = 250, p.adjust.methods = "holm")
set.seed(123)
nct_earlyado_vs_adult <- NCT(age_11_13, age_18_21, test.edges = TRUE, it = 250, p.adjust.methods = "holm")
set.seed(123)
nct_lateado_vs_adult <- NCT(age_14_17, age_18_21, test.edges = TRUE, it = 250, p.adjust.methods = "holm")

# Save off network comparison tests
sink(paste0(analyzed_data_dir, 'network_comp_test.txt')) # Write to file instead of the terminal
cat("Children vs. Early Adolescence\n")
cat("------------------------------\n")
nct_child_vs_earlyado
cat("\n")
cat("Children vs. Late Adolescence\n")
cat("-----------------------------\n")
nct_child_vs_lateado
cat("\n")
cat("Children vs. Adults\n")
cat("-------------------\n")
nct_child_vs_adult
cat("\n")
cat("Early Adolescence vs. Late Adolescence\n")
cat("--------------------------------------\n")
nct_earlyado_vs_lateado
cat("\n")
cat("Early Adolescence vs. Adults\n")
cat("----------------------------\n")
nct_earlyado_vs_adult
cat("\n")
cat("Late Adolescence vs. Adults\n")
cat("---------------------------\n")
nct_lateado_vs_adult
sink() # Stop writing to file

# Perform Jenrich's test for equality of Pearson correlation matrices between each age group
#   NOTE: H0 = two correlation matrices are equal in the population
jen_child_vs_earlyado <- cortest.jennrich(cor(age_8_10), cor(age_11_13), n1 = nrow(age_8_10), n2 = nrow(age_11_13))
jen_child_vs_lateado <- cortest.jennrich(cor(age_8_10), cor(age_14_17), n1 = nrow(age_8_10), n2 = nrow(age_14_17))   
jen_child_vs_adult <- cortest.jennrich(cor(age_8_10), cor(age_18_21), n1 = nrow(age_8_10), n2 = nrow(age_18_21))   
jen_earlyado_vs_lateado <- cortest.jennrich(cor(age_11_13), cor(age_14_17), n1 = nrow(age_11_13), n2 = nrow(age_14_17))   
jen_earlyado_vs_adult <- cortest.jennrich(cor(age_11_13), cor(age_18_21), n1 = nrow(age_11_13), n2 = nrow(age_18_21))   
jen_lateado_vs_adult <- cortest.jennrich(cor(age_14_17), cor(age_18_21), n1 = nrow(age_14_17), n2 = nrow(age_18_21))   
ps <- c(jen_child_vs_earlyado$prob, jen_child_vs_lateado$prob, jen_child_vs_adult$prob,
        jen_earlyado_vs_lateado$prob, jen_earlyado_vs_adult$prob, jen_lateado_vs_adult$prob)
adjust_ps <- p.adjust(ps, method = "holm")

# Save off Jenrich's tests
sink(paste0(analyzed_data_dir, 'jenrich_test.txt')) # Write to file instead of the terminal
cat("Children vs. Early Adolescence\n")
cat("------------------------------\n")
adjust_ps[1]
cat("\n")
cat("Children vs. Late Adolescence\n")
cat("-----------------------------\n")
adjust_ps[2]
cat("\n")
cat("Children vs. Adults\n")
cat("-------------------\n")
adjust_ps[3]
cat("\n")
cat("Early Adolescence vs. Late Adolescence\n")
cat("--------------------------------------\n")
adjust_ps[4]
cat("\n")
cat("Early Adolescence vs. Adults\n")
cat("----------------------------\n")
adjust_ps[5]
cat("\n")
cat("Late Adolescence vs. Adults\n")
cat("---------------------------\n")
adjust_ps[6]
sink() # Stop writing to file

##### Fit and Evaluate Group-Level MDS ##### 

# Create (derived) dissimilarity matrix
diss_matrix <- toolbox_cormat %>% sim2diss()
diss_matrix

# Fit interval MDS (because items are metric)
fit_mds <- mds(diss_matrix, type = "interval") 

# Evaluate fit
plot(fit_mds, plot.type = "Shepard", main = "Interval MDS") # Interval transformation is well-fitting
fit_mds$stress # Overall stress is very low

# Assess goodness-of-fit
stressvec <- NULL
for (i in 1:4) { # 5 variables; can have a max of 4 dimensions
  # Refit MDS with 1-4 dimensions
  stressvec[i] <- mds(diss_matrix, ndim = i, type = "interval")$stress 
}
png(paste0(analyzed_data_dir, 'mds_scree.png'), width = 450, height = 450)
par(family = "Avenir", mgp = c(2.5, 1, 0))
plot(1:4, stressvec, 
     cex = 2.5, pch = 20, type = "b", lwd = 3,
     xlab = "Number of Dimensions", 
     ylab = "Stress-1", 
     main = "MDS Scree Plot (Full Sample)",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1, 
     font.main = 3) # According to elbow criterion, 2 dimensions are appropriate --> can move forward with 2D configuration plot
dev.off() 

# Look at stress-per-point (contributions of each variable to stress)
stress_per_point <- sort(fit_mds$spp, decreasing = TRUE)
stress_per_point # Pattern Comparison contributes the most to stress
plot(fit_mds, plot.type = "stressplot", main = "MDS Stress-Per-Point")
plot(fit_mds, plot.type = "bubbleplot", main = "MDS Bubble Plot")

# Check out configuration
png(paste0(analyzed_data_dir, 'mds_config.png'), width = 450, height = 450)
par(family = "Avenir", mgp = c(2.5, 1, 0))
plot(fit_mds, col = cols, cex = 2.5, label.conf = list(col = cols),
     main = "MDS Configuration (Full Sample)",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1,
     font.main = 3,
     xlim = c(-1, 1), ylim = c(-1, 1))
dev.off()

# Since we are using derived dissimilarities, examine MDS stability through bootstrapping
set.seed(123)
bootstrap <- bootmds(fit_mds, data = toolbox_only, method.dat = "pearson", nrep = 500)
png(paste0(analyzed_data_dir, 'mds_bootstrap.png'), width = 450, height = 450)
par(family = "Avenir", mgp = c(2.5, 1, 0))
plot(bootstrap, col = c('white'), ell = list(col = cols, lwd = 3),
     main = "MDS Configuration (Full Sample)",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1, 
     font.main = 3,
     xlim = c(-1, 1), ylim = c(-1, 1))
points(bootstrap$conf, col = cols, pch = 20, cex = 2.5)
dev.off()

##### Fit and Evaluate Age Group MDS (Adults as Target) #####

# Create (derived) dissimilarity matrices for each age group
child_diss_matrix <- cor(age_8_10) %>% sim2diss()
print("Children", quote = FALSE)
child_diss_matrix
earlyado_diss_matrix <- cor(age_11_13) %>% sim2diss()
print("Early Adolescence", quote = FALSE)
earlyado_diss_matrix
lateado_diss_matrix <- cor(age_14_17) %>% sim2diss()
print("Late Adolescence", quote = FALSE)
lateado_diss_matrix
adult_diss_matrix <- cor(age_18_21) %>% sim2diss()
print("Adults", quote = FALSE)
adult_diss_matrix

# Fit interval MDS (because items are metric) for each age group
fit_child_mds <- mds(child_diss_matrix, type = "interval") 
fit_earlyado_mds <- mds(earlyado_diss_matrix, type = "interval") 
fit_lateado_mds <- mds(lateado_diss_matrix, type = "interval") 
fit_adult_mds <- mds(adult_diss_matrix, type = "interval") 

# Evaluate fit across age
plot(fit_child_mds, plot.type = "Shepard", main = "Interval MDS") # Interval transformation is well-fitting
plot(fit_earlyado_mds, plot.type = "Shepard", main = "Interval MDS") # Interval transformation is well-fitting
plot(fit_lateado_mds, plot.type = "Shepard", main = "Interval MDS") # Interval transformation is well-fitting
plot(fit_adult_mds, plot.type = "Shepard", main = "Interval MDS") # Interval transformation is well-fitting

# Glance at overall stress across age
fit_child_mds$stress
fit_earlyado_mds$stress
fit_lateado_mds$stress
fit_adult_mds$stress

# Save off stress values from full sample and age group MDS solutions
sink(paste0(analyzed_data_dir, 'mds_stress.txt')) # Write to file instead of the terminal
cat("Sample\n")
cat("------\n")
fit_mds$stress
cat("\n\n")
cat("Children\n")
cat("--------\n")
fit_child_mds$stress
cat("\n\n")
cat("Early Adolescents\n")
cat("-----------------\n")
fit_earlyado_mds$stress
cat("\n\n")
cat("Late Adolescents\n")
cat("----------------\n")
fit_lateado_mds$stress
cat("\n\n")
cat("Adults\n")
cat("------\n")
fit_adult_mds$stress
sink() # Stop writing to file

# Align each MDS solution to the adult data using Procrustes: children
fit_child_adult <- smacof::Procrustes(X = fit_adult_mds$conf, Y = fit_child_mds$conf) 
fit_child_adult # Fixed adult, transformed child

# Align each MDS solution to the adult data using Procrustes: early adolescence
fit_earlyado_adult <- smacof::Procrustes(X = fit_adult_mds$conf, Y = fit_earlyado_mds$conf) 
fit_earlyado_adult # Fixed adult, transformed early adolescent

# Align each MDS solution to the adult data using Procrustes: late adolescence
fit_lateado_adult <- smacof::Procrustes(X = fit_adult_mds$conf, Y = fit_lateado_mds$conf) 
fit_lateado_adult # Fixed adult, transformed late adolescent

# Plot Procrustes configurations for each age group in a single plot
png(paste0(analyzed_data_dir, 'group_mds_config.png'), width = 450, height = 450)
par(family = "Avenir", mgp = c(2.5, 1, 0))
plot(fit_child_adult$X, col = age_18_21_purple,
     main = "MDS Configuration (Age Groups)", 
     xlab = "Dimension 1", ylab = "Dimension 2",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1,
     font.main = 3,
     xlim = c(-1, 1), ylim = c(-1, 1), pch = 15, cex = 2)
points(fit_child_adult$Yhat, col = age_8_10_orange, pch = 16, cex = 2)
points(fit_earlyado_adult$Yhat, col = age_11_13_red, pch = 17, cex = 2)
points(fit_lateado_adult$Yhat, col = age_14_17_pink, pch = 18, cex = 2.5)
text(fit_child_adult$X, labels = rownames(fit_child_adult$X), col = age_18_21_purple,
     cex = 1, pos = 2)
text(fit_child_adult$Yhat, labels = rownames(fit_child_adult$Yhat), col = age_8_10_orange,
     cex = 1, pos = 2)
text(fit_earlyado_adult$Yhat, labels = rownames(fit_earlyado_adult$Yhat), col = age_11_13_red,
     cex = 1, pos = 4)
text(fit_lateado_adult$Yhat, labels = rownames(fit_lateado_adult$Yhat), col = age_14_17_pink,
     cex = 1, pos = 4)
legend("bottom", col = c(age_8_10_orange, age_11_13_red, age_14_17_pink, age_18_21_purple), 
       pch = c(16, 17, 18, 15), cex = 1, pt.cex = c(2, 2, 2.5, 2), inset = .045, # Move legend vertically
       legend = c("8-10 Years (Aligned)", "11-13 Years (Aligned)", "14-17 Years (Aligned)", "18-21 Years (Target)"), ncol = 2)
dev.off()

# Save off Procrustes coefficients from above configurations
sink(paste0(analyzed_data_dir, 'group_mds_coef.txt')) # Write to file instead of the terminal
cat("Children vs. Adults\n")
cat("-------------------\n")
fit_child_adult
cat("\n\n")
cat("Early Adolescence vs. Adults\n")
cat("----------------------------\n")
fit_earlyado_adult
cat("\n\n")
cat("Late Adolescence vs. Adults\n")
cat("---------------------------\n")
fit_lateado_adult
sink() # Stop writing to file

# For a finer-grained interpretation, look at the pairwise differences from Procrustes
sink(paste0(analyzed_data_dir, 'group_mds_pair_diff.txt')) # Write to file instead of the terminal
cat("Children vs. Adults\n")
cat("-------------------\n")
fit_child_adult$pairdist %>% round(3)
cat("\n\n")
cat("Early Adolescence vs. Adults\n")
cat("----------------------------\n")
fit_earlyado_adult$pairdist %>% round(3)
cat("\n\n")
cat("Late Adolescence vs. Adults\n")
cat("---------------------------\n")
fit_lateado_adult$pairdist %>% round(3)
sink() # Stop writing to file

##### Compare Age Group Distances in MDS #####

# If the points in the configuration become more distinct with age, then the distances should move away from 0

# Function to calculate Euclidean distance
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

# Function to create matrix of Euclidean distances
mat <- function(matrix) {
  
  # Get x and y (D1 and D2) coordinates
  flanker_x = matrix[1, 1]
  flanker_y = matrix[1, 2]
  dccs_x = matrix[2, 1]
  dccs_y = matrix[2, 2]
  pattern_x = matrix[3, 1]
  pattern_y = matrix[3, 2]
  psmt_x = matrix[4, 1]
  psmt_y = matrix[4, 2]
  list_x = matrix[5, 1]
  list_y = matrix[5, 2]
  
  # Create data frame of pairwise Euclidean distances
  Flnk <- 
    c(0,
      euclidean_distance(flanker_x, flanker_y, dccs_x, dccs_y),
      euclidean_distance(flanker_x, flanker_y, pattern_x, pattern_y),
      euclidean_distance(flanker_x, flanker_y, psmt_x, psmt_y),
      euclidean_distance(flanker_x, flanker_y, list_x, list_y))
  DCCS <- 
    c(euclidean_distance(dccs_x, dccs_y, flanker_x, flanker_y),
      0,
      euclidean_distance(dccs_x, dccs_y, pattern_x, pattern_y),
      euclidean_distance(dccs_x, dccs_y, psmt_x, psmt_y),
      euclidean_distance(dccs_x, dccs_y, list_x, list_y))
  PnCn <- 
    c(euclidean_distance(pattern_x, pattern_y, flanker_x, flanker_y),
      euclidean_distance(pattern_x, pattern_y, dccs_x, dccs_y),
      0,
      euclidean_distance(pattern_x, pattern_y, psmt_x, psmt_y),
      euclidean_distance(pattern_x, pattern_y, list_x, list_y))
  PSMT <- 
    c(euclidean_distance(psmt_x, psmt_y, flanker_x, flanker_y),
      euclidean_distance(psmt_x, psmt_y, dccs_x, dccs_y),
      euclidean_distance(psmt_x, psmt_y, pattern_x, pattern_y),
      0,
      euclidean_distance(psmt_x, psmt_y, list_x, list_y))
  LtSg <- 
    c(euclidean_distance(list_x, list_y, flanker_x, flanker_y),
      euclidean_distance(list_x, list_y, dccs_x, dccs_y),
      euclidean_distance(list_x, list_y, pattern_x, pattern_y),
      euclidean_distance(list_x, list_y, psmt_x, psmt_y),
      0)
  df <- rbind(Flnk, DCCS, PnCn, PSMT, LtSg)
}

# Create matrices of Euclidean distances for each age group
child_dist = mat(fit_child_adult$Yhat)
colnames(child_dist) <- short_task_labels
earlyado_dist = mat(fit_earlyado_adult$Yhat)
colnames(earlyado_dist) <- short_task_labels
lateado_dist = mat(fit_lateado_adult$Yhat)
colnames(lateado_dist) <- short_task_labels
adult_dist = mat(fit_child_adult$X)
colnames(adult_dist) <- short_task_labels

# Create single distance data frame
child_dist_df <- melt(child_dist)
child_dist_df$AgeGroup <- c("8-10 Years")
earlyado_dist_df <- melt(earlyado_dist)
earlyado_dist_df$AgeGroup <- c("11-13 Years")
lateado_dist_df <- melt(lateado_dist)
lateado_dist_df$AgeGroup <- c("14-17 Years")
adult_dist_df <- melt(adult_dist)
adult_dist_df$AgeGroup <- c("18-21 Years")
dist_df <- rbind(child_dist_df, earlyado_dist_df, lateado_dist_df, adult_dist_df)
dist_df <- dist_df[dist_df$value != 0, ]

# Create Task and Cluster columns for distance data frame
dist_df$Task <- paste0(dist_df$Var1, "-", dist_df$Var2)
dups <- c("DCCS-Flnk", "PnCn-Flnk", "PSMT-Flnk", "LtSg-Flnk",
          "PnCn-DCCS", "PSMT-DCCS", "LtSg-DCCS",
          "PSMT-PnCn", "LtSg-PnCn",
          "LtSg-PSMT")
dist_df <- dist_df[!(dist_df$Task %in% dups), ]
dist_df <- dist_df %>%
  mutate(Cluster = case_when(
    Task == "Flnk-DCCS" ~ "Within",
    Task == "Flnk-PnCn" ~ "Within",
    Task == "DCCS-PnCn" ~ "Within",
    Task == "PSMT-LtSg" ~ "Within",
    .default ="Between"
  ))
str(dist_df)
dist_df$AgeGroup <- factor(dist_df$AgeGroup, levels = c("8-10 Years", "11-13 Years", "14-17 Years", "18-21 Years"))
dist_df$Task <- as.factor(dist_df$Task)
dist_df$Cluster <- as.factor(dist_df$Cluster)
dist_df <- dist_df %>% dplyr::rename(Distance = value)
str(dist_df)

# Visually examine age-related changes in distances
summ_age <- dist_df %>%
  group_by(AgeGroup) %>%
  dplyr::summarise(avg = mean(Distance)) 
summ_cluster <- dist_df %>%
  group_by(Cluster) %>%
  dplyr::summarise(avg = mean(Distance)) 
summ_age_cluster <- dist_df %>%
  group_by(AgeGroup, Cluster) %>%
  dplyr::summarise(avg = mean(Distance)) 
ggplot(data = summ_age, aes(x = AgeGroup, y = avg, color = AgeGroup, fill = AgeGroup)) +
  geom_hline(yintercept = seq(.25, 1.25, by = .25), colour = 'grey90') +
  scale_y_continuous(breaks = seq(.25, 1.25, by = .25)) +
  geom_point(size = 6) +
  scale_color_manual(values = c(age_8_10_orange, age_11_13_red, age_14_17_pink, age_18_21_purple)) +
  labs(x = "Age Group", y = "Distance") +
  plot_theme + theme(legend.position = "none")
ggsave(paste0(analyzed_data_dir, 'age_distances.png'), width = 7, height = 7)
ggplot(data = summ_cluster, aes(x = Cluster, y = avg, color = Cluster, fill = Cluster)) +
  geom_hline(yintercept = seq(.25, 1.25, by = .25), colour = 'grey90') +
  scale_y_continuous(breaks = seq(.25, 1.25, by = .25)) +
  geom_point(size = 6) +
  scale_color_manual(values = c(reg_gray, reg_brown)) +
  labs(x = "Cluster", y = "Distance") +
  plot_theme + theme(legend.position = "none")
ggsave(paste0(analyzed_data_dir, 'cluster_distances.png'), width = 4.5, height = 7)
ggplot(data = summ_age_cluster, aes(x = AgeGroup, y = avg, group = Cluster, color = Cluster, fill = Cluster)) +
  geom_hline(yintercept = seq(.25, 1.25, by = .25), colour = 'grey90') +
  scale_y_continuous(breaks = seq(.25, 1.25, by = .25)) +
  geom_point(size = 6) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c(reg_gray, reg_brown)) +
  labs(x = "Age Group", y = "Distance") +
  plot_theme
ggsave(paste0(analyzed_data_dir, 'age_cluster_distances.png'), width = 9, height = 7)

# Prepare to bootstrap SEs
samp_8_10 = floor(nrow(age_8_10) * .95)
samp_11_13 = floor(nrow(age_11_13) * .95)
samp_14_17 = floor(nrow(age_14_17) * .95)
samp_18_21 = floor(nrow(age_18_21) * .95)
boot_dist_df <- data.frame()

# Bootstrap SEs across 1000 iterations
for (i in 1:1000) {
  
  # Create random child dissimilarity matrix (sample 95% of age group with replacement)
  set.seed((i*100) + 1) 
  temp_8_10 <- age_8_10[sample(nrow(age_8_10), samp_8_10, replace = TRUE), ]
  temp_child_diss_matrix <- cor(temp_8_10) %>% sim2diss()
  
  # Create random early adolescent dissimilarity matrix (sample 95% of age group with replacement)
  set.seed((i*100) + 2) 
  temp_11_13 <- age_11_13[sample(nrow(age_11_13), samp_11_13, replace = TRUE), ]
  temp_earlyado_diss_matrix <- cor(temp_11_13) %>% sim2diss()
  
  # Create random late adolescent dissimilarity matrix (sample 95% of age group with replacement)
  set.seed((i*100) + 3) 
  temp_14_17 <- age_14_17[sample(nrow(age_14_17), samp_14_17, replace = TRUE), ]
  temp_lateado_diss_matrix <- cor(temp_14_17) %>% sim2diss()
  
  # Create random adult dissimilarity matrix (sample 95% of age group with replacement)
  set.seed((i*100) + 4)
  temp_18_21 <- age_18_21[sample(nrow(age_18_21), samp_18_21, replace = TRUE), ]
  temp_adult_diss_matrix <- cor(temp_18_21) %>% sim2diss()
  
  # Fit interval MDS for each random subsample
  temp_fit_child_mds <- mds(temp_child_diss_matrix, type = "interval") 
  temp_fit_earlyado_mds <- mds(temp_earlyado_diss_matrix, type = "interval") 
  temp_fit_lateado_mds <- mds(temp_lateado_diss_matrix, type = "interval") 
  temp_fit_adult_mds <- mds(temp_adult_diss_matrix, type = "interval") 
  
  # Fit Procrustes solutions for each random subsample
  temp_fit_child_adult <- smacof::Procrustes(X = temp_fit_adult_mds$conf, Y = temp_fit_child_mds$conf) 
  temp_fit_earlyado_adult <- smacof::Procrustes(X = temp_fit_adult_mds$conf, Y = temp_fit_earlyado_mds$conf) 
  temp_fit_lateado_adult <- smacof::Procrustes(X = temp_fit_adult_mds$conf, Y = temp_fit_lateado_mds$conf) 
  
  # Create distance data frame for random subsample of children
  temp_child_dist = mat(temp_fit_child_adult$Yhat)
  colnames(temp_child_dist) <- short_task_labels
  temp_child_dist_df <- melt(temp_child_dist)
  temp_child_dist_df$AgeGroup <- c("8-10 Years")
  temp_child_dist_df$Run <- c(i)
  
  # Create distance data frame for random subsample of early adolescents
  temp_earlyado_dist = mat(temp_fit_earlyado_adult$Yhat)
  colnames(temp_earlyado_dist) <- short_task_labels
  temp_earlyado_dist_df <- melt(temp_earlyado_dist)
  temp_earlyado_dist_df$AgeGroup <- c("11-13 Years")
  temp_earlyado_dist_df$Run <- c(i)
  
  # Create distance data frame for random subsample of late adolescents
  temp_lateado_dist = mat(temp_fit_lateado_adult$Yhat)
  colnames(temp_lateado_dist) <- short_task_labels
  temp_lateado_dist_df <- melt(temp_lateado_dist)
  temp_lateado_dist_df$AgeGroup <- c("14-17 Years")
  temp_lateado_dist_df$Run <- c(i)
  
  # Create distance data frame for random subsample of adults
  temp_adult_dist = mat(temp_fit_child_adult$X)
  colnames(temp_adult_dist) <- short_task_labels
  temp_adult_dist_df <- melt(temp_adult_dist)
  temp_adult_dist_df$AgeGroup <- c("18-21 Years")
  temp_adult_dist_df$Run <- c(i)
  
  # Combine random subsamples together into single distance data frame
  boot_dist_df <- rbind(boot_dist_df, temp_child_dist_df, temp_earlyado_dist_df, temp_lateado_dist_df, temp_adult_dist_df)
  
}

# Create Task and Cluster columns for bootstrapped distance data frame
boot_dist_df <- boot_dist_df[boot_dist_df$value != 0, ]
boot_dist_df$Task <- paste0(boot_dist_df$Var1, "-", boot_dist_df$Var2)
boot_dist_df <- boot_dist_df[!(boot_dist_df$Task %in% dups), ]
boot_dist_df <- boot_dist_df %>%
  mutate(Cluster = case_when(
    Task == "Flnk-DCCS" ~ "Within",
    Task == "Flnk-PnCn" ~ "Within",
    Task == "DCCS-PnCn" ~ "Within",
    Task == "PSMT-LtSg" ~ "Within",
    .default ="Between"
  ))
str(boot_dist_df)
boot_dist_df$AgeGroup <- factor(boot_dist_df$AgeGroup, levels = c("8-10 Years", "11-13 Years", "14-17 Years", "18-21 Years"))
boot_dist_df$Task <- as.factor(boot_dist_df$Task)
boot_dist_df$Cluster <- as.factor(boot_dist_df$Cluster)
boot_dist_df <- boot_dist_df %>% dplyr::rename(Distance = value)
str(boot_dist_df)
dim(boot_dist_df) # Should be 1000 (iterations) * 10 (task combinations) * 4 (age groups) long

# Visually examine age-related changes in distances with bootstrapped SEs
temp_boot_summ_age <- summarySE(data = boot_dist_df, measurevar = "Distance", groupvars = c("AgeGroup", "Run"))
boot_summ_age <- summarySE(data = temp_boot_summ_age, measurevar = "Distance", groupvars = "AgeGroup")
temp_boot_summ_cluster <- summarySE(data = boot_dist_df, measurevar = "Distance", groupvars = c("Cluster", "Run"))
boot_summ_cluster <- summarySE(data = temp_boot_summ_cluster, measurevar = "Distance", groupvars = "Cluster")
temp_boot_summ_age_cluster <- summarySE(data = boot_dist_df, measurevar = "Distance", groupvars = c("AgeGroup", "Cluster", "Run"))
boot_summ_age_cluster <- summarySE(data = temp_boot_summ_age_cluster, measurevar = "Distance", groupvars = c("AgeGroup", "Cluster"))
# SEs too small to show up on plots --> not including 

# Save off bootstrapped SEs
sink(paste0(analyzed_data_dir, 'boot_ses.txt')) # Write to file instead of the terminal
print(boot_summ_age)
cat("-------\n")
print(boot_summ_cluster)
cat("-------\n")
print(boot_summ_age_cluster)
sink() # Stop writing to file

##### Explore 3 Group-Level Factors with EFA ##### 

# Scree plot - suggests that 1 factor should be extracted
scree(toolbox_cormat, factors = FALSE)

# Parallel analysis - suggests that 2 factors should be extracted
set.seed(123)
efa_pa <- fa.parallel(toolbox_cormat, fa = "both", fm = "ml", n.obs = nrow(toolbox_only))

# Given scree plot and parallel analysis, do not pursue 3-factor EFA, despite MDS configuration suggesting that there may be 3 clusters of fluid ability measures

# Consider different numbers of factors across age - children
age_8_10_cormat <- cor(age_8_10) # Compute Pearson correlations
scree(age_8_10_cormat, factors = FALSE) # Scree plot - suggests that 2 factors should be extracted
set.seed(123)
age_8_10_efa_pa <- fa.parallel(age_8_10_cormat, fa = "both", fm = "ml", n.obs = nrow(age_8_10)) # Parallel analysis - suggests that 2 factors should be extracted

# Consider different numbers of factors across age - early adolescence
age_11_13_cormat <- cor(age_11_13) # Compute Pearson correlations
scree(age_11_13_cormat, factors = FALSE) # Scree plot - suggests that 2 factors should be extracted
set.seed(123)
age_11_13_efa_pa <- fa.parallel(age_11_13_cormat, fa = "both", fm = "ml", n.obs = nrow(age_11_13)) # Parallel analysis - suggests that 2 factors should be extracted

# Consider different numbers of factors across age - late adolescence
age_14_17_cormat <- cor(age_14_17) # Compute Pearson correlations
scree(age_14_17_cormat, factors = FALSE) # Scree plot - suggests that 1 factor should be extracted
set.seed(123)
age_14_17_efa_pa <- fa.parallel(age_14_17_cormat, fa = "both", fm = "ml", n.obs = nrow(age_14_17)) # Parallel analysis - suggests that 2 factors should be extracted

# Consider different numbers of factors across age - adults
age_18_21_cormat <- cor(age_18_21) # Compute Pearson correlations
scree(age_18_21_cormat, factors = FALSE) # Scree plot - suggests that 2 factors should be extracted
set.seed(123)
age_18_21_efa_pa <- fa.parallel(age_18_21_cormat, fa = "both", fm = "ml", n.obs = nrow(age_18_21)) # Parallel analysis - suggests that 2 factors should be extracted

# After breaking down scree plots and parallel analyses by age group --> overall factor structure does not change with age

##### Explore 2 Group-Level Factors with EFA ##### 

# Based MDS configuration and 3-factor EFA evaluation, fit 2-factor EFA with an oblique rotation for better interpretability (oblimin rotation, ML estimation)
#   Use Pearson correlation matrix from earlier because of metric data
#   Allow correlations among factors (we do not assume that the factors are independent)
efa_2_rot <- fa(toolbox_cormat, nfactors = 2, rotate = "oblimin", fm = "ml")
print(efa_2_rot$loadings, cutoff = 0.1) # Flanker, DCCS, and Pattern load onto factor 1, PSMT and List load onto factor 2

# Examine correlation matrix from rotated EFA fit
round(efa_2_rot$Phi, 3) # Factors 1 and 2 are highly correlated

# Amount of explained variance - suggests that 2 factors should be extracted (substantial gains in Cumulative Var from 1 to 2 factors)
efa_2_rot

# Interpretability - suggests that 2 factors should be extracted
png(paste0(analyzed_data_dir, 'efa_2factor_loadings.png'), width = 450, height = 450)
par(family = "Avenir", mgp = c(2.5, 1, 0))
plot(efa_2_rot$loadings, col = cols,
     main = "2-Factor EFA Loadings",
     xlab = "Factor 1: Cognitive Control", 
     ylab = "Factor 2: Memory", 
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1,
     font.main = 3,
     asp = 1, xlim = c(-0.5, 1), ylim = c(-0.5, 1), 
     type = "n", pch = 20, cex = 6)
abline(h = 0, v = 0, col = 'black', lty = 2)
points(efa_2_rot$loadings, col = cols, pch = 20, cex = 6)
dev.off()

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'efa.txt'))

cat("LOADINGS\n")
print(efa_2_rot$loadings, cutoff = 0.1)

cat("\nFACTOR CORRELATIONS\n\n")
round(efa_2_rot$Phi, 3)

cat("\nEXPLAINED VARIANCE\n\n")
efa_2_rot

# Stop writing to file
sink()
closeAllConnections()

##### Confirm 2 Group-Level Factors with CFA ##### 

# Specify CFA model
cfa_2 <- 'CgCn =~ Flnk + DCCS + PnCn
          Mmry =~ PSMT + LtSg'

# Fit CFA
#   Feed in dataframe; correlation matrix is computed internally
#   Treat all endogenous variables as numeric, not ordinal (we only have metric variables here)
fit_cfa_2 <- lavaan::cfa(cfa_2, data = toolbox_only, ordered = FALSE) 

# Plot CFA and show constituent matrices
#   With unstandardized parameters
png(paste0(analyzed_data_dir, 'cfa_2factor_est_loadings.png'), width = 450, height = 550)
par(family = "Avenir")
semPaths(fit_cfa_2, what = "est", edge.label.cex = 0.7, edge.color = 1, esize = 1, 
         sizeMan = 5, asize = 2.5, intercepts = FALSE, rotation = 4, mar = c(1, 5, 1.5, 5), 
         fade = FALSE, nCharNodes = 4, label.color = c(cols, 'grey80', 'grey80'), border.color = c(cols, 'grey80', 'grey80'))
title("2-Factor CFA with Estimated Parameters", cex.main = 1.5, font.main = 3, line = 3, 
      col.main = 'black')
dev.off()

# Inspect the: 
inspect(fit_cfa_2, what = "est")$lambda # Loadings matrix (of indicators on latent variables)
inspect(fit_cfa_2, what = "est")$theta # Error variance-covariance matrix
inspect(fit_cfa_2, what = "est")$psi # Factor variance-covariance matrix

# Plot CFA and show constituent matrices
#   With standardized parameters (i.e., loadings can be compared to each other)
png(paste0(analyzed_data_dir, 'cfa_2factor_std_loadings.png'), width = 450, height = 550)
par(family = "Avenir")
semPaths(fit_cfa_2, what = "std", edge.label.cex = 0.7, edge.color = 1, esize = 1, 
         sizeMan = 5, asize = 2.5, intercepts = FALSE, rotation = 4, mar = c(1, 5, 1.5, 5), 
         fade = FALSE, nCharNodes = 4, label.color = c(cols, 'grey80', 'grey80'), border.color = c(cols, 'grey80', 'grey80'))
title("2-Factor CFA with Standardized Parameters", cex.main = 1.5, font.main = 3, line = 3, 
      col.main = 'black')
dev.off()

# Inspect the: 
inspect(fit_cfa_2, what = "std")$lambda # Loadings matrix (of indicators on latent variables)
inspect(fit_cfa_2, what = "std")$theta # Error variance-covariance matrix
inspect(fit_cfa_2, what = "std")$psi # Factor variance-covariance matrix

# Consider full summary output and goodness-of-fit 
summary(fit_cfa_2, fit.measures = TRUE)
# Chi-square: non-significant p-value suggests that model fits (0.589 - looks good)
# CFI: value >= .95 suggests that model fits (1.000 - looks good)
# RMSEA: value <= .05 suggests that model fits (0.000 - looks good)
# RMSEA upper CI bound: value <= .10 suggests that model fits (0.040 - looks good)
# SRMR: value <= .08 suggests that model fits (0.008 - looks good)
# Overall, BEAUTIFUL fit!

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'sample_cfa.txt'))

cat("(STANDARDIZED) LOADINGS\n\n")
inspect(fit_cfa_2, what = "std")$lambda

cat("\n(STANDARDIZED) VARIANCES\n\n")
inspect(fit_cfa_2, what = "std")$theta

cat("\n(STANDARDIZED) COVARIANCES\n\n")
inspect(fit_cfa_2, what = "std")$psi

cat("\nFULL (UNSTANDARDIZED) SUMMARY OUTPUT\n\n")
summary(fit_cfa_2, fit.measures = TRUE)

# Stop writing to file
sink()

# Extract and save the factor scores
cfa_2_scores <- lavPredict(fit_cfa_2)   
dim(cfa_2_scores)
head(cfa_2_scores)            
cfa_2_scores <- cbind(cfa_2_scores, uncorr_toolbox[, c("Subject", "Age")])
write.csv(cfa_2_scores, paste0(analyzed_data_dir, "cfa_2factor_scores.csv"))

##### Confirm 2 Factors (Age as Multiple Groups) with CFA ##### 

# Rename columns for CFA plot
toolbox_age_group <- uncorr_toolbox[, c("Age_Group", toolbox_z)]
names(toolbox_age_group) <- c("Age Group", short_task_labels)

# Test for differences in CORRELATIONS rather than COVARIANCES among variables (i.e., set the SD of the latent variables to 1)

# Fit invariance models, following from the 2 factor CFA
fit_weak <- lavaan::cfa(cfa_2, data = toolbox_age_group, group = "Age Group", group.equal = "loadings", std.lv = TRUE)
summary(fit_weak) # loadings are constrained to be = across age groups
lavInspect(fit_weak, what = "cov.lv") # loadings are constant across age groups, but correlations (off-diagonal elements) are not
fit_strong <- lavaan::cfa(cfa_2, data = toolbox_age_group, group = "Age Group", group.equal = c("loadings", "lv.covariances"), std.lv = TRUE)
summary(fit_strong) # loadings AND latent variable correlations are constrained to be = across age groups
lavInspect(fit_strong, what = "cov.lv") # loadings AND correlations (off-diagonal elements) are constant across age groups

# Conduct LR-test
#   NOTE: H0 = no difference between models
#         H1 = difference between models
lavTestLRT(fit_weak, fit_strong)
# Weak model: loadings are constant across age groups, correlations can change across age groups
# Strong model: loadings are constant across age groups, correlations are constant across age groups
# No difference between models --> changes in correlations across age is no more likely --> no evidence for age variability in CgCn-Mmry correlations

# Consider full summary outputs and goodness-of-fit 
summary(fit_weak, fit.measures = TRUE) 
# Chi-square: non-significant p-value suggests that model fits (0.500 - looks good)
# CFI: value >= .95 suggests that model fits (1.000 - looks good)
# RMSEA: value <= .05 suggests that model fits (0.000 - looks good)
# RMSEA upper CI bound: value <= .10 suggests that model fits (0.048 - looks good)
# SRMR: value <= .08 suggests that model fits (0.034 - looks good)
# Overall, BEAUTIFUL fit!
summary(fit_strong, fit.measures = TRUE) 
# Chi-square: non-significant p-value suggests that model fits (0.573 - looks good)
# CFI: value >= .95 suggests that model fits (1.000 - looks good)
# RMSEA: value <= .05 suggests that model fits (0.000 - looks good)
# RMSEA upper CI bound: value <= .10 suggests that model fits (0.043 - looks good)
# SRMR: value <= .08 suggests that model fits (0.039 - looks good)
# Overall, BEAUTIFUL fit!

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'multigroup_cfa.txt'))

cat("LIKELIHOOD RATIO TEST\n")
lavTestLRT(fit_weak, fit_strong)

cat("\nFULL (WEAK INVARIANCE) SUMMARY OUTPUT\n\n")
summary(fit_weak, fit.measures = TRUE) 

cat("\nFULL (STRONG INVARIANCE) SUMMARY OUTPUT\n\n")
summary(fit_strong, fit.measures = TRUE) 

# Stop writing to file
sink()

##### Perform Power Analysis for Multigroup CFA ##### 

# Resources:
# MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and determination of 
#   sample size for covariance structure modeling. Psychological Methods, 1, 130-149. 
# Zhang, Z., & Yuan, K.-H. (2018). Practical Statistical Power Analysis Using Webpower and R (Eds). 
#   Granger, IN: ISDSA Press.

# Why use this strategy?
# We only perform global LR-testing:
# lavTestLRT(fit_weak, fit_strong)
# We are not testing for differences in specific parameters. This RMSEA-based power analysis 
# evaluates how likely the LR-test is to detect a meaningful level of model misfit (e.g., RMSEA = .10) 
# given the sample size and model degrees of freedom, under the assumption that a smaller RMSEA value 
# (e.g., .05) represents an acceptable or "close" fit. It quantifies the sensitivity of the overall 
# model-fit test, not the precision of individual parameter estimates.

# Double-check the difference in df
lavTestLRT(fit_weak, fit_strong) # 3

# Double-check the total sample size n
lavInspect(fit_weak, "nobs") |> sum() # 1049

# Set rmsea0 = 0.05, representing an acceptable fit, and vary the misfit rmsea1 systematically
rmsea1_vec <- c(0.08, 0.09, 0.10, 0.11, 0.12)
# e.g., if we are super strict, rmsea1 = 0.08 means that we would consider such a RMSEA as a misfit already 

# Input n and df into power analysis, plus systematic rmsea values
power_vals <- sapply(rmsea1_vec, function(r1) wp.sem.rmsea(n = 1049, df = 3, rmsea0 = 0.05, rmsea1 = r1, type = "close")$power)
names(power_vals) <- paste("RMSEA1", rmsea1_vec)
power_vals |> round(3)
# Verdict: a power analysis for the LR-test (delta df = 3, n = 1049) indicated adequate power 
# (> .80) to detect a meaningful deterioration in fit corresponding to RMSEA >= 0.10, while
# smaller deviations (RMSEA = 0.08) would be detected with moderate power (~.50).

# Save off power analysis
sink(paste0(analyzed_data_dir, 'lr_power.txt')) # Write to file instead of the terminal
cat("LR-Test\n")
cat("-------\n")
lavTestLRT(fit_weak, fit_strong)
cat("\n\n")
cat("rmsea1 Parameter Values\n")
cat("-----------------------\n")
rmsea1_vec
cat("\n\n")
cat("Power Values\n")
cat("------------\n")
power_vals |> round(3)
sink() # Stop writing to file
