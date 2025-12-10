##### Calculate HCPD FC Demographic Data Summaries #####

# Written by: Camille Phaneuf-Hadd (cphaneuf@g.harvard.edu)
# Last updated: 12/4/25
# Question: who participated in HCPD?

# Inputs: demographic data
# Computes: 
# - visual summary of sample
# - statistical summary of sample
# Outputs: 
# - png plot of age-sex distribution into results/demog directory
# - txt files of demographic stats into results/demog directory

##### Set up Script ##### 

# Load needed libraries
require(pacman) # for p_load()
p_load(tidyverse, # for df manipulation
       dplyr) # for %>% and other operators

# Load shared HCPD FC functions
source("utilities.R")

# Set paths
getwd()
processed_data_dir <- "../data/"
analyzed_data_dir <- "../results/demog/"

# Read in demographic data
demog <- read.csv(paste0(processed_data_dir, "uncorr_toolbox.csv"))

# Get floored age and age groups
demog$AgeFloor <- floor(demog$Age)
demog <- demog %>% mutate(AgeGroup = case_when((AgeFloor <= 10) ~ "8-10 Years",
                                               (AgeFloor >= 11 & AgeFloor <= 13) ~ "11-13 Years",
                                               (AgeFloor >= 14 & AgeFloor <= 17) ~ "14-17 Years",
                                               (AgeFloor >= 18) ~ "18-21 Years"))

# Create data frames for each age group
age_8_10 <- demog[demog$AgeGroup == "8-10 Years", ]
age_11_13 <- demog[demog$AgeGroup == "11-13 Years", ]
age_14_17 <- demog[demog$AgeGroup == "14-17 Years", ]
age_18_21 <- demog[demog$AgeGroup == "18-21 Years", ]

##### Save Age-Sex Distribution ##### 

demog <- demog %>% mutate(Sex = case_when((Sex == "F") ~ "Female",
                                          (Sex == "M") ~ "Male")) 
ggplot(data = demog, aes(x = AgeFloor, fill = Sex, color = Sex)) +
  scale_x_continuous(breaks = c(8:21)) +
  geom_hline(yintercept = seq(0, 100, by = 25), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) + 
  scale_fill_manual(values = c(light_gray, reg_gray)) + 
  scale_color_manual(values = c(light_gray, reg_gray)) + 
  geom_histogram(binwidth = 1, alpha = .75) +
  geom_vline(xintercept = c(10.5, 13.5, 17.5), colour = 'black', linetype = 'dashed') +
  labs(x = "Age (Years)", y = "Number of Participants") +
  plot_theme
ggsave(paste0(analyzed_data_dir, "sample_dist.png"), plot = last_plot(), width = 9, height = 7)

##### Save Age 8-10 Summaries ##### 

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'age_8-10.txt'))

cat("8-10 Years\n")
cat("----------\n")
cat("Total N:", length(age_8_10$Subject), "\n")
cat("Total Female:", sum(age_8_10$Sex == "F"), "\n")
cat("Total Male:", sum(age_8_10$Sex == "M"), "\n")
cat("Mean Age:", mean(age_8_10$Age), "\n")
cat("Standard Deviation Age:", sd(age_8_10$Age), "\n")

# Stop writing to file
sink()

##### Save Age 11-13 Summaries ##### 

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'age_11-13.txt'))

cat("11-13 Years\n")
cat("-----------\n")
cat("Total N:", length(age_11_13$Subject), "\n")
cat("Total Female:", sum(age_11_13$Sex == "F"), "\n")
cat("Total Male:", sum(age_11_13$Sex == "M"), "\n")
cat("Mean Age:", mean(age_11_13$Age), "\n")
cat("Standard Deviation Age:", sd(age_11_13$Age), "\n")

# Stop writing to file
sink()

##### Save Age 14-17 Summaries ##### 

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'age_14-17.txt'))

cat("14-17 Years\n")
cat("-----------\n")
cat("Total N:", length(age_14_17$Subject), "\n")
cat("Total Female:", sum(age_14_17$Sex == "F"), "\n")
cat("Total Male:", sum(age_14_17$Sex == "M"), "\n")
cat("Mean Age:", mean(age_14_17$Age), "\n")
cat("Standard Deviation Age:", sd(age_14_17$Age), "\n")

# Stop writing to file
sink()

##### Save Age 18-21 Summaries ##### 

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'age_18-21.txt'))

cat("18-21 Years\n")
cat("-----------\n")
cat("Total N:", length(age_18_21$Subject), "\n")
cat("Total Female:", sum(age_18_21$Sex == "F"), "\n")
cat("Total Male:", sum(age_18_21$Sex == "M"), "\n")
cat("Mean Age:", mean(age_18_21$Age), "\n")
cat("Standard Deviation Age:", sd(age_18_21$Age), "\n")

# Stop writing to file
sink()

##### Save Demographic Summaries ##### 

# Reporting the same information as Tervo-Clemmens et al. (2023) Supplementary Table S1
# (https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-42540-8/MediaObjects/41467_2023_42540_MOESM1_ESM.pdf)

N = length(demog$Subject)

# Write to file instead of the terminal
sink(paste0(analyzed_data_dir, 'sample_demog.txt'))

cat("All Participants\n")
cat("----------------\n")
cat("Total N:", length(demog$Subject), "\n")

cat("\nAGE\n")
cat("---\n")
cat("Youngest Age:", min(demog$Age), "\n")
cat("Oldest Age:", max(demog$Age), "\n")
cat("Mean Age:", mean(demog$Age), "\n")
cat("Standard Deviation Age:", sd(demog$Age), "\n")
cat("Percent 8-10 Years:", (length(age_8_10$Subject) / N) * 100, "\n")
cat("Percent 11-13 Years:", (length(age_11_13$Subject) / N) * 100, "\n")
cat("Percent 14-17 Years:", (length(age_14_17$Subject) / N) * 100, "\n")
cat("Percent 18-21 Years:", (length(age_18_21$Subject) / N) * 100, "\n")

cat("\nSEX\n")
cat("---\n")
cat("Percent Female:", (sum(demog$Sex == "Female") / N) * 100, "\n")
cat("Percent Male:", (sum(demog$Sex == "Male") / N) * 100, "\n")

cat("\nRACE\n")
cat("----\n")
cat("Percent Asian:", (sum(demog$Race == "Asian") / N) * 100, "\n")
cat("Percent American Indian/Alaska Native:", (sum(demog$Race == "American Indian/Alaska Native") / N) * 100, "\n")
cat("Percent Black or African American:", (sum(demog$Race == "Black or African American") / N) * 100, "\n")
cat("Percent Hawaiian or Pacific Islander:", (sum(demog$Race == "Hawaiian or Pacific Islander") / N) * 100, "\n")
cat("Percent White:", (sum(demog$Race == "White") / N) * 100, "\n")
cat("Percent More than one race:", (sum(demog$Race == "More than one race") / N) * 100, "\n")
cat("Percent Unknown or not reported:", (sum(demog$Race == "Unknown or not reported") / N) * 100, "\n")

cat("\nETHNICITY\n")
cat("---------\n")
cat("Percent Hispanic or Latino:", (sum(demog$Ethnicity == "Hispanic or Latino") / N) * 100, "\n")
cat("Percent Not Hispanic or Latino:", (sum(demog$Ethnicity == "Not Hispanic or Latino") / N) * 100, "\n")
cat("Percent Unknown or not reported:", ((sum(demog$Ethnicity == "Unknown or not reported") + sum(demog$Ethnicity == "unknown or not reported")) / N) * 100, "\n")

cat("\nFAMILY INCOME\n")
cat("------------\n")
cat("Percent <25k:", (sum(demog$Income == "<25k", na.rm = TRUE) / N) * 100, "\n")
cat("Percent 25k-49k:", (sum(demog$Income == "25k-49k", na.rm = TRUE) / N) * 100, "\n")
cat("Percent 50k-74k:", (sum(demog$Income == "50k-74k", na.rm = TRUE) / N) * 100, "\n")
cat("Percent 75k-99k:", (sum(demog$Income == "75k-99k", na.rm = TRUE) / N) * 100, "\n")
cat("Percent 100k-199k:", (sum(demog$Income == "100k-199k", na.rm = TRUE) / N) * 100, "\n")
cat("Percent 200k+:", (sum(demog$Income == "200k+", na.rm = TRUE) / N) * 100, "\n")
cat("Percent Unknown or not reported:", (sum(is.na(demog$Income)) / N) * 100, "\n")

cat("\nGUARDIAN HIGHEST EDUCATION\n")
cat("--------------------------\n")
cat("Percent Incomplete High School:", (sum(demog$Education == "Incomplete High School", na.rm = TRUE) / N) * 100, "\n")
cat("Percent High School:", (sum(demog$Education == "High School", na.rm = TRUE) / N) * 100, "\n")
cat("Percent 1-3 Years of College:", (sum(demog$Education == "1-3 Years of College", na.rm = TRUE) / N) * 100, "\n")
cat("Percent Bachelor's:", (sum(demog$Education == "Bachelor's", na.rm = TRUE) / N) * 100, "\n")
cat("Percent Postgraduate:", (sum(demog$Education == "Postgraduate", na.rm = TRUE) / N) * 100, "\n")
cat("Percent Unknown or not reported:", (sum(is.na(demog$Education)) / N) * 100, "\n")

# Stop writing to file
sink()
