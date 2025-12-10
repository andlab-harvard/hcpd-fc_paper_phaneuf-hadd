##### Define Shared Variables and Functions for HCPD FC #####

# Written by: Camille Phaneuf-Hadd (cphaneuf@g.harvard.edu)
# Last updated: 12/6/25

##### Plotting Utilities ##### 

# Set plotting color scheme
flanker_blue <- "#0C2E73"
dccs_blue <- "#1A6591"
pattern_blue <- "#8BBCD0"
psmt_green <- "#2C6159"
list_green <- "#57A384"
comp_yellow <- "#C3A52F"
reg_gray <- "#595959"
light_gray <- "#AFABAB"
reg_brown <- "#945200"
age_8_10_orange <- "#F96302"
age_11_13_red <- "#DC2004"
age_14_17_pink <- "#F11D80"
age_18_21_purple <- "#9C2F71"

# Set common color list
cols <- c(flanker_blue, dccs_blue, pattern_blue, psmt_green, list_green)

# Set plotting theme
plot_theme <- theme(title = element_text(size = 24, face = "bold", family = "Avenir"),
                    plot.title = element_text(hjust = .5),
                    axis.title.x = element_text(size = 24, family = "Avenir"),
                    axis.title.y = element_text(size = 24, family = "Avenir"),
                    axis.text.x = element_text(size = 16, colour = "black", family = "Avenir"),
                    axis.text.y = element_text(size = 16, colour = "black", family = "Avenir"),
                    legend.text = element_text(size = 16, colour = "black", family = "Avenir"),
                    legend.position = "right",
                    legend.key = element_rect(fill = "transparent", color = NA),
                    strip.text.x = element_text(size = 16, colour = "black", family = "Avenir"),
                    strip.text.y = element_text(size = 16, colour = "black", family = "Avenir"),
                    panel.grid.major = element_blank(), # Remove grid marks
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    axis.line = element_line(colour = "black"))

# Set other plotting variables
node_fontsize <- 1

##### Task Utilities ##### 

# Set common data sub-setting and/or labeling variables
short_task_labels <- c("Flnk", "DCCS", "PnCn", "PSMT", "LtSg")
long_task_labels <- c("Flanker", "DCCS", "PatternComparison", "PSMT", "ListSorting")
toolbox_z <- c("Flanker_z", "DCCS_z", "Pattern_z", "PSMT_z", "List_z")
no_flanker_labels <- c("Flanker", "Remaining Tasks")
no_dccs_labels <- c("DCCS", "Remaining Tasks")
no_pattern_labels <- c("Pattern Comparison", "Remaining Tasks")
no_psmt_labels <- c("PSMT", "Remaining Tasks")
no_list_labels <- c("List Sorting", "Remaining Tasks")

##### GAM Utilities ##### 

# Credit: Katherine Grisanzio 
# (code https://osf.io/fvy8d for paper https://andl.wjh.harvard.edu/files/2023/11/Grisanzio-et-al-2023-1.pdf)

# Define function for dealing with float variable type
too_small <- function(x) abs(x) < 10^-15

# Define function for notating significance of first derivative, computed using gratia::derivatives
clip_on_siggratia <- function(ci){
  # Does confidence interval include 0? Confidence interval includes 0 when:
  #   - lower interval is negative and upper interval is positive (or vice versa)
  #   - lower interval and upper interval are very close to 0
  not_sig <- (ci$.lower_ci * ci$.upper_ci < 0) | (too_small(ci$.lower_ci) & too_small(ci$.upper_ci))
  ci$sig <- 1 # Mark all points as significant regions of age-related change
  ci$sig[not_sig] <- 0 # Fix points that are actually non-significant regions of age-related change
  return(ci)
}
