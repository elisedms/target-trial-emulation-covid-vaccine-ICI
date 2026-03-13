##################################
# Libraries and utils
##################################

###############
#Libraries
library(openxlsx)
library(tidyverse)
library(survival)
library(survminer)
library(cowplot)

###############
#Set seed
set.seed(28)

###############
#Global parameters
label_order = c("Vaccine", "No vaccine")
#Blue for vaccine - red for no vaccine
color_for_plots = c("#00649A", "#c1121f")
#Number of repetitions of bootstrap
number_bootstrap = 1000

###############
#Utils functions
source("utils_function.R")

###############
#Survival analyses

# NSCLC patients : Original analyses
same_time_zero = FALSE
remove_historic_control = FALSE
title_txt = "Initial analyses\n"
source("nsclc_original.R")
plot_survival_curvesA = plot_survival_curves_nsclc_initial

# NSCLC patients : Patients enter at first ICI treatment line (and first systemic treatment) only
same_time_zero = TRUE
remove_historic_control = FALSE
title_txt = "After alignment of time zero\nat first ICI treatment"
source("nsclc_original.R")
plot_survival_curvesB = plot_survival_curves_nsclc_initial

# NSCLC patients : Patients enter at first ICI treatment line only + removal of historical controls
same_time_zero = TRUE
remove_historic_control = TRUE
title_txt = "After exclusion of patients\ndiagnosed before 2020"
source("nsclc_original.R")
plot_survival_curvesC = plot_survival_curves_nsclc_initial

# NSCLC patients : Target trial emulation (cloning, censoring and weighting)
source("nsclc_target_trial.R")
plot_survival_curvesD = plot_survival_curves_nsclc

###############
#Final plot assembly
plot_tot = cowplot::plot_grid(plot_survival_curvesA,
                              plot_survival_curvesB,
                              plot_survival_curvesC,
                              plot_survival_curvesD,
                   labels = c("A","B","C","D"), nrow = 2)

#Save plot
pdf("Fig_survival_curves.pdf",
    width = unit(6.6,"cm"),
    height = unit(8, "cm"))
  plot_tot
dev.off()