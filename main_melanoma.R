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
#Overall survival
###############

# Melanoma patients : Original analyses
remove_historic_control = FALSE
title_txt = "Initial analyses\n"
outcome = "OS"
source("melanoma_original.R")
plot_survival_curvesA = plot_survival_curves_melanoma_initial

#Melanoma patients : already only at first ICI treatment line

#Melanoma patients : Removal of historical control
remove_historic_control = TRUE
outcome = "OS"
title_txt = "After exclusion of patients\ndiagnosed before 2020"
source("melanoma_original.R")
plot_survival_curvesC = plot_survival_curves_melanoma_initial

# Melanoma patients : Target trial emulation
source("melanoma_target_trial.R")
plot_survival_curvesD = plot_survival_curves_melanoma

###############
#Progression-free survival
###############

# Melanoma patients : Original analyses
remove_historic_control = FALSE
title_txt = "Initial analyses\n"
outcome = "PFS"
source("melanoma_original.R")
plot_survival_curvesA_pfs = plot_survival_curves_melanoma_initial

#Melanoma patients : already only at first ICI treatment line

#Melanoma patients : Remove of historical control (year starts 2020)
remove_historic_control = TRUE
outcome = "PFS"
title_txt = "After exclusion of patients\ndiagnosed before 2020"
source("melanoma_original.R")
plot_survival_curvesC_pfs = plot_survival_curves_melanoma_initial

# Melanoma patients : Target trial emulation
outcome = "PFS"
source("melanoma_target_trial.R")
plot_survival_curvesD_pfs = plot_survival_curves_melanoma

###############
#Final plot assembly
###############

plot_tot = cowplot::plot_grid(plot_survival_curvesA,
                              plot_survival_curvesC,
                              plot_survival_curvesD,
                              plot_survival_curvesA_pfs,
                              plot_survival_curvesC_pfs,
                              plot_survival_curvesD_pfs,
                   labels = c("A","B","C", "D", "E","F"), nrow = 2)
plot_tot

#Save plot
pdf("Fig_survival_curves_melanoma.pdf",
    width = unit(1.5*6.6,"cm"),
    height = unit(8, "cm"))
  plot_tot
dev.off()