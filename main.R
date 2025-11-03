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

#Target trial emulation
source("nsclc.R") #For NSCLC patients
source("melanoma.R") #For melanoma patients - outcome overall survival
source("melanoma_pfs.R") #For melanoma patients - outcome progression free survival

###############
#Final plot assembly
plot_tot = cowplot::plot_grid(plot_survival_curves_nsclc, plot_survival_curves_melanoma, plot_survival_curves_melanoma_pfs,
                   labels = c("A","B","C"), nrow = 1)
plot_tot

#Save plot
pdf("Fig_survival_curves.pdf",
    width = unit(10,"cm"),
    height = unit(4, "cm"))
  plot_tot
dev.off()