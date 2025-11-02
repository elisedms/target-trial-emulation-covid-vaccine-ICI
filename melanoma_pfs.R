#########################################
# TARGET TRIAL EMULATION : MELANOMA PATIENTS
# OUTCOME : Progression free survival
#######################################

##################################
# Renaming and formatting of variables
##################################

#Load data 
#PFS is a composite outcome -- that includes death
df =  read.xlsx(xlsxFile = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09655-y/MediaObjects/41586_2025_9655_MOESM3_ESM.xlsx",
          sheet = "1e Raw Data")

#For simplicity to run the rest of the code -- we rename PFS variables as OS
df = df %>%
  #Renaming variables with special characters in name
  #For the sake of simplicity, we rename as status_os and delay_os the PFS variables (easier for the codes to run afterwards)
  rename(status_os = `PFS_Code`,
         delay_os = `PFS_`,
         group = `Group_1=vaccine`
  ) %>%
  #Age as numeric
  #Three patients had age written as >89 ; we impute by 90 years as written in the manuscript
  mutate(Age_at_ICI_start = ifelse(Age_at_ICI_start == ">89", "90", Age_at_ICI_start)) %>%
  mutate(Age_at_ICI_start = as.numeric(Age_at_ICI_start)) %>%
  #Relabel vaccine group
  mutate(group = ifelse(group == 1, "Observed vaccine", "Observed no vaccine"))

##################################
# Restriction to eligible patients
##################################
df_restricted = df %>% 
  #During the years when the vaccine was available
  #Vaccine available December 2020 -- with the grace period of 100 days we want to include patients starting September 2020
  #Given that we only have the year of ICI start available we have to include starting either 2020 or 2021.
  #We chose to keep 2020 to increase the number of patients.
  #This will be discussed on the target trial emulation protocol.
  filter(Treatment_Year_ICI >= 2020) %>%
  #Add a patient ID
  rownames_to_column("row_number")

print(dim(df_restricted))
summary(factor(df_restricted$group))

##################################
# Cloning and censoring
##################################

df_restricted_clone = clone_censor(df_restricted)
#Check
summary(factor(df_restricted_clone$clone))

##################################
# Weighting
##################################

#In theory weights should be updated every month during the grace period -- accounting for time-varying confounding covariates
#But time-varying confounding factors were not available in the dataset
#So we adjusted on baseline confounding factors only and assume no time-varying confounding factors

#Variables adjusted for : we adjust on patient socio-demographic characteristics, comorbid conditions, and the year of treatment

variables_used_for_adjustment = c("Gender",
                                  "Age_at_ICI_start",
                                  "ECOG", 
                                  "Had_Heart_Disease",
                                  "Had_Immune_Deficiency",
                                  "Had_Diabetes",
                                  "Had_Respiratory_Condition",
                                  "Previous_history_of_malignancy_at_ICI_start",
                                  "English_as_primary_language",
                                  "Had_CKD",
                                  "CNS_disease")
weighting = weight(df_restricted_clone,
                   variables_used_for_adjustment,
                   print = TRUE, #To get plots for quality check
                   color_for_plots)
df_restricted_clone_long = weighting[[1]] #Dataset with weight - long format
df_restricted_clone = weighting[[2]] #Dataset with weight - short format
#Plots for quality check
#Propensity score overlap
print(weighting[[3]])
#Distribution of the weights
print(weighting[[4]])

#####################################################
# Estimate survival curves
####################################################

#Set order
df_restricted_clone_long$clone = factor(df_restricted_clone_long$clone,
                                        levels = label_order)
df_restricted_clone$clone = factor(df_restricted_clone$clone,
                                   levels = label_order)

surv = estimate_survival(df_restricted_clone_long, 
                         df_restricted_clone,
                         color_for_plots,
                         ylab_txt = "PFS (%)",
                         print = TRUE,
                         time_vector = c(12,24,36),
                         title_txt = "Stage IV melanoma - PFS")

plot_survival_curves = surv[[1]]
plot_survival_tables = surv[[2]]
surv_point_estimates = surv[[3]]

###########################################################
# We use non-parametric bootstrapping with 1000 repetitions
# for the confidence intervals
##########################################################

#Non-parametric bootsrapping
sample_size = nrow(df_restricted)

#Vector to store the boostrap estimates
surv_estimates_bootstrap = c()

for (bootstrap_number in 1:number_bootstrap){
  
  # Create bootsrap dataset
  id_sample <- sample(1:sample_size, size = sample_size, replace = TRUE)
  df_bootstrap <- suppressMessages(df_restricted %>%
    right_join(data.frame(row_number = as.character(id_sample))) %>%
    as.data.frame())
  
  #Clone and censor
  df_bootstrap_clone = clone_censor(df_bootstrap)
  #Get weights
  df_bootstrap_clone_list = weight(df_bootstrap_clone,
                                   variables_used_for_adjustment,
                                   print = FALSE)
  df_bootstrap_clone_long = df_bootstrap_clone_list[[1]]
  df_bootstrap_clone = df_bootstrap_clone_list[[2]]
  #Estimate Kaplan-Meier curves and 3-year OS
  df_bootstrap_clone_long$clone = factor(df_bootstrap_clone_long$clone,
                                         levels = label_order)
  df_bootstrap_clone$clone = factor(df_bootstrap_clone$clone,
                                    levels = label_order)
  surv = estimate_survival(df_bootstrap_clone_long,
                           df_bootstrap_clone,
                           print = FALSE,
                           time_vector = c(12,24,36))

  surv_estimates_bootstrap = bind_rows(surv_estimates_bootstrap, surv[[1]] %>% mutate(rep_boot = bootstrap_number))
  
}

#Estimate 95% confidence intervals (from percentiles)
#And get labels at three years 
list_ci = estimate_ci(surv_point_estimates, surv_estimates_bootstrap,
                      outcome = "PFS")

label_ci_all = list_ci[[1]]
label_ci_diff = list_ci[[2]]
label_tot = list_ci[[3]]

##Print results
print("Stage IV melanoma - PFS")
print(label_ci_all)
print(paste0("Difference : \n", label_ci_diff))

###########################################################
# Final tuning of the survival curves plot
##########################################################

#Add the label for 3-year PFS
plot_survival_curves = plot_survival_curves +
  annotate(geom = "text", x = 1, y = 0.15, label = label_tot, hjust = 0, size = 3)

#Add error bars at 12, 24, 36 months
plot_survival_curves = plot_survival_curves +
  geom_errorbar(data = label_ci_all, inherit.aes = FALSE,
                aes(x = time,
                    ymin = CI_low, ymax = CI_high,
                    color = clone, group = clone),
                position=position_dodge(width=0.5),
                width = 1)

#Trick to align the risk table and the curves vertically
plot_survival_curves_align <- cowplot::align_plots(plot_survival_curves,
                                                   plot_survival_tables + theme(plot.title = element_text(size=10)) +
                                                     theme(legend.position = "none"),
                                                   align = 'v', axis ='l')
plot_survival_curves_melanoma_pfs <- cowplot::plot_grid(plot_survival_curves_align[[1]],
                                                    plot_survival_curves_align[[2]],
                                                    ncol = 1, rel_heights = c(6,2))

print(plot_survival_curves_melanoma_pfs)
